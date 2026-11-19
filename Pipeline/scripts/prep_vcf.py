import argparse
import os
import re
import shutil
import sys
import tempfile
import numpy as np
import pandas as pd
import pysam

# Add local directory to path for importing pipeline utilities
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pgs_harmonise import (clean_chrom, is_palindromic, match_variants, resolve_palindromic, reverse_complement)
from prep_vcf_safety import (
    exact_effect_alt_index, first_nonmissing, is_concrete_allele, is_reference_block,
    recode_gt_to_effect, reference_block_status, select_target_records,
    reference_block_effect_role, safe_info_get,
    infer_ref_for_effect_alt_from_fasta, match_effect_only_concrete_record,
    inspect_vcf_genotype_convention, infer_reference_allele_pair_from_fasta,
)

# Standard quality control flags to be appended to notes when a variant fails QC
QC_NOTES = ("LowRead", "LowQual", "IncorrectRatio")


# ==============================================================================
# 1. METADATA & CONFIGURATION LOADING
# ==============================================================================

def load_variant_actions(actions_file):
    """
    Loads manual variant overrides (Keep, Drop, Force Impute) 
    configured by the user via the Streamlit UI.
    """
    actions_dict = {}
    if actions_file and os.path.exists(actions_file) and os.path.getsize(actions_file) > 0:
        try:
            df = pd.read_csv(actions_file, sep='\t')
            if 'Action' in df.columns:
                chr_col = next((c for c in df.columns if c.lower() in ['chr', 'chr_name', 'hm_chr']), None)
                pos_col = next((c for c in df.columns if c.lower() in ['pos', 'hm_pos', 'chr_position']), None)
                if chr_col and pos_col:
                    for _, row in df.iterrows():
                        c = clean_chrom(row[chr_col])
                        p = str(row[pos_col]).strip()
                        actions_dict[f"{c}:{p}"] = str(row['Action']).strip()
        except Exception: 
            pass
    return actions_dict

def get_harmonized_known_ref(row, *, is_hg19, chrom, pos):
    """
    Extracts the strict matched reference allele from the harmonization audit.
    Applies strand flipping if the original match was on the opposite strand.
    """
    try:
        # Ignore weak or ambiguous matches (match_strength != 0)
        if int(float(row.get('match_strength', ''))) != 0:
            return None
    except (TypeError, ValueError):
        return None

    source = str(row.get('source', '')).strip().lower()
    matched_pos = str(row.get('matched_pos', '')).strip().replace('.0', '')
    matched_chrom = clean_chrom(row.get('matched_chrom', ''))
    
    if matched_chrom != clean_chrom(chrom):
        return None
        
    if source == 'primary':
        if matched_pos != str(pos) or matched_chrom != clean_chrom(chrom):
            return None
    elif source != 'fallback':
        return None
        
    known_ref = str(row.get('matched_ref', '')).strip().upper()
    if not known_ref:
        return None
        
    flipped = str(row.get('flipped', '')).strip().lower() in {'true', '1', 'yes'}
    if flipped:
        known_ref = reverse_complement(known_ref)
        
    return known_ref or None

def get_somalier_ancestry(master_somalier_file, sample_id):
    """Retrieves the predicted ancestry for the current sample to pick the correct Allele Frequency."""
    df = pd.read_table(master_somalier_file)
    id_col = '#sample_id' if '#sample_id' in df.columns else 'sample_id'
    patient_row = df[df[id_col] == sample_id]
    
    if patient_row.empty: 
        return "EUR"
        
    if 'predicted_ancestry' in df.columns: 
        return str(patient_row['predicted_ancestry'].iloc[0]).strip()
        
    ancestry_cols = [c for c in df.columns if c.endswith('_prob')]
    if ancestry_cols: 
        return patient_row[ancestry_cols].iloc[0].idxmax().split('_')[0]
        
    return "EUR"


# ==============================================================================
# 2. POPULATION ALLELE FREQUENCY HANDLING
# ==============================================================================

def load_cached_frequencies(cache_file_path, ancestry, build_query_col):
    """
    Loads precomputed Allele Frequencies (AF) from the cached reference distributions
    for the sample's specific ancestry.
    """
    if not os.path.exists(cache_file_path) or cache_file_path == "NONE": 
        return {}
        
    cached_df = pd.read_csv(cache_file_path, sep=None, engine='python')
    if cached_df.empty: 
        return {}
        
    # Attempt to locate coordinate columns
    chr_col = next((c for c in cached_df.columns if c.lower() in ['chr', '#chrom', 'chrom']), None)
    ref_col = next((c for c in cached_df.columns if c.lower() == 'ref'), None)
    alt_col = next((c for c in cached_df.columns if c.lower() == 'alt'), None)
    pos_col = next((c for c in cached_df.columns if c.lower() == 'pos'), None)
    id_col = next((c for c in cached_df.columns if c.lower() == 'id'), None)
    
    # Locate the appropriate AF column for the patient's ancestry
    target_af = f"{ancestry}_af".lower()
    af_col = next((c for c in cached_df.columns if c.lower() == target_af), None)
    if not af_col: 
        af_col = next((c for c in cached_df.columns if c.lower() == "eur_af"), None)
    if not af_col:
        non_af = {chr_col, pos_col, ref_col, alt_col, id_col, 'obs_ct'}
        af_cols = [c for c in cached_df.columns if c not in non_af and c is not None]
        af_col = af_cols[-1] if af_cols else cached_df.columns[-1]

    af_dict = {}
    
    # Index standard VCF-like formats
    if chr_col and ref_col and alt_col and (pos_col or id_col):
        for _, row in cached_df.iterrows():
            c = str(row[chr_col]).strip()
            if c.lower().startswith('chr'): c = c[3:]
            if c.endswith('.0'): c = c[:-2]
            
            p = ""
            if pos_col and pd.notna(row[pos_col]): p = str(row[pos_col]).strip()
            elif id_col and pd.notna(row[id_col]): p = str(row[id_col]).strip()
            if p.endswith('.0'): p = p[:-2]
            if ':' in p: p = p.split(':')[1]
            
            key = f"chr{c}:{p}"
            r, a = str(row[ref_col]).strip().upper(), str(row[alt_col]).strip().upper()
            
            try: 
                val = float(row[af_col])
            except (ValueError, TypeError): 
                val = 0.0
                
            if key not in af_dict: af_dict[key] = []
            af_dict[key].append((r, a, val))
        return af_dict

    # Index custom query-based formats
    if build_query_col in cached_df.columns: query_col = build_query_col
    elif "Original_Query" in cached_df.columns: query_col = "Original_Query"
    else: query_col = cached_df.columns[0]
    
    for _, row in cached_df.iterrows():
        try: val = float(row[af_col])
        except: val = 0.0
        
        q = str(row[query_col])
        parts = q.split('_')
        if len(parts) >= 2:
            key, ra = parts[0], parts[1].split('>')
            if len(ra) == 2:
                if key not in af_dict: af_dict[key] = []
                af_dict[key].append((ra[0], ra[1], val))
    return af_dict

def get_robust_af(row, af_dict, assume_forward, effect_is_alt):
    """
    Safely retrieves the population Allele Frequency for a target variant,
    accounting for potential strand flips and ambiguous multi-allelic sites.
    """
    chrom = str(row['hm_chr_norm'])
    pos = str(row['active_pos'])
    eff = str(row.get('effect_allele', '')).strip().upper()
    oth = str(row.get('clean_other', 'N')).strip().upper()
    if oth in ('', 'N', '.'): oth = ''
    
    key = f"chr{chrom}:{pos}"
    records = af_dict.get(key, [])
    if not records: 
        return 0.0

    if effect_is_alt:
        # Require strict REF matching when effect is explicitly forced to be ALT
        for r, a, af in records:
            ref = str(r).strip().upper()
            alts = [x.strip().upper() for x in str(a).split(',') if x]
            if oth and ref != oth:
                continue
            if eff in alts:
                return af
        return 0.0

    # Match across all alleles for normal harmonization
    for r, a, af in records:
        alts = [x for x in a.split(',') if x]
        allele_set = [r] + alts
        
        matches, ambiguous, relaxed = match_variants(
            eff, oth, allele_set, allow_complement=not assume_forward
        )
        ambiguous = len(matches) > 1
        
        if matches and not relaxed and not ambiguous and len(matches) == 1:
            idx, _flipped = matches[0]
            target = allele_set[idx]
            
            if target in alts: return af
            if target == r: return 1.0 - af if af > 0 else 0.0

    return 0.0


# ==============================================================================
# 3. VARIANT QC & CONFIDENCE CHECKING
# ==============================================================================

def check_lcr(prs_df, lcr_bed, pos_col='active_pos'):
    """Flags variants that fall into known Low Complexity Regions (LCR) for quality tracking."""
    lcr_flags = []
    try:
        lcr = pysam.TabixFile(lcr_bed)
        lcr_contigs = set(lcr.contigs)
    except Exception: 
        return [False] * len(prs_df)
        
    for _, row in prs_df.iterrows():
        norm_c, pos = str(row['hm_chr_norm']), int(row[pos_col])
        c1, c2 = f"chr{norm_c}", norm_c
        chrom = c1 if c1 in lcr_contigs else (c2 if c2 in lcr_contigs else c1)
        
        try: 
            lcr_flags.append(len(list(lcr.fetch(chrom, pos - 1, pos))) > 0)
        except ValueError: 
            lcr_flags.append(False)
            
    return lcr_flags

def _infer_depths_from_dp(target_count, total_dp):
    """Estimates allele depths for homozygous and perfectly balanced heterozygous calls."""
    if target_count == 0: return total_dp, 0
    elif target_count == 2: return 0, total_dp
    elif target_count == 1: return total_dp // 2, total_dp - (total_dp // 2)
    return None, None

def get_ref_alt_depths(variant, sample_id, match_vcf_idx):
    """
    Safely extracts Reference and Alternate allele read depths from various VCF formats.
    Checks multiple standard fields (AD, ADALL, TAR/TIR, RO/AO, DP4).
    """
    sample_data = variant.samples[sample_id]
    eff_idx = match_vcf_idx if match_vcf_idx > 0 else 1
    alt_array_idx = eff_idx - 1
    
    def parse_ad(ad_field):
        if ad_field is None: return None
        if isinstance(ad_field, str): ad_field = ad_field.split(',')
        if isinstance(ad_field, (list, tuple)) and len(ad_field) > eff_idx and ad_field[0] is not None and ad_field[eff_idx] is not None:
            try: return int(ad_field[0]), int(ad_field[eff_idx])
            except (ValueError, TypeError): return None
        return None
        
    if 'ADALL' in sample_data:
        res = parse_ad(sample_data['ADALL'])
        if res: return res
    if 'AD' in sample_data:
        res = parse_ad(sample_data['AD'])
        if res: return res
    if 'TAR' in sample_data and 'TIR' in sample_data:
        try: return int(sample_data['TAR'][0]) if sample_data['TAR'] else 0, int(sample_data['TIR'][alt_array_idx]) if sample_data['TIR'] else 0
        except (ValueError, TypeError, IndexError): pass
    if 'RO' in sample_data and 'AO' in sample_data:
        try: return int(sample_data['RO']), int(sample_data['AO'][alt_array_idx] if isinstance(sample_data['AO'], (tuple, list)) else sample_data['AO'])
        except (ValueError, TypeError, IndexError): pass
    
    dp4 = safe_info_get(variant, 'DP4')
    if dp4 is not None and len(dp4) == 4:
        try: return int(dp4[0]) + int(dp4[1]), int(dp4[2]) + int(dp4[3])
        except (ValueError, TypeError): pass

    # If concrete depths are missing, infer from DP and Genotype
    target_count = sample_data.get('GT', []).count(match_vcf_idx) if match_vcf_idx > 0 and sample_data.get('GT') else 0
    dp_val = first_nonmissing(
        sample_data.get('DP'), sample_data.get('MIN_DP'),
        safe_info_get(variant, 'DP'), safe_info_get(variant, 'MinDP'),
    )
    if dp_val is not None:
        try:
            return _infer_depths_from_dp(target_count, int(dp_val))
        except (ValueError, TypeError):
            pass
            
    return None, None

def variant_qa_checker(variant, sample_id, qa_idx, notes_list, depth=15, min_gq=50.0, min_af=0.20):
    """
    Validates a concrete variant against QC thresholds (Minimum Depth, Genotype Quality, Allele Fraction).
    Appends failure reasons to the notes_list.
    """
    if depth <= 0: return True
    if sample_id not in variant.samples: return False
    
    sample_data = variant.samples[sample_id]
    GT = sample_data.get('GT')
    
    # Check for empty reference blocks
    is_ref_block = False
    if variant.alts is None: is_ref_block = True
    elif variant.alts: is_ref_block = any(a in ('<*>', '<NON_REF>') for a in variant.alts)
    
    if (GT is None or None in GT) and not is_ref_block: 
        return False

    ref_dp, alt_dp = get_ref_alt_depths(variant, sample_id, qa_idx)
    passed = True
    
    # 1. Depth Check
    dp_val = first_nonmissing(
        sample_data.get('DP'), sample_data.get('MIN_DP'),
        safe_info_get(variant, 'DP'), safe_info_get(variant, 'MinDP'),
    )
    if dp_val is None and ref_dp is not None and alt_dp is not None:
        dp_val = ref_dp + alt_dp
        
    if dp_val is not None:
        try:
            if int(dp_val) < depth:
                if "LowRead" not in notes_list: notes_list.append("LowRead")
                passed = False
        except (ValueError, TypeError): pass

    # 2. Genotype Quality Check (GQ, or inferred from GL/PL)
    gq_val = first_nonmissing(sample_data.get('GQ'), sample_data.get('RGQ'))
    if gq_val is None and 'GL' in sample_data:
        gls = sample_data['GL']
        if gls and len(gls) > 1:
            sorted_gls = sorted([g for g in gls if g is not None], reverse=True)
            if len(sorted_gls) > 1: gq_val = abs(sorted_gls[0] - sorted_gls[1]) * 10.0
    elif gq_val is None and 'PL' in sample_data:
        pls = sample_data['PL']
        if pls and len(pls) > 1:
            sorted_pls = sorted([p for p in pls if p is not None])
            if len(sorted_pls) > 1: gq_val = abs(sorted_pls[1] - sorted_pls[0])

    if gq_val is not None:
        if isinstance(gq_val, (tuple, list)): gq_val = gq_val[0]
        try:
            if float(gq_val) < min_gq:
                if "LowQual" not in notes_list: notes_list.append("LowQual")
                passed = False
        except (ValueError, TypeError): pass

    # 3. Allele Fraction Ratio Check
    if ref_dp is not None and alt_dp is not None and (ref_dp + alt_dp) > 0:
        af = alt_dp / (ref_dp + alt_dp)
        max_af = 1.0 - min_af
        target_count = GT.count(qa_idx) if GT else 0
        
        if qa_idx > 0 and GT is not None and None not in GT:
            if target_count == 0 and af > min_af:
                if "IncorrectRatio" not in notes_list: notes_list.append("IncorrectRatio")
                passed = False
            elif target_count == 1 and (af < min_af or af > max_af):
                if "IncorrectRatio" not in notes_list: notes_list.append("IncorrectRatio")
                passed = False
            elif target_count == 2 and af < max_af:
                if "IncorrectRatio" not in notes_list: notes_list.append("IncorrectRatio")
                passed = False
                
    return passed


# ==============================================================================
# 4. MAIN PIPELINE EXECUTION
# ==============================================================================

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--build', type=str, default='GRCh38')
    parser.add_argument('--ref_fasta', type=str, default='', help='Reference FASTA for reference-block and absent-standard-VCF allele resolution')
    parser.add_argument('--input_is_gvcf', type=str, default='False', help='True when the original source input is a gVCF')
    parser.add_argument('--vcf', required=True)
    parser.add_argument('--prs', required=True)
    parser.add_argument('--somalier', required=True)
    parser.add_argument('--lcr', required=True)
    parser.add_argument('--af_cache_file', required=True)
    parser.add_argument('--af_mode', type=str, default="1000G")
    parser.add_argument('--audit_file', type=str, default="")
    parser.add_argument('--sample', required=True)
    parser.add_argument('--out_vcf', required=True)
    parser.add_argument('--out_chroms', required=True)
    parser.add_argument('--min_depth', type=int, default=15)
    parser.add_argument('--min_gq', type=float, default=50.0)
    parser.add_argument('--min_af', type=float, default=0.20)
    parser.add_argument('--variant_actions', type=str, default="")
    parser.add_argument('--assume_forward', type=str, default="False")
    parser.add_argument('--effect_is_alt', type=str, default="False")
    parser.add_argument('--skip_beagle', type=str, default="False")
    parser.add_argument('--pal_low', type=float, default=0.40)
    parser.add_argument('--pal_high', type=float, default=0.60)
    args = parser.parse_args()

    # Apply configuration toggles
    assume_forward = args.assume_forward.lower() == "true"
    effect_is_alt = args.effect_is_alt.lower() == "true"
    skip_beagle = args.skip_beagle.lower() == "true"
    input_is_gvcf = str(args.input_is_gvcf).strip().lower() in {"true", "1", "yes"}
    is_hg19 = str(args.build).strip().lower() in {'hg19', 'grch37', '37'}
    
    pos_col = 'chr_position' if is_hg19 else 'hm_pos'
    query_id_col = 'query_id_hg19' if is_hg19 else 'query_id_hg38'

    actions_dict = load_variant_actions(args.variant_actions)
    harmonization_audit_df = None
    
    # Load and validate the prior harmonization audit logs
    if args.audit_file and os.path.exists(args.audit_file):
        try:
            audit_df = pd.read_csv(args.audit_file, sep='\t', dtype=str)
            harmonization_audit_df = audit_df.copy()
            if (audit_df['status'].astype(str).str.lower() == 'kept').any():
                required_audit = {'matched_ref','matched_alts','matched_build','match_strength','matched_pos','matched_chrom'}
                missing_audit = required_audit - set(audit_df.columns)
                if missing_audit:
                    raise ValueError(
                        "Harmonization audit is from an older format and lacks "
                        f"{sorted(missing_audit)}. Rerun build_reference_cache.py before prep_vcf.py."
                    )
                    
            # Propagate Dropped variants into the actions dictionary
            dropped_df = audit_df[audit_df['status'] == 'Dropped']
            for _, row in dropped_df.iterrows():
                chrom = clean_chrom(row['chr'])
                pos = int(row['pos_hg19']) if is_hg19 else int(row['pos_hg38'])
                action_key = f"{chrom}:{pos}"
                if action_key not in actions_dict:
                    actions_dict[action_key] = "Drop"
        except Exception as e:
            if "older format" in str(e): raise
            print(f"Warning: Failed to load audit file: {e}")

    # Load and normalize the target PGS panel
    prs_df = pd.read_csv(args.prs, sep='\t', comment='#')
    
    if 'hm_chr' not in prs_df.columns and 'chr' in prs_df.columns: prs_df['hm_chr'] = prs_df['chr']
    if pos_col not in prs_df.columns and 'pos' in prs_df.columns: prs_df[pos_col] = prs_df['pos']
    if 'hm_chr' not in prs_df.columns: sys.exit(1)

    prs_df['hm_chr_norm'] = prs_df['hm_chr'].apply(clean_chrom)
    prs_df = prs_df[prs_df['hm_chr_norm'] != ""]
    prs_df['active_pos'] = pd.to_numeric(prs_df[pos_col], errors='coerce').fillna(0).astype(int)
    prs_df = prs_df[prs_df['active_pos'] > 0]

    # Extract target allele frequencies
    if 'effect_allele_frequency' in prs_df.columns: prs_df['eaf_clean'] = pd.to_numeric(prs_df['effect_allele_frequency'], errors='coerce')
    elif 'hm_effect_allele_frequency' in prs_df.columns: prs_df['eaf_clean'] = pd.to_numeric(prs_df['hm_effect_allele_frequency'], errors='coerce')
    elif 'allelefrequency_effect' in prs_df.columns: prs_df['eaf_clean'] = pd.to_numeric(prs_df['allelefrequency_effect'], errors='coerce')
    else: prs_df['eaf_clean'] = np.nan

    # Normalize Effect and Other alleles based on genome build conventions
    if is_hg19:
        other_source = 'other_allele' if 'other_allele' in prs_df.columns else 'hm_other_allele'
        effect_source = 'effect_allele' if 'effect_allele' in prs_df.columns else 'alt'
    else:
        other_source = 'hm_inferOtherAllele' if 'hm_inferOtherAllele' in prs_df.columns else ('hm_other_allele' if 'hm_other_allele' in prs_df.columns else 'other_allele')
        effect_source = 'hm_effect_allele' if 'hm_effect_allele' in prs_df.columns else ('effect_allele' if 'effect_allele' in prs_df.columns else 'alt')
        
    if other_source in prs_df.columns:
        prs_df['clean_other'] = prs_df[other_source].fillna('').astype(str).str.split('/').str[0].str.upper()
        if not is_hg19:
            if 'hm_other_allele' in prs_df.columns and other_source != 'hm_other_allele':
                hm_other = prs_df['hm_other_allele'].fillna('').astype(str).str.upper()
                prs_df['clean_other'] = prs_df['clean_other'].where(prs_df['clean_other'].str.strip() != '', hm_other)
            if 'other_allele' in prs_df.columns:
                original_other = prs_df['other_allele'].fillna('').astype(str).str.upper()
                prs_df['clean_other'] = prs_df['clean_other'].where(prs_df['clean_other'].str.strip() != '', original_other)
    elif 'ref' in prs_df.columns:
        prs_df['clean_other'] = prs_df['ref'].fillna('').astype(str).str.upper()
    else:
        prs_df['clean_other'] = ''
        
    if effect_source in prs_df.columns:
        selected_effect = prs_df[effect_source].fillna('').astype(str).str.upper()
        if (not is_hg19) and 'effect_allele' in prs_df.columns and effect_source == 'hm_effect_allele':
            original_effect = prs_df['effect_allele'].fillna('').astype(str).str.upper()
            selected_effect = selected_effect.where(selected_effect.str.strip() != '', original_effect)
        prs_df['effect_allele'] = selected_effect
    else:
        raise ValueError(f"No effect allele column for build {args.build}")

    # Transfer harmonization audit data directly to the active variant dataframe
    for col in ('matched_ref','matched_alts','matched_build','match_strength','matched_pos','matched_chrom','matched_effect_allele','source','flipped'):
        prs_df[col] = ''
        
    if harmonization_audit_df is not None and not harmonization_audit_df.empty:
        audit = harmonization_audit_df.copy()
        pos_audit_col = 'pos_hg19' if is_hg19 else 'pos_hg38'
        required = {'chr', pos_audit_col, 'status'}
        if required.issubset(audit.columns):
            audit = audit[audit['status'].astype(str).str.lower() == 'kept'].copy()
            audit['_key'] = audit.apply(
                lambda r: f"{clean_chrom(r.get('chr', ''))}:{str(r.get(pos_audit_col, '')).replace('.0', '')}:"
                          f"{str(r.get('target_effect_allele', '')).strip().upper()}:"
                          f"{str(r.get('target_other_allele', '')).strip().upper()}", axis=1
            )
            audit = audit.drop_duplicates('_key', keep='first').set_index('_key')
            
            prs_df['_key'] = prs_df.apply(
                lambda r: f"{r['hm_chr_norm']}:{int(r['active_pos'])}:"
                          f"{str(r['effect_allele']).strip().upper()}:"
                          f"{str(r['clean_other']).strip().upper()}", axis=1
            )
            
            map_cols = {'matched_ref':'matched_ref','matched_alts':'matched_alts','matched_build':'matched_build','match_strength':'match_strength','matched_pos':'matched_pos','matched_chrom':'matched_chrom','matched_effect_allele':'matched_effect_allele','source':'source','flipped':'flipped'}
            for dest, src in map_cols.items():
                if src in audit.columns:
                    prs_df[dest] = prs_df['_key'].map(audit[src]).fillna('')
            prs_df.drop(columns=['_key'], inplace=True)

    # Establish Pop AF
    ancestry = get_somalier_ancestry(args.somalier, args.sample)
    gnomad_afs = load_cached_frequencies(args.af_cache_file, ancestry, build_query_col=query_id_col)

    if args.af_mode == "PGS_FILE":
        if 'allelefrequency_effect' in prs_df.columns: prs_df['Pop_AF'] = pd.to_numeric(prs_df['allelefrequency_effect'], errors='coerce').fillna(0.0)
        elif 'hm_effect_allele_frequency' in prs_df.columns: prs_df['Pop_AF'] = pd.to_numeric(prs_df['hm_effect_allele_frequency'], errors='coerce').fillna(0.0)
        else: prs_df['Pop_AF'] = 0.0
    elif args.af_mode == "NONE":
        prs_df['Pop_AF'] = 0.0
    else:
        prs_df['Pop_AF'] = prs_df.apply(lambda r: get_robust_af(r, gnomad_afs, assume_forward, effect_is_alt), axis=1)
                
    prs_df['LCR_Flag'] = check_lcr(prs_df, args.lcr, pos_col='active_pos')
    prs_df['Palindromic'] = [is_palindromic(eff, oth) for eff, oth in zip(prs_df['effect_allele'].astype(str), prs_df['clean_other'].astype(str))]

    # Scan the input to see how it represents homozygous reference sites (e.g. gVCF blocks vs standard)
    with pysam.VariantFile(args.vcf) as convention_vcf:
        convention_samples = list(convention_vcf.header.samples)
        if not convention_samples:
            raise ValueError("Input VCF has no sample column")
            
        actual_sample_id = convention_samples[0]
        if input_is_gvcf:
            has_explicit_hom_ref, has_reference_blocks = False, True
        else:
            has_explicit_hom_ref, has_reference_blocks = inspect_vcf_genotype_convention(
                convention_vcf, actual_sample_id
            )
            
    input_is_gvcf = input_is_gvcf or has_reference_blocks
    assume_absent_reference = (not input_is_gvcf) and (not has_explicit_hom_ref)

    mode_text = (
        "gVCF: use explicit variants/reference blocks; no absent-site synthesis"
        if input_is_gvcf else
        "standard VCF with explicit 0/0 records: absent sites remain missing for imputation"
        if has_explicit_hom_ref else
        "standard variant-only VCF with no 0/0 records: absent target sites assume reference"
    )
    print(f"Input genotype convention: {mode_text}")

    # Set up VCF I/O headers
    vcf_in = pysam.VariantFile(args.vcf)
    ancestry = get_somalier_ancestry(args.somalier, args.sample)
    reference_fasta = None
    
    if args.ref_fasta and args.ref_fasta.upper() != 'NONE':
        if not os.path.isfile(args.ref_fasta):
            raise FileNotFoundError(f'Reference FASTA not found: {args.ref_fasta}')
        reference_fasta = pysam.FastaFile(args.ref_fasta)

    clean_header = pysam.VariantHeader()
    clean_header.add_line('##fileformat=VCFv4.2')
    clean_header.add_line('##FORMAT=<ID=RGT,Number=1,Type=String,Description="Raw genotype assumption before effect-allele recoding">')
    for contig in vcf_in.header.contigs.values(): clean_header.add_line(str(contig.header_record))
    
    clean_header.add_line(f'##SomalierAncestry={ancestry}')
    policy_tag = (
        "GVCF_REFERENCE_BLOCKS_ONLY" if input_is_gvcf else
        "EXPLICIT_HOM_REF_RECORDS_MISSING_SITES_IMPUTED" if has_explicit_hom_ref else
        "NO_HOM_REF_RECORDS_ABSENT_SITES_ASSUME_REFERENCE"
    )
    
    clean_header.add_line(f'##InputGenotypeConvention={policy_tag}')
    clean_header.add_line('##ALT=<ID=NON_REF,Description="Unresolved allele orientation; not for imputation or scoring">')
    clean_header.add_line('##INFO=<ID=EFF,Number=1,Type=Float,Description="Effect weight">')
    clean_header.add_line('##INFO=<ID=AF,Number=1,Type=Float,Description="Population Allele Frequency">')
    clean_header.add_line('##INFO=<ID=LCR,Number=0,Type=Flag,Description="Low Complexity Region">')
    clean_header.add_line('##INFO=<ID=PAL,Number=0,Type=Flag,Description="Palindromic variant">')
    clean_header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    clean_header.add_line('##FORMAT=<ID=DS,Number=1,Type=Float,Description="Dosage">')
    clean_header.add_line('##FORMAT=<ID=NOTE,Number=1,Type=String,Description="Variant Notes">')
    clean_header.add_sample(actual_sample_id)

    vcf_out = pysam.VariantFile(args.out_vcf, 'w', header=clean_header)
    chroms_to_impute = set()
    vcf_contigs = set(vcf_in.header.contigs.keys())

    # ==========================================================================
    # CORE PROCESSING LOOP
    # ==========================================================================
    for _, row in prs_df.iterrows():
        norm_c, pos = str(row['hm_chr_norm']), int(row['active_pos'])
        action = actions_dict.get(f"{norm_c}:{pos}", "Keep")
        if action == "Drop": continue
        
        c_chr = f"chr{norm_c}"
        if c_chr in vcf_contigs: chrom = c_chr
        elif norm_c in vcf_contigs: chrom = norm_c
        else: continue

        effect_allele = str(row['effect_allele']).strip().upper()
        other_allele = str(row['clean_other']).strip().upper()
        if other_allele.strip().upper() in ('N', '', '.', 'NAN', 'NONE', '<NA>'): 
            other_allele = ''
        expected_ref = other_allele

        dosage, actual_gt, notes_list = None, None, []
        assumed_reference = False
        assumption_orientation_unresolved = False
        if row.get('LCR_Flag', False): notes_list.append("LCR")
        pal_flag = bool(row.get('Palindromic', False))
        saw_record, qc_failed = False, False

        selected_rec = None
        if action == "Force Impute":
            saw_record = False
            qc_failed = False
            notes_list.append("Force_Impute_Requested")
        else:
            try:
                overlapping_records = [
                    candidate
                    for candidate in vcf_in.fetch(region=f"{chrom}:{pos}-{pos}")
                    if candidate.start <= pos - 1 < candidate.stop
                ]
            except (ValueError, OSError):
                overlapping_records = []

            record_kind, target_records = select_target_records(
                overlapping_records, pos, other_allele
            )

            # Route processing based on the type of record discovered at the locus
            if record_kind == "concrete":
                saw_record = True
                if len(target_records) != 1:
                    notes_list.append("AmbiguousConcreteRecordsAtPosition")
                    qc_failed = True
                else:
                    selected_rec = target_records[0]
                    
            elif record_kind == "reference_block":
                saw_record = True
                if len(target_records) != 1:
                    notes_list.append("AmbiguousReferenceBlocks")
                    qc_failed = True
                else:
                    block = target_records[0]
                    block_ok, block_notes = reference_block_status(
                        block, actual_sample_id,
                        min_depth=args.min_depth, min_gq=args.min_gq,
                    )
                    notes_list.extend(block_notes)
                    if block_ok:
                        role, role_note = reference_block_effect_role(
                            reference_fasta, chrom, pos, effect_allele, other_allele,
                            assume_forward=assume_forward,
                            effect_is_alt=effect_is_alt,
                            known_ref_allele=get_harmonized_known_ref(
                                row, is_hg19=is_hg19, chrom=norm_c, pos=pos
                            ),
                        )
                        notes_list.append(role_note)
                        
                        if role is None:
                            qc_failed = True
                        else:
                            if not other_allele:
                                inferred_ref = infer_ref_for_effect_alt_from_fasta(
                                    reference_fasta, chrom, pos, effect_allele
                                ) if effect_is_alt and role == "other" else None
                                
                                if inferred_ref is None:
                                    notes_list.append("EffectOnlyReferenceBlockAlleleUnresolved")
                                    qc_failed = True
                                    role = None
                                else:
                                    expected_ref = inferred_ref
                                    notes_list.extend([
                                        "EffectOnlyMatch",
                                        "OtherAlleleInferredFromReferenceSequence",
                                    ])

                            if role is not None:
                                raw_gt = block.samples[actual_sample_id].get("GT")
                                if raw_gt is None or len(raw_gt) == 0:
                                    block_gt = (0, 0)
                                else:
                                    block_gt = tuple(raw_gt)
                                    
                                ploidy = len(block_gt)
                                if role == "effect":
                                    actual_gt = tuple(1 for _ in block_gt)
                                    dosage = float(ploidy)
                                    notes_list.append("ReferenceBlock_EffectAlleleIsReference")
                                else:
                                    actual_gt = tuple(0 for _ in block_gt)
                                    dosage = 0.0
                                    notes_list.append("ReferenceBlock_OtherAlleleIsReference")
                    else:
                        qc_failed = True
                        
            elif record_kind == "overlap_only":
                saw_record = True
                notes_list.append("OverlappingVariant_NotAtTargetStart")
                qc_failed = True
                
            elif record_kind == "missing" and assume_absent_reference:
                assumed_reference = True
                notes_list.extend([
                    "AbsentFromStandardVCF_AssumedReference",
                    "InputVCF_NoHomozygousReferenceRecords",
                ])

                harmonized_known_ref = get_harmonized_known_ref(
                    row, is_hg19=is_hg19, chrom=norm_c, pos=pos
                )
                synthetic_pair = infer_reference_allele_pair_from_fasta(
                    reference_fasta, chrom, pos, effect_allele, other_allele,
                    assume_forward=assume_forward,
                    effect_is_alt=effect_is_alt,
                    known_ref_allele=harmonized_known_ref,
                )

                if synthetic_pair is not None:
                    synthetic_ref, synthetic_alt = synthetic_pair
                    synthetic_rec = vcf_in.header.new_record(
                        contig=chrom, start=pos - 1, alleles=(synthetic_ref, synthetic_alt),
                    )
                    synthetic_rec.samples[actual_sample_id]["GT"] = (0, 0)

                    selected_rec = synthetic_rec
                    saw_record = True

                    if not other_allele:
                        expected_ref = synthetic_ref
                        notes_list.append("OtherAlleleInferredFromReferenceSequence")
                else:
                    assumption_orientation_unresolved = True
                    notes_list.append(
                        "AssumedReference_AlleleOrientationUnresolved"
                    )

            # Analyze the selected record to determine standard Genotype coding
            if selected_rec is not None:
                rec = selected_rec
                actual_vcf_ref = str(rec.ref or "").strip().upper()
                allele_set = [actual_vcf_ref] + [
                    str(a).strip().upper() for a in (rec.alts or ())
                ]

                match_vcf_idx = -1
                flipped = False
                other_mismatch = False
                effect_only = not other_allele
                other_index = -1

                if effect_only:
                    result = match_effect_only_concrete_record(
                        rec, effect_allele,
                        allow_complement=not assume_forward,
                        effect_is_alt=effect_is_alt,
                    )
                    idx, effect_flipped, inferred_other_idx, inferred_other, reason = result
                    if idx is None:
                        notes_list.append(reason)
                        qc_failed = True
                    else:
                        match_vcf_idx = idx
                        flipped = bool(effect_flipped)
                        other_index = int(inferred_other_idx)
                        expected_ref = str(inferred_other).strip().upper()
                        notes_list.extend([
                            "EffectOnlyMatch",
                            "OtherAlleleInferredFromVCF" if not effect_is_alt
                            else "OtherAlleleInferredFromVCF_REF",
                        ])

                elif effect_is_alt:
                    exact_idx = exact_effect_alt_index(
                        rec, effect_allele, other_allele
                    )
                    if exact_idx is None:
                        notes_list.append("OtherAlleleMismatch")
                        qc_failed = True
                    else:
                        match_vcf_idx = exact_idx
                        other_index = 0

                else:
                    matches, ambiguous, relaxed = match_variants(
                        effect_allele, other_allele, allele_set,
                        allow_complement=not assume_forward,
                    )
                    ambiguous = len(matches) > 1
                    
                    if relaxed:
                        notes_list.append("OtherAlleleMismatch")
                        qc_failed = True
                    elif not matches:
                        notes_list.append("UnmatchedAllele")
                        qc_failed = True
                    elif len(matches) > 1 or ambiguous:
                        if not is_palindromic(effect_allele, other_allele):
                            notes_list.append("AmbiguousAlleleMatch")
                            qc_failed = True
                        else:
                            pal_flag = True
                            pop_af = float(row.get('Pop_AF', 0.0))
                            eaf = row.get('eaf_clean')
                            eaf = None if pd.isna(eaf) else float(eaf)
                            
                            candidate_afs = [
                                {
                                    'idx': idx,
                                    'flipped': is_flip,
                                    'pop_af': (pop_af if not is_flip else 1.0 - pop_af)
                                    if pop_af > 0.0 else None,
                                }
                                for idx, is_flip in matches
                            ]
                            sel_idx, sel_flip, reason = resolve_palindromic(
                                eaf, candidate_afs,
                                af_low=args.pal_low, af_high=args.pal_high,
                            )
                            notes_list.append(reason)
                            if sel_idx is not None:
                                match_vcf_idx, flipped = sel_idx, sel_flip
                            else:
                                qc_failed = True
                    else:
                        match_vcf_idx, flipped = matches[0]
                        try:
                            other_native = (
                                reverse_complement(other_allele)
                                if flipped else other_allele
                            )
                            other_index = allele_set.index(
                                str(other_native).strip().upper()
                            )
                        except ValueError:
                            notes_list.append("UnmatchedAllele")
                            qc_failed = True

                # Perform QC mapping once the allele is resolved
                if match_vcf_idx >= 0:
                    notes_list = [n for n in notes_list if n not in QC_NOTES]
                    raw_gt = rec.samples[actual_sample_id].get('GT', (None, None))
                    
                    if raw_gt is None or len(raw_gt) == 0 or any(g is None for g in raw_gt):
                        notes_list.append("MissingGT")
                        qc_failed = True
                    elif not variant_qa_checker(
                        rec, actual_sample_id, match_vcf_idx, notes_list,
                        depth=args.min_depth, min_gq=args.min_gq, min_af=args.min_af,
                    ):
                        qc_failed = True
                    else:
                        if flipped:
                            notes_list.append("StrandFlipped")
                        if other_mismatch:
                            notes_list.append("OtherAlleleMismatch")

                        if effect_is_alt:
                            other_index = 0
                            if "Harmonization_Bypassed" not in notes_list:
                                notes_list.append("Harmonization_Bypassed")
                        
                        if other_index < 0:
                            notes_list.append("UnmatchedAllele")
                            qc_failed = True
                        else:
                            actual_gt, dosage, gt_problem = recode_gt_to_effect(
                                raw_gt, match_vcf_idx, other_index
                            )
                            if gt_problem:
                                notes_list.append(gt_problem)
                                qc_failed = True

        # Generate the output standardized VCF record
        out_alt = effect_allele
        out_ref = expected_ref
        
        if assumption_orientation_unresolved:
            if reference_fasta is None:
                raise ValueError(f"Cannot represent unresolved site without target FASTA: {chrom}:{pos}")
            refs = set(reference_fasta.references)
            c_norm = chrom[3:] if chrom.lower().startswith('chr') else chrom
            fasta_contig = next((c for c in (chrom, c_norm, f'chr{c_norm}') if c in refs), None)
            if fasta_contig is None:
                raise ValueError(f"Target FASTA contig missing for unresolved site: {chrom}:{pos}")
            out_ref = reference_fasta.fetch(fasta_contig, pos - 1, pos).upper()
            if not is_concrete_allele(out_ref):
                raise ValueError(f"Invalid target FASTA REF for unresolved site: {chrom}:{pos}")
            out_alt = '<NON_REF>'
        
        if out_ref == out_alt:
            if selected_rec is not None and selected_rec.alts:
                concrete_alts = [a for a in selected_rec.alts if is_concrete_allele(a)]
                if concrete_alts:
                    out_ref = str(concrete_alts[0]).strip().upper()
                else:
                    out_ref = 'N'
            else:
                out_ref = 'N'
                
        if out_ref == '': out_ref = 'N'
        if out_alt == '': out_alt = 'N'
        
        if out_ref == out_alt:
            if out_alt == 'N': out_ref = 'A'
            else: out_ref = 'N'

        new_rec = vcf_out.new_record(contig=chrom, start=pos - 1, stop=pos, alleles=(out_ref, out_alt))
        if assumed_reference:
            new_rec.samples[actual_sample_id]["RGT"] = "0/0"
            
        new_rec.info['EFF'] = float(row.get('effect_weight', 0.0))
        new_rec.info['AF'] = float(row['Pop_AF'])
        
        if row.get('LCR_Flag', False): 
            new_rec.info['LCR'] = True
        if pal_flag: 
            new_rec.info['PAL'] = True
        if pal_flag and "Palindromic" not in notes_list: 
            notes_list.append("Palindromic")

        if assumption_orientation_unresolved:
            notes_list.append("HarmonizationOrientationUnresolved_NoImputation")
            new_rec.samples[actual_sample_id]["GT"] = (None, None)
            new_rec.samples[actual_sample_id]["DS"] = None
            new_rec.samples[actual_sample_id]["NOTE"] = ",".join(notes_list)

        elif dosage is not None:
            new_rec.samples[actual_sample_id]["GT"] = actual_gt
            new_rec.samples[actual_sample_id]["DS"] = dosage
            new_rec.samples[actual_sample_id]["NOTE"] = (",".join(notes_list) if notes_list else None)
        else:
            new_rec.samples[actual_sample_id]['GT'] = (None, None)
            if not qc_failed:
                if saw_record:
                    failure_notes = (
                        "UnmatchedAllele", "OtherAlleleMismatch", "MissingGT",
                        "UnexpectedOtherALTGenotype", "AmbiguousConcreteRecordsAtPosition",
                        "AmbiguousReferenceBlocks", "ReferenceBlock_", "OverlappingVariant_",
                    )
                    if not any(any(token in n for token in failure_notes) for n in notes_list):
                        if "AmbiguousPalindromic" not in "".join(notes_list):
                            notes_list.append("UnmatchedAllele")
                elif action != "Force Impute": 
                    notes_list.append("MissingFromVCF")
            if "Target_Imputation" not in notes_list: 
                notes_list.append("Target_Imputation")
                
            new_rec.samples[actual_sample_id]['NOTE'] = ",".join(notes_list)
            chroms_to_impute.add(chrom)

        vcf_out.write(new_rec)
        
    vcf_out.close(); vcf_in.close()
    if reference_fasta is not None:
        reference_fasta.close()
        
    # Standardize VCF header descriptions
    out_dir = os.path.dirname(os.path.abspath(args.out_vcf))
    fd, tmp_path = tempfile.mkstemp(
        dir=out_dir, prefix=os.path.basename(args.out_vcf) + ".", suffix=".tmp"
    )
    try:
        with open(args.out_vcf, 'r') as f_in, os.fdopen(fd, 'w') as f_out:
            for line in f_in:
                if line.startswith('##FORMAT=<ID=DS'): line = re.sub(r'Description=Dosage', 'Description="Dosage"', line)
                elif line.startswith('##FORMAT=<ID=GT'): line = re.sub(r'Description=Genotype', 'Description="Genotype"', line)
                elif line.startswith('##FORMAT=<ID=NOTE'): line = re.sub(r'Description=Variant Notes(?: and Filters)?', 'Description="Variant Notes"', line)
                f_out.write(line)
        try:
            os.chmod(tmp_path, 0o644)   # Adjust permissions safely
        except OSError:
            pass
        os.replace(tmp_path, args.out_vcf)
    except BaseException:
        if os.path.exists(tmp_path):
            os.unlink(tmp_path)
        raise

    # Write target chromosomes for beagle imputation
    with open(args.out_chroms, 'w') as f:
        if not skip_beagle:
            for c in sorted(list(chroms_to_impute)): 
                f.write(f"{c}\n")

if __name__ == "__main__":
    main()