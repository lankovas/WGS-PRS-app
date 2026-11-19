import argparse
import os
import pandas as pd
import numpy as np
import pysam
import math
import sys
import scipy.stats as stats  
import json

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pgs_harmonise import is_palindromic


# ==============================================================================
# 1. HELPER FUNCTIONS
# ==============================================================================

def clean_chrom(c):
    """Normalizes chromosome names (e.g., 'chr1' or '1.0' becomes '1')."""
    if pd.isna(c): return ""
    c = str(c).strip()
    if c.lower().startswith('chr'): c = c[3:]
    if c.endswith('.0'): c = c[:-2]
    return c
    
def load_variant_actions(actions_file):
    """Loads manual variant overrides (e.g., 'Keep', 'Drop', 'Force Impute') set by the user."""
    actions_dict = {}
    if actions_file and os.path.exists(actions_file):
        try:
            df = pd.read_csv(actions_file, sep='\t')
            if 'Action' in df.columns:
                chr_col = 'chr' if 'chr' in df.columns else ('hm_chr' if 'hm_chr' in df.columns else None)
                pos_col = 'pos' if 'pos' in df.columns else ('hm_pos' if 'hm_pos' in df.columns else ('chr_position' if 'chr_position' in df.columns else None))
                if chr_col and pos_col:
                    for _, row in df.iterrows():
                        c = clean_chrom(row[chr_col])
                        p = str(row[pos_col]).strip()
                        actions_dict[f"{c}:{p}"] = str(row['Action']).strip()
        except Exception:
            pass
    return actions_dict

def get_somalier_ancestry(master_somalier_file, sample_id):
    """Extracts the predicted ancestry, confidence, and notes for the patient from Somalier output."""
    df = pd.read_table(master_somalier_file)
    id_col = '#sample_id' if '#sample_id' in df.columns else 'sample_id'
    patient_row = df[df[id_col] == sample_id]
    
    notes = "N/A"
    if 'ancestry_notes' in df.columns and not patient_row.empty:
        notes = str(patient_row['ancestry_notes'].iloc[0]).strip()
    
    if patient_row.empty: 
        return "EUR", 0.0, notes
        
    ancestry_cols = [c for c in df.columns if c.endswith('_prob')]
    if ancestry_cols: 
        probs = patient_row[ancestry_cols].iloc[0]
        max_prob = float(probs.max())
        return probs.idxmax().split('_')[0], max_prob, notes
            
    if 'predicted_ancestry' in df.columns: 
        return str(patient_row['predicted_ancestry'].iloc[0]).strip(), 1.0, notes
        
    return "EUR", 0.0, notes

def load_cached_frequencies(cache_file_path, ancestry, build_query_col):
    """Loads population Allele Frequencies (AF) specific to the patient's ancestry."""
    if not os.path.exists(cache_file_path) or cache_file_path == "NONE": return {}
    cached_df = pd.read_csv(cache_file_path, sep=None, engine='python')
    if cached_df.empty: return {}
    
    chr_col = next((c for c in cached_df.columns if c.lower() in ['chr', '#chrom', 'chrom']), None)
    ref_col = next((c for c in cached_df.columns if c.lower() == 'ref'), None)
    alt_col = next((c for c in cached_df.columns if c.lower() == 'alt'), None)
    pos_col = next((c for c in cached_df.columns if c.lower() == 'pos'), None)
    id_col = next((c for c in cached_df.columns if c.lower() == 'id'), None)
    
    # Locate the correct AF column for the patient's ancestry (fallback to EUR if missing)
    target_af = f"{ancestry}_af".lower()
    af_col = next((c for c in cached_df.columns if c.lower() == target_af), None)
    if not af_col: af_col = next((c for c in cached_df.columns if c.lower() == "eur_af"), None)
    if not af_col:
        non_af = {chr_col, pos_col, ref_col, alt_col, id_col, 'obs_ct'}
        af_cols = [c for c in cached_df.columns if c not in non_af and c is not None]
        af_col = af_cols[-1] if af_cols else cached_df.columns[-1]

    af_dict = {}
    
    # Standard coordinate-based matching
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
            
            r, a = str(row[ref_col]).strip().upper(), str(row[alt_col]).strip().upper()
            key = f"chr{c}:{p}_{r}>{a}"
            try: val = float(row[af_col])
            except (ValueError, TypeError): val = 0.0
            af_dict[key] = val
        return af_dict

    # Query ID string matching
    if build_query_col in cached_df.columns: query_col = build_query_col
    elif "Original_Query" in cached_df.columns: query_col = "Original_Query"
    else: query_col = cached_df.columns[0]
            
    for _, row in cached_df.iterrows():
        try: val = float(row[af_col])
        except: val = 0.0
        q = str(row[query_col])
        af_dict[q] = val
        
    return af_dict

def check_lcr(prs_df, lcr_bed, pos_col='active_pos'):
    """Checks if variants fall inside Low Complexity Regions (LCR) for QC flagging."""
    lcr_flags = []
    try:
        lcr = pysam.TabixFile(lcr_bed)
        lcr_contigs = set(lcr.contigs)
    except Exception:
        return [False] * len(prs_df)
        
    for _, row in prs_df.iterrows():
        norm_c = str(row['hm_chr_norm'])
        c1, c2 = f"chr{norm_c}", norm_c
        chrom = c1 if c1 in lcr_contigs else (c2 if c2 in lcr_contigs else c1)
        pos = int(row[pos_col])
        try: lcr_flags.append(len(list(lcr.fetch(chrom, pos - 1, pos))) > 0)
        except ValueError: lcr_flags.append(False)
    return lcr_flags

def get_match_idx(rec, expected_alt):
    """Finds the index of the expected effect allele within the VCF record's REF/ALTs."""
    exp = expected_alt.upper()
    if rec.ref.upper() == exp: return 0
    if rec.alts:
        for i, alt in enumerate(rec.alts):
            if str(alt).upper() == exp: return i + 1
    return None

def check_bool(val):
    if pd.isna(val): return False
    return str(val).strip().upper() in ['TRUE', 'T', '1', 'YES', 'Y']


# ==============================================================================
# 2. CORE PRS SCORING LOGIC
# ==============================================================================

def score_and_create_clean_vcf(vcf_in_path, vcf_out_path, prs_df, is_raw_mode, actions_dict):
    """
    Calculates the raw Polygenic Risk Score for the patient, generating a final, 
    clean VCF containing only the variants used in the score.
    Returns tracking metrics including imputed variance for confidence intervals.
    """
    vcf_in = pysam.VariantFile(vcf_in_path)
    actual_sample_id = list(vcf_in.header.samples)[0]
    vcf_contigs = set(vcf_in.header.contigs.keys())
    
    # Load entire patient VCF into memory for rapid access
    vcf_dict = {}
    for rec in vcf_in: vcf_dict[f"{rec.chrom}:{rec.pos}"] = rec
    
    # Set up final annotated VCF headers
    clean_header = pysam.VariantHeader()
    clean_header.add_line('##fileformat=VCFv4.2')
    for contig in vcf_in.header.contigs.values(): clean_header.add_line(str(contig.header_record))
        
    clean_header.add_line('##INFO=<ID=EFF,Number=1,Type=Float,Description="Effect weight">')
    clean_header.add_line('##INFO=<ID=AF,Number=1,Type=Float,Description="Population Allele Frequency">')
    clean_header.add_line('##INFO=<ID=LCR,Number=0,Type=Flag,Description="Falls within Low Complexity Region">')
    clean_header.add_line('##INFO=<ID=PRS_IMPUTED,Number=0,Type=Flag,Description="Variant was Imputed from AF or Beagle">')
    clean_header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    clean_header.add_line('##FORMAT=<ID=DS,Number=1,Type=Float,Description="Dosage">')
    clean_header.add_line('##FORMAT=<ID=NOTE,Number=1,Type=String,Description="Variant Notes">')
    clean_header.add_sample(actual_sample_id)
    
    vcf_out = pysam.VariantFile(vcf_out_path, 'w', header=clean_header)
    
    # Tracking accumulators
    total_raw_prs, min_raw_prs, max_raw_prs = 0.0, 0.0, 0.0
    variants_imputed, variants_beagle, variants_af = 0, 0, 0
    imputed_variance = 0.0
    
    for _, row in prs_df.iterrows():
        norm_c = str(row['hm_chr_norm'])
        pos = int(row['active_pos'])
        action = actions_dict.get(f"{norm_c}:{pos}", "Keep")
        
        if action == "Drop": continue
            
        c_chr = f"chr{norm_c}"
        if c_chr in vcf_contigs: chrom = c_chr
        elif norm_c in vcf_contigs: chrom = norm_c
        else: continue 
            
        ref, expected_alt = str(row['clean_other']).strip().upper(), str(row['effect_allele']).strip().upper()
        weight, af = float(row.get('effect_weight', 0.0)), float(row['Pop_AF'])
        is_lcr = row.get('LCR_Flag', False)
        
        # Panel trait scoring types
        is_dominant = row.get('is_dominant', False)
        is_recessive = row.get('is_recessive', False)
        has_dosage = row.get('has_dosage_weights', False)
        w0, w1, w2 = float(row.get('dosage_0_weight', 0.0)), float(row.get('dosage_1_weight', 0.0)), float(row.get('dosage_2_weight', 0.0))
        
        dosage, actual_gt, notes_str = 0.0, (None, None), ""
        rec = vcf_dict.get(f"{chrom}:{pos}")
        
        # 1. Fetch Genotype and Imputation Information from the VCF record
        if rec is not None:
            actual_gt = rec.samples[actual_sample_id].get('GT', (None, None))
            notes_raw = rec.samples[actual_sample_id].get('NOTE')
            if isinstance(notes_raw, (tuple, list)): notes_str = ",".join(str(n) for n in notes_raw if n is not None)
            else: notes_str = str(notes_raw) if notes_raw is not None else ""

            match_idx = get_match_idx(rec, expected_alt)
            if match_idx is not None:
                ds_val = rec.samples[actual_sample_id].get('DS')
                if isinstance(ds_val, (tuple, list)): ds_val = ds_val[0]
                
                # Assign dosage based on expected allele match
                if ds_val is not None: 
                    alt_ds = float(ds_val)
                    if match_idx == 0: dosage = 2.0 - alt_ds
                    elif match_idx == 1: dosage = alt_ds
                    else: dosage = 0.0
                elif actual_gt != (None, None) and None not in actual_gt: 
                    dosage = float(sum([1 for x in actual_gt if x == match_idx]))

        is_beagle = "Imputed_Beagle" in notes_str
        is_missing = (actual_gt == (None, None) or None in actual_gt)
        is_af = "AF_Imputed" in notes_str or (rec is None) or (is_missing and not is_beagle)

        score = 0.0
        
        # 2. Calculate raw variant score based on data availability
        if is_af and not is_raw_mode:
            # Fallback to Allele Frequency expected value
            f = af
            p0 = (1 - f)**2
            p1 = 2 * f * (1 - f)
            p2 = f**2
            
            if has_dosage: score = w0 * p0 + w1 * p1 + w2 * p2
            elif is_dominant: score = weight * (p1 + p2)
            elif is_recessive: score = weight * p2
            else: score = 2.0 * f * weight
            
            dosage = 2.0 * af 
            if "AF_Imputed" not in notes_str: notes_str = (notes_str + ",AF_Imputed").strip(',')
            
        elif is_af and is_raw_mode:
            # Drop missing variants entirely in strict raw mode
            dosage = 0.0
            if "Dropped_RawMode" not in notes_str: notes_str = (notes_str + ",Dropped_RawMode").strip(',')
            
        else:
            # Direct calculation from VCF or Beagle dosage
            hard_dos = round(dosage)
            if has_dosage:
                if hard_dos == 0: score = w0
                elif hard_dos == 1: score = w1
                elif hard_dos == 2: score = w2
            elif is_dominant: score = weight if hard_dos >= 1 else 0.0
            elif is_recessive: score = weight if hard_dos == 2 else 0.0
            else: score = dosage * weight

        # 3. Create the clean output VCF record
        new_rec = vcf_out.new_record(contig=chrom, start=pos-1, stop=pos, alleles=(ref, expected_alt))
        new_rec.info['EFF'] = weight
        new_rec.info['AF'] = af
        if is_lcr: new_rec.info['LCR'] = True
        
        new_rec.samples[actual_sample_id]['GT'] = actual_gt
        new_rec.samples[actual_sample_id]['DS'] = dosage
        if notes_str: new_rec.samples[actual_sample_id]['NOTE'] = notes_str

        total_raw_prs += score
        
        # 4. Calculate tracking metrics and confidence interval bounds
        if is_af and not is_raw_mode:
            new_rec.info['PRS_IMPUTED'] = True
            variants_imputed += 1
            variants_af += 1
            
            # Theoretical variance of imputed variants for CI calculation
            f = af
            p0 = (1 - f)**2
            p1 = 2 * f * (1 - f)
            p2 = f**2
            
            if has_dosage:
                ex = w0 * p0 + w1 * p1 + w2 * p2
                ex2 = (w0**2) * p0 + (w1**2) * p1 + (w2**2) * p2
                imputed_variance += max(0, ex2 - (ex**2))
            elif is_dominant:
                ex = weight * (p1 + p2)
                ex2 = (weight**2) * (p1 + p2)
                imputed_variance += max(0, ex2 - (ex**2))
            elif is_recessive:
                ex = weight * p2
                ex2 = (weight**2) * p2
                imputed_variance += max(0, ex2 - (ex**2))
            else:
                imputed_variance += (2.0 * af * (1.0 - af)) * (weight ** 2)
                
            if weight < 0: min_raw_prs += 2.0 * weight 
            else: max_raw_prs += 2.0 * weight 
            
        elif is_beagle:
            new_rec.info['PRS_IMPUTED'] = True
            variants_imputed += 1
            variants_beagle += 1
            lower_gt, upper_gt = math.floor(dosage), math.ceil(dosage)
            if weight < 0:
                min_raw_prs += upper_gt * weight
                max_raw_prs += lower_gt * weight
            else:
                min_raw_prs += lower_gt * weight
                max_raw_prs += upper_gt * weight
                
        else:
            min_raw_prs += score
            max_raw_prs += score
            
        vcf_out.write(new_rec)
            
    vcf_out.close()
    vcf_in.close()
    return total_raw_prs, variants_imputed, variants_beagle, variants_af, imputed_variance, min_raw_prs, max_raw_prs


# ==============================================================================
# 3. MAIN EXECUTION PIPELINE
# ==============================================================================

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--build', type=str, default='GRCh38')
    parser.add_argument('--vcf', required=True)
    parser.add_argument('--prs', required=True)
    parser.add_argument('--somalier', required=True)
    parser.add_argument('--lcr', required=True)
    parser.add_argument('--af_cache_file', required=True)
    parser.add_argument('--sample', required=True) 
    parser.add_argument('--out', required=True)
    parser.add_argument('--variant_actions', type=str, default="")
    parser.add_argument('--audit_file', type=str, default="")
    parser.add_argument('--ref_stats_json', default="NONE")
    parser.add_argument('--af_mode', default="1000G") 
    parser.add_argument('--seq_type', type=str, default="WGS")
    parser.add_argument('--ref_cache_dir', type=str, default="NONE")
    parser.add_argument('--somalier_raw', type=str, default="NONE") 
    args = parser.parse_args()

    # Determine coordinate conventions based on genome build
    is_hg19 = args.build in ['hg19', 'GRCh37']
    pos_col = 'chr_position' if is_hg19 else 'hm_pos'
    query_id_col = 'query_id_hg19' if is_hg19 else 'query_id_hg38'

    actions_dict = load_variant_actions(args.variant_actions)

    # 1. Propagate previously dropped variants from the harmonization audit log
    if args.audit_file and os.path.exists(args.audit_file):
        try:
            audit_df = pd.read_csv(args.audit_file, sep='\t')
            dropped_df = audit_df[audit_df['status'] == 'Dropped']
            
            def safe_int(val):
                try: return int(float(val))
                except (ValueError, TypeError): return 0
                
            for _, row in dropped_df.iterrows():
                chrom = clean_chrom(row['chr'])
                p19 = safe_int(row.get('pos_hg19', 0))
                p38 = safe_int(row.get('pos_hg38', 0))
                if p19 > 0: actions_dict[f"{chrom}:{p19}"] = "Drop"
                if p38 > 0: actions_dict[f"{chrom}:{p38}"] = "Drop"
        except Exception as e:
            print(f"Warning: Failed to load audit file: {e}")

    # 2. Load and normalize the target PGS panel
    prs_df = pd.read_csv(args.prs, sep='\t', comment='#')
    
    # Calculate additive weights natively from OR/HR if missing
    if 'effect_weight' not in prs_df.columns:
        if 'or' in prs_df.columns: prs_df['effect_weight'] = np.log(pd.to_numeric(prs_df['or'], errors='coerce'))
        elif 'hr' in prs_df.columns: prs_df['effect_weight'] = np.log(pd.to_numeric(prs_df['hr'], errors='coerce'))

    prs_df['is_dominant'] = prs_df.get('is_dominant', pd.Series([False]*len(prs_df))).apply(check_bool)
    prs_df['is_recessive'] = prs_df.get('is_recessive', pd.Series([False]*len(prs_df))).apply(check_bool)
    if 'dosage_0_weight' in prs_df.columns and 'dosage_1_weight' in prs_df.columns and 'dosage_2_weight' in prs_df.columns:
        prs_df['has_dosage_weights'] = prs_df.apply(lambda r: pd.notna(r.get('dosage_0_weight')) and pd.notna(r.get('dosage_1_weight')) and pd.notna(r.get('dosage_2_weight')), axis=1)
    else:
        prs_df['has_dosage_weights'] = False
        
    prs_df['is_additive'] = ~(prs_df['is_dominant'] | prs_df['is_recessive'] | prs_df['has_dosage_weights'])

    if 'hm_chr' not in prs_df.columns and 'chr' in prs_df.columns: prs_df['hm_chr'] = prs_df['chr']
    if pos_col not in prs_df.columns and 'pos' in prs_df.columns: prs_df[pos_col] = prs_df['pos']
    
    if pos_col not in prs_df.columns and is_hg19:
        print(f"Error: Build {args.build} requested, but '{pos_col}' is missing from the PGS file.")
        sys.exit(1)
        
    prs_df['hm_chr_norm'] = prs_df['hm_chr'].apply(clean_chrom)
    prs_df = prs_df[prs_df['hm_chr_norm'] != ""] 
    
    prs_df['active_pos'] = pd.to_numeric(prs_df[pos_col], errors='coerce').fillna(0).astype(int)
    prs_df = prs_df[prs_df['active_pos'] > 0]
    
    # 3. Retrieve ancestry data to inform allele frequency calculations
    ancestry, ancestry_confidence, ancestry_notes = get_somalier_ancestry(args.somalier, args.sample)
    gnomad_afs = load_cached_frequencies(args.af_cache_file, ancestry, build_query_col=query_id_col)
    
    is_raw_mode = (args.af_mode == "NONE")

    if 'hm_inferotherallele' in prs_df.columns: prs_df['clean_other'] = prs_df['hm_inferotherallele'].fillna('').astype(str).str.split('/').str[0].str.upper()
    elif 'other_allele' in prs_df.columns: prs_df['clean_other'] = prs_df['other_allele'].fillna('').astype(str).str.upper()
    elif 'ref' in prs_df.columns: prs_df['clean_other'] = prs_df['ref'].fillna('').astype(str).str.upper()
    else: prs_df['clean_other'] = 'N'
    
    prs_df.loc[prs_df['clean_other'] == '', 'clean_other'] = 'N'
    
    if 'effect_allele' not in prs_df.columns and 'alt' in prs_df.columns:
        prs_df['effect_allele'] = prs_df['alt']

    # Attach AF expectations to each panel row
    if args.af_mode == "PGS_FILE":
        if 'allelefrequency_effect' in prs_df.columns: prs_df['Pop_AF'] = pd.to_numeric(prs_df['allelefrequency_effect'], errors='coerce').fillna(0.0)
        elif 'hm_effect_allele_frequency' in prs_df.columns: prs_df['Pop_AF'] = pd.to_numeric(prs_df['hm_effect_allele_frequency'], errors='coerce').fillna(0.0)
        else: prs_df['Pop_AF'] = 0.0
    elif not is_raw_mode:
        prs_df['var_id_key'] = "chr" + prs_df['hm_chr_norm'] + ":" + prs_df['active_pos'].astype(str) + "_" + prs_df['clean_other'] + ">" + prs_df['effect_allele'].astype(str).str.upper()
        prs_df['Pop_AF'] = prs_df.apply(lambda r: gnomad_afs.get(r.get(query_id_col, r.get('var_id_key')), 0.0), axis=1)
    else:
        prs_df['Pop_AF'] = 0.0

    prs_df['Palindromic'] = [is_palindromic(eff, oth) for eff, oth in zip(prs_df.get('effect_allele', pd.Series(['']*len(prs_df))), prs_df.get('clean_other', pd.Series(['']*len(prs_df))))]
    prs_df['LCR_Flag'] = check_lcr(prs_df, args.lcr, pos_col='active_pos')
    
    out_vcf = args.out.replace('.csv', '_scored_variants.vcf')
    
    # 4. Execute the main patient scoring loop
    patient_raw_prs, variants_imputed, variants_beagle, variants_af, imputed_variance, min_raw_prs, max_raw_prs = score_and_create_clean_vcf(args.vcf, out_vcf, prs_df, is_raw_mode, actions_dict)
    
    pop_mu, pop_sigma, z_score, z_score_percentile, z_score_sd = ["N/A"] * 5
    ci_95_lower, ci_95_upper, min_z_score, max_z_score = ["N/A"] * 4
    z2_score, z2_mean, z2_sd = ["N/A"] * 3
    raw_prs_sd = 0.0
    
    # 5. Compute contextual population statistics and Z-Scores
    if not is_raw_mode:
        active_prs_df = prs_df[~prs_df.apply(lambda row: actions_dict.get(f"{row['hm_chr_norm']}:{row['active_pos']}") == "Drop", axis=1)].copy()
        
        pop_mu_val, pop_sigma_val = None, None
        
        # Load empirical reference distribution stats if available
        if args.ref_stats_json != "NONE" and os.path.exists(args.ref_stats_json):
            try:
                with open(args.ref_stats_json, 'r') as f:
                    stats_dict = json.load(f)
                target_pop = ancestry if ancestry in stats_dict else "EUR"
                if target_pop in stats_dict:
                    pop_mu_val = stats_dict[target_pop]["mean"]
                    pop_sigma_val = stats_dict[target_pop]["sd"]
            except Exception: pass

        # Fallback to computing theoretical population stats directly from the panel and AFs
        if pop_mu_val is None or pop_sigma_val is None:
            pop_mu_val, pop_var_val = 0.0, 0.0
            adjust_all = True
        else:
            pop_var_val = pop_sigma_val ** 2
            adjust_all = False
            
        for _, row in active_prs_df.iterrows():
            if adjust_all or not row['is_additive']:
                f = float(row['Pop_AF'])
                w = float(row.get('effect_weight', 0.0))
                p0, p1, p2 = (1 - f)**2, 2 * f * (1 - f), f**2
                
                if row.get('has_dosage_weights', False):
                    w0, w1, w2 = float(row.get('dosage_0_weight', 0.0)), float(row.get('dosage_1_weight', 0.0)), float(row.get('dosage_2_weight', 0.0))
                    ex = w0 * p0 + w1 * p1 + w2 * p2
                    ex2 = (w0**2) * p0 + (w1**2) * p1 + (w2**2) * p2
                elif row.get('is_dominant', False):
                    ex = w * (p1 + p2)
                    ex2 = (w**2) * (p1 + p2)
                elif row.get('is_recessive', False):
                    ex = w * p2
                    ex2 = (w**2) * p2
                else:
                    ex = 2.0 * f * w
                    ex2 = (w**2) * p1 + (4.0 * w**2) * p2
                    
                pop_mu_val += ex
                pop_var_val += max(0, ex2 - (ex**2))
                
        pop_sigma_val = math.sqrt(pop_var_val) if pop_var_val > 0 else 1.0
        pop_mu, pop_sigma = pop_mu_val, pop_sigma_val

        # Primary categorical Z1 Score
        raw_prs_sd = math.sqrt(imputed_variance)
        z_score_val = (patient_raw_prs - pop_mu) / pop_sigma
        z_score = z_score_val
        z_score_sd = raw_prs_sd / pop_sigma
        z_score_percentile = stats.norm.cdf(z_score_val) * 100
        ci_95_lower = z_score_val - (1.96 * z_score_sd)
        ci_95_upper = z_score_val + (1.96 * z_score_sd)
        min_z_score = (min_raw_prs - pop_mu) / pop_sigma
        max_z_score = (max_raw_prs - pop_mu) / pop_sigma
        total_eval_variants = len(active_prs_df)

        # 6. Compute Admixture-Adjusted Z2 Score using Principal Component linear regression
        if args.seq_type.startswith("Whole Genome") and args.ref_cache_dir != "NONE":
            p_sscore_path = os.path.join(args.ref_cache_dir, "ref_scores_primary.sscore")
            f_sscore_path = os.path.join(args.ref_cache_dir, "ref_scores_fallback.sscore")
            
            ref_scores = None
            for path in [p_sscore_path, f_sscore_path]:
                if os.path.exists(path):
                    df_s = pd.read_csv(path, sep='\t')
                    df_s.columns = [c.replace('#', '') for c in df_s.columns]
                    if ref_scores is None:
                        ref_scores = df_s[['IID', 'SCORE1_SUM']].copy()
                    else:
                        ref_scores = pd.merge(ref_scores, df_s[['IID', 'SCORE1_SUM']], on='IID', how='outer')
                        ref_scores['SCORE1_SUM'] = ref_scores['SCORE1_SUM_x'].fillna(0) + ref_scores['SCORE1_SUM_y'].fillna(0)
                        ref_scores = ref_scores[['IID', 'SCORE1_SUM']]
            
            if ref_scores is not None and not ref_scores.empty and args.somalier_raw != "NONE" and os.path.exists(args.somalier_raw):
                somalier_df = pd.read_table(args.somalier_raw)
                id_col = '#sample_id' if '#sample_id' in somalier_df.columns else 'sample_id'
                
                merged_ref = ref_scores.merge(somalier_df, left_on='IID', right_on=id_col, how='inner')
                patient_row = somalier_df[somalier_df[id_col] == args.sample]
                
                if not merged_ref.empty and not patient_row.empty and 'PC1' in merged_ref.columns:
                    X = merged_ref[['PC1', 'PC2', 'PC3', 'PC4', 'PC5']].fillna(0).values
                    X = np.c_[np.ones(X.shape[0]), X]
                    Y = merged_ref['SCORE1_SUM'].values
                    beta, _, _, _ = np.linalg.lstsq(X, Y, rcond=None)
                    
                    residuals = Y - X.dot(beta)
                    Y_var = residuals ** 2
                    beta_var, _, _, _ = np.linalg.lstsq(X, Y_var, rcond=None)
                    
                    p_pcs = patient_row[['PC1', 'PC2', 'PC3', 'PC4', 'PC5']]
                    if not p_pcs.isna().all().all():
                        X_p = np.concatenate(([1], p_pcs.fillna(0).iloc[0].values))
                        
                        pred_mean = X_p.dot(beta)
                        pred_var = X_p.dot(beta_var)
                        pred_sd = math.sqrt(pred_var) if pred_var > 0 else np.std(residuals)
                        
                        if pred_sd > 0:
                            z2_val = (patient_raw_prs - pred_mean) / pred_sd
                            z2_score = z2_val
                            z2_mean = pred_mean
                            z2_sd = pred_sd
                            
    else:
        total_eval_variants = sum(1 for _, row in prs_df.iterrows() if actions_dict.get(f"{row['hm_chr_norm']}:{row['active_pos']}") != "Drop")

    # 7. Write Final Result CSV Data
    pal_count = prs_df.get('Palindromic', pd.Series([False]*len(prs_df))).sum()
    dropped_vars = len(prs_df) - total_eval_variants
    
    with open(args.out, 'w') as f:
        f.write("Metric,Value\n")
        f.write(f"Ancestry,{ancestry}\n")
        f.write(f"Ancestry Confidence,{ancestry_confidence:.2f}\n")
        f.write(f"Ancestry Notes,{ancestry_notes}\n")
        f.write(f"Population Mean,{pop_mu}\n")
        f.write(f"Population SD,{pop_sigma}\n")
        f.write(f"Raw PRS,{patient_raw_prs:.5f}\n")
        f.write(f"SD of Raw PRS,{raw_prs_sd}\n")
        f.write(f"Percentile,{z_score_percentile}\n")
        f.write(f"Categorical Z-Score (Z1),{z_score}\n")
        f.write(f"SD of Z-Score,{z_score_sd}\n")
        f.write(f"PC-Adjusted Z-Score (Z2),{z2_score}\n")
        f.write(f"Z2 Predicted Mean,{z2_mean}\n")
        f.write(f"Z2 Predicted SD,{z2_sd}\n")
        f.write(f"Z-Score Lower 95% Uncertainty Interval,{ci_95_lower}\n")
        f.write(f"Z-Score Upper 95% Uncertainty Interval,{ci_95_upper}\n")
        f.write(f"Z-Score Deterministic Minimum,{min_z_score}\n")
        f.write(f"Z-Score Deterministic Maximum,{max_z_score}\n")
        f.write(f"Total Variants in Panel,{len(prs_df)}\n")
        f.write(f"Variants Dropped by Harmonization,{dropped_vars}\n")
        f.write(f"Total Evaluated Variants,{total_eval_variants}\n")
        f.write(f"Number of Imputed Variants,{variants_imputed}\n")
        f.write(f"Number of Imputed Variants by Beagle,{variants_beagle}\n")
        f.write(f"Number of Imputed Variants by AF,{variants_af}\n")
        f.write(f"Number of Variants in Low-Confidence Regions,{prs_df['LCR_Flag'].sum()}\n")
        f.write(f"Total Palindromic Variants in Panel,{pal_count}\n")

if __name__ == "__main__":
    main()