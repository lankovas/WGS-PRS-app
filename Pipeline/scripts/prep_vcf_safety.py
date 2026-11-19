"""Safety helpers for PGS prep VCF allele/genotype handling.

These helpers safely extract and validate genetic information (alleles, genotypes)
from VCF and gVCF files without crashing when expected data is missing or malformed.
"""
from __future__ import annotations

import math
from typing import Any, Iterable

# Placeholders used in VCF files when a real alternate allele isn't present
_PLACEHOLDER_ALLELES = {".", "*", "<*>", "<NON_REF>", "<NONREF>"}


# ==============================================================================
# 1. VALUE EXTRACTION & TYPE SAFETY
# ==============================================================================

def _first(value: Any) -> Any:
    """Extracts the first valid item from a list/tuple, or returns the value if it's already a scalar."""
    if isinstance(value, (tuple, list)):
        for item in value:
            if item is not None:
                return item
        return None
    return value


def _finite_float(value: Any) -> float | None:
    """Safely converts a value to a float, returning None if it's invalid, infinity, or missing."""
    value = _first(value)
    if value is None:
        return None
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def first_nonmissing(*values: Any) -> Any:
    """Returns the first valid, non-empty value from a sequence of possibilities."""
    for value in values:
        if isinstance(value, (tuple, list)):
            value = next((item for item in value if item is not None), None)
        if value is None:
            continue
        if isinstance(value, str) and value.strip() in {"", ".", "NA", "NAN", "NONE"}:
            continue
        try:
            if isinstance(value, float) and math.isnan(value):
                continue
        except TypeError:
            pass
        return value
    return None


def safe_info_get(record: Any, key: str) -> Any:
    """Safely reads INFO fields from a VCF record without crashing if the field is undeclared."""
    try:
        header = getattr(record, "header", None)
        info_header = getattr(header, "info", None) if header is not None else None
        if info_header is not None and key not in info_header:
            return None
        info = getattr(record, "info", None)
        if info is None:
            return None
        return info.get(key)
    except (AttributeError, KeyError, TypeError, ValueError):
        return None


# ==============================================================================
# 2. ALLELE PROPERTIES & SEQUENCE UTILITIES
# ==============================================================================

def is_concrete_allele(allele: Any) -> bool:
    """Checks if the allele is an actual DNA sequence (A, C, G, T) rather than a symbolic placeholder."""
    if allele is None:
        return False
    allele = str(allele).strip().upper()
    if not allele or allele in _PLACEHOLDER_ALLELES:
        return False
    if allele.startswith("<") and allele.endswith(">"):
        return False
    if "[" in allele or "]" in allele:
        return False
    return all(base in "ACGT" for base in allele)


def _reverse_complement(allele: str) -> str:
    """Returns the reverse complement of a DNA sequence (e.g., 'A' becomes 'T')."""
    allele = str(allele or "").strip().upper()
    if not allele or any(base not in "ACGTN" for base in allele):
        return ""
    return allele.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def infer_ref_for_effect_alt_from_fasta(
    fasta: Any, chrom: str, pos: int, effect_allele: str
) -> str | None:
    """Looks up the reference allele from the FASTA genome file when only the alternate is known."""
    effect = str(effect_allele or "").strip().upper()
    if len(effect) != 1 or not is_concrete_allele(effect) or fasta is None:
        return None

    refs = set(getattr(fasta, "references", ()))
    chrom_text = str(chrom).strip()
    chrom_norm = chrom_text[3:] if chrom_text.lower().startswith("chr") else chrom_text
    contig = next((c for c in (chrom_text, chrom_norm, f"chr{chrom_norm}") if c in refs), None)
    if contig is None:
        return None

    try:
        ref = fasta.fetch(contig, int(pos) - 1, int(pos)).upper()
    except (KeyError, OSError, ValueError, TypeError):
        return None

    if len(ref) != 1 or not is_concrete_allele(ref) or ref == effect:
        return None
    return ref


# ==============================================================================
# 3. gVCF REFERENCE BLOCK HANDLING
# ==============================================================================

def is_reference_block(record: Any, target_pos: int, other_allele: str | None) -> bool:
    """Determines if a VCF record represents a 'reference block' (a region matching the reference genome without variants)."""
    target0 = int(target_pos) - 1
    try:
        if not (int(record.start) <= target0 < int(record.stop)):
            return False
    except (AttributeError, TypeError, ValueError):
        return False

    alts = list(record.alts or ())
    if any(is_concrete_allele(alt) for alt in alts):
        return False

    if alts and not any(str(a).strip().upper() in _PLACEHOLDER_ALLELES for a in alts):
        return False

    if int(record.pos) != int(target_pos):
        return int(record.stop) > int(record.start) + 1

    other = str(other_allele or "").strip().upper()
    if other and str(record.ref or "").strip().upper() == other:
        return True
    return int(record.stop) > int(record.start) + 1


def select_target_records(
    overlapping_records: Iterable[Any],
    target_pos: int,
    other_allele: str | None,
) -> tuple[str, list[Any]]:
    """Prioritizes actual variants over reference blocks when multiple records overlap a position."""
    records = list(overlapping_records)
    concrete = [
        rec for rec in records
        if int(rec.pos) == int(target_pos)
        and any(is_concrete_allele(a) for a in (rec.alts or ()))
    ]
    if concrete:
        return "concrete", concrete
    blocks = [rec for rec in records if is_reference_block(rec, target_pos, other_allele)]
    if blocks:
        return "reference_block", blocks
    if records:
        return "overlap_only", []
    return "missing", []


def reference_block_status(
    record: Any,
    sample_id: str,
    min_depth: int = 15,
    min_gq: float = 50.0,
) -> tuple[bool, list[str]]:
    """Validates if a reference block has high enough quality (depth and genotype quality) to be trusted."""
    if sample_id not in record.samples:
        return False, ["ReferenceBlock_MissingSample"]

    sample = record.samples[sample_id]
    gt = sample.get("GT")
    
    inferred_gt = False
    if gt is None or len(gt) == 0 or any(a is None for a in gt):
        inferred_gt = True
    elif any(a != 0 for a in gt):
        return False, ["ReferenceBlock_NonReferenceGT"]

    depth = None
    for value in (
        sample.get("MIN_DP"), sample.get("DP"),
        safe_info_get(record, "MinDP"), safe_info_get(record, "DP"),
    ):
        depth = _finite_float(value)
        if depth is not None:
            break

    if int(min_depth) > 0:
        if depth is None:
            return False, ["ReferenceBlock_MissingDepth"]
        if depth < int(min_depth):
            return False, ["ReferenceBlock_LowDepth"]

    gq = None
    for value in (sample.get("GQ"), sample.get("RGQ")):
        gq = _finite_float(value)
        if gq is not None:
            break

    # Estimate quality from likelihoods if direct GQ is missing
    if gq is None:
        gl = sample.get("GL")
        if gl:
            values = [_finite_float(x) for x in gl]
            values = sorted((x for x in values if x is not None), reverse=True)
            if len(values) >= 2:
                gq = (values[0] - values[1]) * 10.0
    if gq is None:
        pl = sample.get("PL")
        if pl:
            values = [_finite_float(x) for x in pl]
            values = sorted(x for x in values if x is not None)
            if len(values) >= 2:
                gq = values[1] - values[0]

    if gq is not None and gq < float(min_gq):
        return False, ["ReferenceBlock_LowGQ"]

    notes = ["ReferenceBlock_HomozygousReference"]
    if inferred_gt:
        notes.append("ReferenceBlock_InferredMissingGT")
    if gq is None:
        notes.append("ReferenceBlock_GQ_Unavailable")
    return True, notes


# ==============================================================================
# 4. EFFECT MATCHING & SCORING LOGIC
# ==============================================================================

def exact_effect_alt_index(record: Any, effect_allele: str, other_allele: str) -> int | None:
    """Finds the alternate allele index ONLY if the target effect strictly matches the alternate allele."""
    effect = str(effect_allele or "").strip().upper()
    other = str(other_allele or "").strip().upper()
    ref = str(record.ref or "").strip().upper()
    
    if not effect or not other or effect == other or ref != other:
        return None

    indices = [
        i + 1 for i, alt in enumerate(record.alts or ())
        if is_concrete_allele(alt) and str(alt).strip().upper() == effect
    ]
    return indices[0] if len(indices) == 1 else None


def match_effect_only_concrete_record(
    record: Any,
    effect_allele: str,
    *,
    allow_complement: bool = True,
    effect_is_alt: bool = False,
) -> tuple[int, bool, int, str, str] | tuple[None, None, None, None, str]:
    """Matches a VCF record to the panel effect allele even if the 'other' allele is unknown."""
    effect = str(effect_allele or "").strip().upper()
    if not effect:
        return None, None, None, None, "EffectOnlyEffectAlleleMissing"

    ref = str(getattr(record, "ref", "") or "").strip().upper()
    alts = [str(a).strip().upper() for a in (getattr(record, "alts", None) or ())]
    alleles = [ref] + alts

    if effect_is_alt:
        matches = [i for i, alt in enumerate(alts, start=1)
                   if is_concrete_allele(alt) and alt == effect]
        if len(matches) != 1 or not is_concrete_allele(ref):
            return None, None, None, None, "EffectOnlyMatchAmbiguousOrAbsent"
        return matches[0], False, 0, ref, "EffectOnlyMatch"

    if len(alleles) != 2 or not all(is_concrete_allele(a) for a in alleles):
        return None, None, None, None, "EffectOnlyMultiallelicUnresolved"

    candidates = []
    for flipped, query in ((False, effect), (True, _reverse_complement(effect))):
        if flipped and not allow_complement:
            continue
        if query and query in alleles:
            candidates.append((alleles.index(query), flipped))

    candidates = list(dict.fromkeys(candidates))
    if not candidates:
        return None, None, None, None, "UnmatchedAllele"
    if len(candidates) != 1:
        return None, None, None, None, "EffectOnlyMatchAmbiguous"

    effect_index, flipped = candidates[0]
    other_index = 1 - effect_index
    inferred_other = alleles[other_index]
    if flipped:
        inferred_other = _reverse_complement(inferred_other)
        if not inferred_other:
            return None, None, None, None, "EffectOnlyOtherAlleleComplementInvalid"

    return effect_index, flipped, other_index, inferred_other, "EffectOnlyMatch"


def recode_gt_to_effect(
    gt: Iterable[int | None] | None,
    effect_index: int,
    other_index: int,
) -> tuple[tuple[int, ...] | None, float | None, str | None]:
    """Translates standard genotypes (e.g. 0/1) into an effect dosage count (e.g., 0, 1, or 2)."""
    if gt is None:
        return None, None, "MissingGT"
    gt = tuple(gt)
    if not gt or any(a is None for a in gt):
        return None, None, "MissingGT"
    if effect_index == other_index:
        return None, None, "UnmatchedAllele"
    
    allowed = {int(effect_index), int(other_index)}
    if any(int(a) not in allowed for a in gt):
        return None, None, "UnexpectedOtherALTGenotype"
        
    out_gt = tuple(1 if int(a) == int(effect_index) else 0 for a in gt)
    return out_gt, float(sum(out_gt)), None


def reference_block_effect_role(
    fasta: Any,
    chrom: str,
    pos: int,
    effect_allele: str,
    other_allele: str,
    *,
    assume_forward: bool = False,
    effect_is_alt: bool = False,
    known_ref_allele: str | None = None,
) -> tuple[str | None, str]:
    """Determines whether a homozygous reference block counts as having the 'effect' or the 'other' allele."""
    effect = str(effect_allele or "").strip().upper()
    other = str(other_allele or "").strip().upper()

    if not effect:
        return None, "ReferenceBlock_EffectAlleleUnavailable"

    if effect_is_alt and not other:
        return "other", "ReferenceBlock_EffectIsALT_MissingOtherAssumption"

    if not other:
        return None, "ReferenceBlock_OtherAlleleUnavailable"

    if effect == other:
        return None, "ReferenceBlock_AllelePairUnavailable"

    # Strict FASTA verification
    if known_ref_allele is not None:
        known_ref = str(known_ref_allele).strip().upper()
        if fasta is None:
            return None, "ReferenceBlock_KnownREF_NoReferenceFasta"
        refs = set(getattr(fasta, "references", ()))
        chrom_text = str(chrom).strip()
        chrom_norm = chrom_text[3:] if chrom_text.lower().startswith("chr") else chrom_text
        contig = next((c for c in (chrom_text, chrom_norm, f"chr{chrom_norm}") if c in refs), None)
        if contig is None:
            return None, "ReferenceBlock_ReferenceContigMissing"
        try:
            observed = fasta.fetch(contig, int(pos) - 1, int(pos) - 1 + len(known_ref)).upper()
        except (KeyError, OSError, ValueError, TypeError):
            return None, "ReferenceBlock_KnownREF_FetchFailed"
        if observed != known_ref:
            return None, "ReferenceBlock_KnownREF_ReferenceMismatch"
            
        candidates = [(effect, "effect"), (other, "other")]
        roles = {role for allele, role in candidates if allele and allele == known_ref}
        if effect_is_alt:
            roles = {role for role in roles if role == "other"}
        if len(roles) == 1:
            return next(iter(roles)), "ReferenceBlock_HarmonizedREFVerified"
        return None, "ReferenceBlock_KnownREF_NotInAllelePair"

    if effect_is_alt and fasta is None:
        return "other", "ReferenceBlock_EffectIsALT_Assumption"
    if fasta is None:
        return None, "ReferenceBlock_NoReferenceFasta"

    refs = set(getattr(fasta, "references", ()))
    chrom_text = str(chrom).strip()
    chrom_norm = chrom_text[3:] if chrom_text.lower().startswith("chr") else chrom_text
    contig = next((c for c in (chrom_text, chrom_norm, f"chr{chrom_norm}") if c in refs), None)
    if contig is None:
        return None, "ReferenceBlock_ReferenceContigMissing"

    start = int(pos) - 1

    def matches_ref(allele: str) -> bool:
        if not allele:
            return False
        try:
            seq = fasta.fetch(contig, start, start + len(allele)).upper()
        except (KeyError, OSError, ValueError):
            return False
        return seq == allele

    if effect_is_alt:
        if matches_ref(other):
            return "other", "ReferenceBlock_AllelePairReferenceVerified"
        return None, "ReferenceBlock_EffectIsALT_ReferenceMismatch"

    if len(effect) != len(other):
        return None, "ReferenceBlock_IndelOrientationUnresolved"

    possible_roles = set()
    if matches_ref(effect):
        possible_roles.add("effect")
    if matches_ref(other):
        possible_roles.add("other")

    if (
        not assume_forward
        and len(effect) == len(other) == 1
        and effect in "ACGT"
        and other in "ACGT"
    ):
        effect_rc = _reverse_complement(effect)
        other_rc = _reverse_complement(other)
        if matches_ref(effect_rc):
            possible_roles.add("effect")
        if matches_ref(other_rc):
            possible_roles.add("other")

    if len(possible_roles) == 1:
        role = next(iter(possible_roles))
        return role, "ReferenceBlock_AllelePairReferenceVerified"
    if len(possible_roles) > 1:
        return None, "ReferenceBlock_AmbiguousStrandOrAlleles"
    return None, "ReferenceBlock_AllelePairReferenceUnverified"


def inspect_vcf_genotype_convention(records: Iterable[Any], sample_id: str) -> tuple[bool, bool]:
    """Scans the VCF to see if it uses explicit reference blocks (gVCF style) or standard calls."""
    has_hom_ref = False
    has_blocks = False
    
    for record in records:
        try:
            sample = record.samples[sample_id]
        except (KeyError, TypeError, AttributeError):
            sample = None
            
        if sample is not None and not has_hom_ref:
            gt = sample.get("GT")
            if gt is not None and tuple(gt) == (0, 0):
                has_hom_ref = True

        if not has_blocks:
            alts = list(getattr(record, "alts", None) or ())
            has_placeholder = any(str(a).strip().upper() in _PLACEHOLDER_ALLELES for a in alts)
            has_nonref_symbol = any(str(a).strip().upper() in {"<NON_REF>", "<NONREF>", "<*>"} for a in alts)
            has_concrete = any(is_concrete_allele(a) for a in alts)
            
            try:
                spans_multiple_bases = int(record.stop) > int(record.start) + 1
            except (AttributeError, TypeError, ValueError):
                spans_multiple_bases = False
                
            if (has_nonref_symbol and not has_concrete) or (not alts and spans_multiple_bases) or (has_placeholder and not has_concrete and spans_multiple_bases):
                has_blocks = True

        if has_hom_ref and has_blocks:
            break
            
    return has_hom_ref, has_blocks


def infer_reference_allele_pair_from_fasta(
    fasta: Any,
    chrom: str,
    pos: int,
    effect_allele: str,
    other_allele: str,
    *,
    assume_forward: bool = False,
    effect_is_alt: bool = False,
    known_ref_allele=None,
) -> tuple[str, str] | None:
    """Uses the FASTA genome file to reconstruct a valid (REF, ALT) string pair for a missing site."""
    effect = str(effect_allele or "").strip().upper()
    other = str(other_allele or "").strip().upper()
    if fasta is None or not is_concrete_allele(effect):
        return None

    refs = set(getattr(fasta, "references", ()))
    chrom_text = str(chrom).strip()
    chrom_norm = chrom_text[3:] if chrom_text.lower().startswith("chr") else chrom_text
    contig = next((c for c in (chrom_text, chrom_norm, f"chr{chrom_norm}") if c in refs), None)
    if contig is None:
        return None

    def matches_fasta(allele: str) -> bool:
        if not is_concrete_allele(allele):
            return False
        try:
            observed = fasta.fetch(contig, int(pos) - 1, int(pos) - 1 + len(allele)).upper()
        except (KeyError, OSError, ValueError, TypeError):
            return False
        return observed == allele
        
    if known_ref_allele is not None:
        known_ref = str(known_ref_allele).strip().upper()
        allele_pairs = [(effect, other)]
        resolved = set()
        
        for oriented_effect, oriented_other in allele_pairs:
            if effect_is_alt:
                if known_ref == oriented_other:
                    resolved.add((known_ref, oriented_effect))
            else:
                if known_ref == oriented_effect:
                    resolved.add((known_ref, oriented_other))
                elif known_ref == oriented_other:
                    resolved.add((known_ref, oriented_effect))
                    
        if len(resolved) != 1:
            return None
        ref_allele, alt_allele = next(iter(resolved))
        if not is_concrete_allele(ref_allele) or not is_concrete_allele(alt_allele) or not matches_fasta(ref_allele):
            return None
        return ref_allele, alt_allele
    
    if effect_is_alt:
        if other:
            if not is_concrete_allele(other) or other == effect or not matches_fasta(other):
                return None
            return other, effect
        inferred_ref = infer_ref_for_effect_alt_from_fasta(fasta, contig, pos, effect)
        return (inferred_ref, effect) if inferred_ref is not None else None

    if not other or not is_concrete_allele(other) or other == effect:
        return None

    orientations = (False,) if assume_forward else (False, True)
    candidates: set[tuple[str, str]] = set()
    
    for flipped in orientations:
        e = _reverse_complement(effect) if flipped else effect
        o = _reverse_complement(other) if flipped else other
        
        if not e or not o or not is_concrete_allele(e) or not is_concrete_allele(o):
            continue
        if matches_fasta(e):
            candidates.add((e, o))
        if matches_fasta(o):
            candidates.add((o, e))
            
    if len(candidates) != 1:
        return None
    return next(iter(candidates))