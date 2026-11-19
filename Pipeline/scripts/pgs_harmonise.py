"""Strict allele matching and ambiguity resolution for PRS.

This module provides tools to safely match genetic variants (alleles) between 
a scoring file and a reference genome/VCF. It strictly enforces that the 
"other" allele matches the variant identity, avoiding unsafe fallback assumptions. 
It also contains logic to resolve palindromic SNPs (e.g. A/T or C/G) by comparing 
allele frequencies.
"""
from __future__ import annotations

import math
from typing import Any, Iterable

# Mapping used to flip DNA strands to their reverse complement
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")

# Standard representations of missing or unavailable alleles in VCFs
MISSING_ALLELES = {"", "N", ".", "NAN", "NONE", "<NA>"}


# ==============================================================================
# 1. STRING CLEANING & NORMALIZATION
# ==============================================================================

def _clean_allele(value: Any) -> str:
    """Standardizes allele strings to uppercase and converts missing values to empty strings."""
    if value is None:
        return ""
    value = str(value).strip().upper()
    return "" if value in MISSING_ALLELES else value


def reverse_complement(allele: Any) -> str:
    """Returns the reverse complement of a DNA sequence (e.g., 'A' becomes 'T').
    Leaves symbolic alleles (like <DEL> or <NON_REF>) intact.
    """
    allele = _clean_allele(allele)
    if not allele:
        return ""
        
    # Ignore symbolic representations
    if (allele.startswith("<") and allele.endswith(">")) or any(
        token in allele for token in ("[", "]", "*")
    ):
        return allele
        
    # Ignore non-standard bases
    if any(base not in "ACGTN" for base in allele):
        return allele
        
    return allele.translate(COMPLEMENT)[::-1]


def clean_chrom(chrom: Any) -> str:
    """Normalizes chromosome labels to a clean integer/string format (e.g. 'chr1' -> '1', 1.0 -> '1')."""
    try:
        import pandas as pd
        if pd.isna(chrom):
            return ""
    except (ImportError, TypeError, ValueError):
        pass
        
    chrom = str(chrom).strip()
    if chrom.lower().startswith("chr"):
        chrom = chrom[3:]
    if chrom.endswith(".0"):
        chrom = chrom[:-2]
        
    return chrom


# ==============================================================================
# 2. ALLELE MATCHING & IDENTIFICATION
# ==============================================================================

def is_palindromic(effect: Any, other: Any) -> bool:
    """Checks if a variant is palindromic (A/T or C/G), meaning its orientation cannot be determined by bases alone."""
    effect, other = _clean_allele(effect), _clean_allele(other)
    if len(effect) != 1 or len(other) != 1:
        return False
    if effect not in "ACGT" or other not in "ACGT":
        return False
    return reverse_complement(effect) == other


def match_variants(
    effect: Any,
    other: Any,
    allele_set: Iterable[str],
    *,
    allow_complement: bool = True,
):
    """Matches the requested effect/other alleles against the actual alleles present in a VCF record.
    
    Returns a tuple containing:
      - A list of successful matches (each match is a tuple: (index_in_allele_set, was_flipped))
      - A boolean indicating if the match is ambiguous (multiple valid orientations)
      - A boolean indicating if a 'relaxed' match occurred (only the effect allele matched, other was missing)
    """
    effect, other = _clean_allele(effect), _clean_allele(other)
    alleles = [_clean_allele(a) for a in allele_set]
    if not effect or not alleles:
        return [], False, bool(other)

    orientations = (False, True) if allow_complement else (False,)
    
    # 1. If an 'other' allele is provided, BOTH alleles must match the VCF exactly.
    if other:
        pair_matches: list[tuple[int, bool]] = []
        for flipped in orientations:
            e = reverse_complement(effect) if flipped else effect
            o = reverse_complement(other) if flipped else other
            
            if not e or not o or e == o:
                continue
            if e not in alleles or o not in alleles:
                continue
                
            ie, io = alleles.index(e), alleles.index(o)
            # Ensure at least one of the matched alleles is the Reference (index 0)
            if ie == 0 or io == 0:
                pair_matches.append((ie, flipped))

        pair_matches = list(dict.fromkeys(pair_matches))
        if pair_matches:
            return pair_matches, len(pair_matches) > 1, False
        return [], False, True

    # 2. If NO 'other' allele is provided, match purely on the effect allele.
    effect_matches: list[tuple[int, bool]] = []
    for flipped in orientations:
        e = reverse_complement(effect) if flipped else effect
        if e and e in alleles:
            effect_matches.append((alleles.index(e), flipped))
            
    effect_matches = list(dict.fromkeys(effect_matches))
    return effect_matches, len(effect_matches) > 1, False


# ==============================================================================
# 3. ALLELE FREQUENCY RESOLUTION FOR AMBIGUOUS VARIANTS
# ==============================================================================

def _valid_af(value: Any, af_low: float, af_high: float) -> float | None:
    """Validates an Allele Frequency, returning None if it falls inside the ambiguous 'middle' threshold."""
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
        
    if not math.isfinite(value) or not 0.0 <= value <= 1.0:
        return None
    # Reject frequencies too close to 50/50 (e.g. 40% to 60%) as they are unreliable for resolving orientation
    if af_low <= value <= af_high:
        return None
        
    return value


def resolve_palindromic(
    eaf: Any,
    candidate_afs: list[dict[str, Any]],
    label: str = "AmbiguousPalindromic",
    af_low: float = 0.4,
    af_high: float = 0.6,
):
    """Resolves the strand orientation of a palindromic SNP by comparing its 
    expected allele frequency (EAF) to the population frequencies of the two candidate orientations."""
    
    # Check if the expected allele frequency (EAF) is valid for resolution
    try:
        eaf = float(eaf)
    except (TypeError, ValueError):
        return None, None, f"{label}_NoEAF"
        
    if not math.isfinite(eaf) or not 0.0 <= eaf <= 1.0:
        return None, None, f"{label}_NoEAF"
    if af_low <= eaf <= af_high:
        return None, None, f"{label}_5050"

    side = eaf > 0.5
    valid = []
    
    # Filter candidates to only those with valid, extreme-enough population AFs
    for candidate in candidate_afs:
        af = _valid_af(candidate.get("pop_af"), af_low, af_high)
        if af is not None:
            valid.append((candidate, af))
            
    if not valid:
        return None, None, f"{label}_NoPopAF"
        
    # Find candidates whose frequency 'side' (major vs minor) matches the expected EAF
    agreeing = [(candidate, af) for candidate, af in valid if (af > 0.5) == side]
    if not agreeing:
        return None, None, f"{label}_Unresolved"
        
    # Pick the candidate whose frequency is mathematically closest to the expected EAF
    best_dist = min(abs(eaf - af) for _, af in agreeing)
    best = [(candidate, af) for candidate, af in agreeing if abs(abs(eaf - af) - best_dist) <= 1e-9]
    
    if len(best) != 1:
        return None, None, f"{label}_PopConflict"
        
    chosen = best[0][0]
    return chosen.get("idx"), chosen.get("flipped"), f"{label}_Resolved"


def resolve_candidates_global(
    eaf: Any,
    candidates: list[dict[str, Any]],
    label: str = "AmbiguousVariant",
    af_low: float = 0.4,
    af_high: float = 0.6,
    preferred_pop: str | None = None,
):
    """Resolves allele orientations by comparing expected allele frequency across multiple global populations.
    
    Note: This function assumes all candidates represent orientations of the SAME variant record, 
    not different records that happen to overlap.
    """
    try:
        eaf = float(eaf)
    except (TypeError, ValueError):
        return None, f"{label}_NoEAF"
        
    if not math.isfinite(eaf) or not 0.0 <= eaf <= 1.0:
        return None, f"{label}_NoEAF"
    if af_low <= eaf <= af_high:
        return None, f"{label}_5050"
        
    side = eaf > 0.5

    # 1. Attempt to resolve using a preferred population (e.g. if the patient's ancestry is known)
    if preferred_pop:
        preferred = []
        for candidate in candidates:
            af = _valid_af(candidate.get("pop_afs", {}).get(preferred_pop), af_low, af_high)
            if af is not None:
                preferred.append((candidate, af))
                
        if preferred:
            agreeing = [(c, af) for c, af in preferred if (af > 0.5) == side]
            if not agreeing:
                return None, f"{label}_Unresolved_{preferred_pop}"
                
            best_dist = min(abs(eaf - af) for _, af in agreeing)
            best = [(c, af) for c, af in agreeing if abs(abs(eaf - af) - best_dist) <= 1e-9]
            if len(best) != 1:
                return None, f"{label}_PopConflict_{preferred_pop}"
            return best[0][0]["key"], f"{label}_Resolved_{preferred_pop}"

    # 2. Fallback: Aggregate and score across all available populations
    scored = []
    for candidate in candidates:
        afs = [
            af for value in candidate.get("pop_afs", {}).values()
            if (af := _valid_af(value, af_low, af_high)) is not None
        ]
        if not afs:
            continue
            
        agree = sum(1 for af in afs if (af > 0.5) == side)
        if agree == 0:
            continue
            
        mean_dist = sum(abs(eaf - af) for af in afs) / len(afs)
        scored.append((agree, mean_dist, candidate["key"]))
        
    if not scored:
        return None, f"{label}_Unresolved"
        
    # Select the candidate that agrees with the most populations, breaking ties by closest frequency
    max_agree = max(item[0] for item in scored)
    top = [item for item in scored if item[0] == max_agree]
    min_dist = min(item[1] for item in top)
    best = [item for item in top if abs(item[1] - min_dist) <= 1e-9]
    
    if len(best) != 1:
        return None, f"{label}_PopConflict"
        
    return best[0][2], f"{label}_Resolved_Global"