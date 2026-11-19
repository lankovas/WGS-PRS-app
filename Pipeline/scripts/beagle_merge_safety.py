"""Pure allele/dosage safety helpers used by merge_post_beagle.py.

This module provides tools to validate genetic alleles and safely remap 
imputed genotypes and dosages from Beagle back to their original reference 
orientations.
"""
from __future__ import annotations
import math
from typing import Any

# Standard representations of missing or symbolic alleles in VCF files
_PLACEHOLDERS = {"", ".", "*", "<*>", "<NON_REF>", "<NONREF>", "N"}


# ==============================================================================
# 1. ALLELE VALIDATION
# ==============================================================================

def is_concrete_allele(value: Any) -> bool:
    """Checks if an allele is a real DNA sequence (A, C, G, T) rather than a symbolic placeholder."""
    if value is None:
        return False
    allele = str(value).strip().upper()
    
    if allele in _PLACEHOLDERS:
        return False
    if allele.startswith("<") and allele.endswith(">"):
        return False
    if "[" in allele or "]" in allele:
        return False
        
    return bool(allele) and all(base in "ACGT" for base in allele)


def concrete_biallelic_pair(alleles: Any) -> tuple[str, str] | None:
    """Validates that a variant has exactly one valid REF and one valid ALT allele."""
    if not alleles or len(alleles) != 2:
        return None
        
    ref, alt = (str(a).strip().upper() for a in alleles)
    if not is_concrete_allele(ref) or not is_concrete_allele(alt) or ref == alt:
        return None
        
    return ref, alt


# ==============================================================================
# 2. TYPE CONVERSION & UTILITIES
# ==============================================================================

def scalar(value: Any) -> Any:
    """Extracts the first item if the input is a list or tuple; otherwise returns the value."""
    if isinstance(value, (tuple, list)):
        return value[0] if value else None
    return value


def finite_float(value: Any) -> float | None:
    """Safely converts a value to a float, returning None for invalid or infinite values."""
    value = scalar(value)
    if value is None:
        return None
        
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
        
    return number if math.isfinite(number) else None


def missing_gt(gt: Any) -> bool:
    """Checks if a Genotype (GT) field is empty or contains missing calls (e.g., './.')."""
    return gt is None or len(gt) == 0 or any(a is None for a in gt)


# ==============================================================================
# 3. DOSAGE & GENOTYPE REMAPPING
# ==============================================================================

def remap_beagle_call(
    prepared_pair: tuple[str, str],
    beagle_pair: tuple[str, str],
    gt: tuple[int | None, ...] | None,
    ds: float | None,
) -> tuple[tuple[int, ...], float] | None:
    """Maps Beagle's imputed Genotype (GT) and Dosage (DS) back to the original ALT allele.

    Beagle sometimes swaps the REF and ALT alleles during imputation based on 
    population frequencies. If a swap is detected, this function flips the genotype 
    (e.g., 0/0 becomes 1/1) and inverses the dosage (e.g., DS=1.8 becomes DS=0.2) 
    so it matches the original pipeline's orientation.
    """
    p_ref, p_alt = (str(a).strip().upper() for a in prepared_pair)
    b_ref, b_alt = (str(a).strip().upper() for a in beagle_pair)

    # Detect if Beagle inverted the allele pair
    if (b_ref, b_alt) == (p_ref, p_alt):
        swapped = False
    elif (b_ref, b_alt) == (p_alt, p_ref):
        swapped = True
    else:
        # Reject if Beagle introduced completely different alleles
        return None

    if missing_gt(gt) or any(int(a) not in (0, 1) for a in gt):
        return None

    gt = tuple(int(a) for a in gt)
    ploidy = len(gt)
    
    # Resolve dosage (calculate from GT if missing)
    ds = finite_float(ds)
    if ds is None:
        ds = float(sum(gt))
        
    if ds < -1e-6 or ds > ploidy + 1e-6:
        return None
        
    # Clamp dosage to strict boundaries [0, ploidy]
    ds = min(float(ploidy), max(0.0, ds))

    # Apply the mathematical flip if Beagle swapped REF and ALT
    if swapped:
        gt = tuple(1 - allele for allele in gt)
        ds = float(ploidy) - ds

    return gt, ds