#!/usr/bin/env python3
"""Allele-safe merger for prepared PRS VCFs and Beagle imputation output.

This script takes the missing genotypes that were sent to Beagle for imputation,
extracts Beagle's predictions, verifies that Beagle didn't alter the variant 
identity, and safely merges the imputed genotypes back into the main patient VCF.
Variants that Beagle failed to impute securely are routed to Allele Frequency fallback.
"""
from __future__ import annotations

import math
import os
import re
import sys
from collections import defaultdict
from typing import Any

import pysam

# Import local safety helpers to ensure alleles and dosages are correctly oriented
from beagle_merge_safety import (
    concrete_biallelic_pair,
    finite_float as _finite_float,
    missing_gt as _missing_gt,
    remap_beagle_call,
)


# ==============================================================================
# 1. HELPER FUNCTIONS
# ==============================================================================

def _normalise_chrom(chrom: Any) -> str:
    """Ensures all chromosome labels consistently use the 'chr' prefix."""
    text = str(chrom).strip()
    return text if text.lower().startswith("chr") else f"chr{text}"


def _record_pair(record: Any) -> tuple[str, str] | None:
    """Extracts a valid, concrete (REF, ALT) pair from a VCF record."""
    return concrete_biallelic_pair(getattr(record, "alleles", None))


def _safe_info_get(record: Any, key: str) -> Any:
    """Reads INFO fields safely without crashing if the field is missing from the VCF header."""
    try:
        header = getattr(record, "header", None)
        info_header = getattr(header, "info", None) if header is not None else None
        if info_header is not None and key not in info_header:
            return None
        info = getattr(record, "info", None)
        return info.get(key) if info is not None else None
    except (AttributeError, KeyError, TypeError, ValueError):
        return None


def _target_key(chrom: Any, pos: Any, ref: Any, alt: Any) -> tuple[str, int, str, str]:
    """Generates a unique identifier tuple for a specific variant."""
    return (
        _normalise_chrom(chrom),
        int(pos),
        str(ref).strip().upper(),
        str(alt).strip().upper(),
    )


# ==============================================================================
# 2. MAIN BEAGLE MERGER LOGIC
# ==============================================================================

class CanRiskBeagleMerger:
    """Orchestrates the safe merging of Beagle imputation results into the primary VCF."""

    def __init__(self, initial_vcf, beagle_vcf, out_vcf, dr2_threshold=0.9):
        self.initial_vcf = initial_vcf
        self.beagle_vcf = beagle_vcf
        self.out_vcf = out_vcf
        self.dr2_threshold = float(dr2_threshold)
        
        if not 0.0 <= self.dr2_threshold <= 1.0:
            raise ValueError("dr2_threshold must be between 0 and 1")

        self.targets_by_coord: dict[tuple[str, int], set[tuple[str, str]]] = defaultdict(set)
        self.beagle_data: dict[tuple[str, int, str, str], dict[str, Any]] = {}
        self.rejections: dict[tuple[str, int, str, str], str] = {}
        self.counts: dict[str, int] = defaultdict(int)

        # Execute Pipeline
        self.find_imputation_targets()
        self.extract_beagle_data()
        self.merge_and_write()

    def find_imputation_targets(self):
        """Scans the initial VCF to find all variants that require imputation (missing GT)."""
        with pysam.VariantFile(self.initial_vcf) as vcf:
            samples = list(vcf.header.samples)
            if len(samples) != 1:
                raise ValueError(
                    f"Expected one sample in prepared VCF; found {samples}"
                )
            self.sample = samples[0]

            for rec in vcf:
                gt = rec.samples[self.sample].get("GT", (None, None))
                if not _missing_gt(gt):
                    continue

                pair = _record_pair(rec)
                if pair is None:
                    self.counts["nonconcrete_or_multiallelic_target"] += 1
                    continue

                coord = (_normalise_chrom(rec.chrom), int(rec.pos))
                self.targets_by_coord[coord].add(pair)

        self.counts["target_pairs"] = sum(
            len(pairs) for pairs in self.targets_by_coord.values()
        )
        print(
            f"Beagle merger: {self.counts['target_pairs']} missing target allele pairs."
        )

    def extract_beagle_data(self):
        """Extracts and validates Beagle's imputed calls for the missing variants."""
        if not self.targets_by_coord:
            return

        if not os.path.isfile(self.beagle_vcf) or os.path.getsize(self.beagle_vcf) == 0:
            print("Notice: Beagle VCF is missing/empty; AF fallback will be used.")
            return

        candidates: dict[tuple[str, int, str, str], list[dict[str, Any]]] = defaultdict(list)

        try:
            with pysam.VariantFile(self.beagle_vcf) as bvcf:
                samples = list(bvcf.header.samples)
                if len(samples) != 1:
                    raise ValueError(
                        f"Expected one sample in Beagle VCF; found {samples}"
                    )
                bsample = samples[0]

                for rec in bvcf:
                    coord = (_normalise_chrom(rec.chrom), int(rec.pos))
                    target_pairs = self.targets_by_coord.get(coord)
                    if not target_pairs:
                        continue

                    self.counts["records_at_target_positions"] += 1
                    bpair = _record_pair(rec)
                    if bpair is None:
                        self.counts["unsupported_beagle_records"] += 1
                        continue

                    for pref, palt in target_pairs:
                        key = _target_key(coord[0], coord[1], pref, palt)
                        gt_raw = rec.samples[bsample].get("GT", (None, None))
                        ds_raw = _finite_float(rec.samples[bsample].get("DS"))

                        # Verify that Beagle didn't mutate the allele pair and correct orientation if needed
                        mapped = remap_beagle_call((pref, palt), bpair, gt_raw, ds_raw)
                        if mapped is None:
                            if bpair not in ((pref, palt), (palt, pref)):
                                self.counts["allele_pair_mismatches"] += 1
                            else:
                                self.rejections[key] = "Beagle_MissingOrInvalidGT_DS"
                            continue

                        mapped_gt, mapped_ds = mapped
                        
                        # Validate Beagle's imputation confidence (DR2)
                        dr2 = _finite_float(_safe_info_get(rec, "DR2"))
                        if dr2 is None:
                            dr2 = 0.0

                        ploidy = len(mapped_gt)
                        high_conf = (
                            dr2 >= self.dr2_threshold
                            or mapped_ds < 0.01
                            or mapped_ds > ploidy - 0.01
                        )
                        if not high_conf:
                            self.rejections[key] = "Beagle_LowDR2"
                            self.counts["low_confidence"] += 1
                            continue

                        candidates[key].append({
                            "GT": mapped_gt,
                            "DS": mapped_ds,
                            "DR2": dr2,
                            "source_pair": bpair,
                        })

        except (OSError, ValueError, TypeError) as exc:
            print(
                f"WARNING: Could not parse Beagle VCF ({exc}); AF fallback will be used."
            )
            return

        for key, values in candidates.items():
            if len(values) != 1:
                self.rejections[key] = "Beagle_DuplicateMatchingRecords"
                self.counts["duplicate_matching_records"] += len(values)
                continue
            self.beagle_data[key] = values[0]

        print(
            "Beagle candidate audit: "
            f"accepted={len(self.beagle_data)}, "
            f"allele_pair_mismatches={self.counts['allele_pair_mismatches']}, "
            f"low_confidence={self.counts['low_confidence']}, "
            f"duplicates={self.counts['duplicate_matching_records']}"
        )

    def merge_and_write(self):
        """Writes the final output VCF, merging Beagle predictions or defaulting to AF imputation."""
        with pysam.VariantFile(self.initial_vcf) as invcf:
            header = invcf.header
            if "END" not in header.info:
                header.info.add("END", 1, "Integer", "Stop position of the interval")
            samples = list(header.samples)
            if len(samples) != 1:
                raise ValueError(f"Expected one sample in prepared VCF; found {samples}")
            sample = samples[0]

            with pysam.VariantFile(self.out_vcf, "w", header=header) as out:
                for rec in invcf:
                    if "END" in rec.info:
                        try:
                            rec.info["END"] = int(rec.info["END"])
                        except (ValueError, TypeError):
                            del rec.info["END"]

                    gt_current = rec.samples[sample].get("GT", (None, None))
                    coord = (_normalise_chrom(rec.chrom), int(rec.pos))

                    # Process missing genotypes that were flagged for imputation
                    if coord in self.targets_by_coord and _missing_gt(gt_current):
                        pair = _record_pair(rec)
                        
                        # Clean existing notes
                        notes_raw = rec.samples[sample].get("NOTE")
                        if isinstance(notes_raw, (tuple, list)):
                            notes = [str(n) for n in notes_raw if n not in (None, "", ".")]
                        else:
                            notes = [
                                n for n in str(notes_raw or "").split(",")
                                if n and n != "."
                            ]
                        notes = [n for n in notes if n != "Target_Imputation"]

                        key = (
                            _target_key(coord[0], coord[1], pair[0], pair[1])
                            if pair is not None else None
                        )
                        data = self.beagle_data.get(key) if key is not None else None

                        # Assign Beagle Imputed Data
                        if data is not None:
                            try:
                                rec.samples[sample]["GT"] = data["GT"]
                                rec.samples[sample]["DS"] = float(data["DS"])
                                notes.append("Imputed_Beagle")
                                self.counts["merged_beagle"] += 1
                            except (ValueError, TypeError, KeyError) as exc:
                                # Reset fields and fallback to AF if assignment fails
                                rec.samples[sample]["GT"] = (None, None)
                                rec.samples[sample]["DS"] = None
                                notes.append("Beagle_MergeAssignmentFailed")
                                notes.append("AF_Imputed")
                                self.counts["assignment_failures"] += 1
                                print(f"WARNING: Merge assignment failed at {coord}: {exc}")
                        
                        # Fallback to Allele Frequency if Beagle was rejected or unavailable
                        else:
                            rec.samples[sample]["GT"] = (None, None)
                            rec.samples[sample]["DS"] = None
                            notes.append("AF_Imputed")
                            if key is not None and key in self.rejections:
                                notes.append(self.rejections[key])
                            elif pair is None:
                                notes.append("Beagle_TargetNotBiallelic")
                            else:
                                notes.append("Beagle_NoMatchingAllelesAtTarget")
                            self.counts["af_fallback"] += 1

                        rec.samples[sample]["NOTE"] = ",".join(dict.fromkeys(notes)) if notes else None

                    out.write(rec)

        # Standardize formatting descriptions to prevent parsing warnings in some downstream tools
        with open(self.out_vcf, "r", encoding="utf-8") as src:
            lines = src.readlines()
        with open(self.out_vcf, "w", encoding="utf-8") as dst:
            for line in lines:
                if line.startswith("##FORMAT=<ID=DS"):
                    line = re.sub(r"Description=Dosage", 'Description="Dosage"', line)
                elif line.startswith("##FORMAT=<ID=GT"):
                    line = re.sub(r"Description=Genotype", 'Description="Genotype"', line)
                elif line.startswith("##FORMAT=<ID=NOTE"):
                    line = re.sub(
                        r'Description=Variant Notes(?: and Filters)?',
                        'Description="Variant Notes"',
                        line,
                    )
                dst.write(line)

        print(
            "Beagle merge complete: "
            f"merged={self.counts['merged_beagle']}, "
            f"AF_fallback={self.counts['af_fallback']}, "
            f"assignment_failures={self.counts['assignment_failures']}"
        )


# ==============================================================================
# 3. CLI ENTRY POINT
# ==============================================================================

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print(
            "Usage: python merge_post_beagle.py "
            "<initial.vcf> <beagle.vcf[.gz]> <out.vcf> [dr2_threshold]"
        )
        raise SystemExit(1)

    try:
        threshold = float(sys.argv[4]) if len(sys.argv) > 4 else 0.9
    except ValueError:
        threshold = 0.9

    CanRiskBeagleMerger(sys.argv[1], sys.argv[2], sys.argv[3], threshold)