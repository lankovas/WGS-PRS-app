#!/usr/bin/env python3
"""Build a harmonized PRS reference cache and ancestry-specific raw-score stats.

This script processes a Polygenic Risk Score (PRS) panel against a reference
genome database (PLINK 2 pgen/pvar/psam files) to resolve ambiguous alleles,
match variants correctly, and compute population-specific baseline scores.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import os
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

# Load pipeline-specific helper utilities
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pgs_harmonise import (  # noqa: E402
    clean_chrom,
    is_palindromic,
    match_variants,
    resolve_candidates_global,
    reverse_complement,
)

# Standard representations of missing alleles
MISSING_ALLELES = {"", "N", ".", "NAN", "NONE", "<NA>"}
CACHE_SCHEMA_VERSION = 2


# ==============================================================================
# 1. UTILITY & CLEANUP FUNCTIONS
# ==============================================================================

def log(message: str) -> None:
    """Print timestamped log messages to standard output."""
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {message}", flush=True)


def clean_pos(value: Any) -> str:
    """Safely converts a coordinate value to an integer string."""
    try:
        if pd.isna(value):
            return ""
        number = int(float(value))
        return str(number) if number > 0 else ""
    except (TypeError, ValueError, OverflowError):
        return ""


def clean_allele(value: Any) -> str:
    """Standardizes allele strings to uppercase and handles missing values."""
    if value is None or pd.isna(value):
        return ""
    value = str(value).strip().upper()
    return "" if value in MISSING_ALLELES else value


def check_bool(value: Any) -> bool:
    """Safely parses a variety of boolean string representations."""
    if value is None or pd.isna(value):
        return False
    return str(value).strip().upper() in {"TRUE", "T", "1", "YES", "Y"}


def _split_comma_field(value: Any) -> list[str]:
    """Splits a comma-separated string into a list of cleaned components."""
    if value is None or pd.isna(value):
        return []
    return [part.strip() for part in str(value).split(",")]


# ==============================================================================
# 2. FILE LOADING & PARSING
# ==============================================================================

def load_variant_actions(actions_file: str) -> dict[str, str]:
    """Loads manual variant overrides (e.g. Keep, Drop) configured by the user."""
    actions: dict[str, str] = {}
    if not actions_file or not os.path.exists(actions_file) or os.path.getsize(actions_file) == 0:
        return actions
    try:
        df = pd.read_csv(actions_file, sep="\t", dtype=str)
        if "Action" not in df.columns:
            return actions
            
        chr_col = next((c for c in df.columns if c.lower() in {"chr", "chr_name", "hm_chr"}), None)
        pos_col = next((c for c in df.columns if c.lower() in {"pos", "hm_pos", "chr_position"}), None)
        
        if not chr_col or not pos_col:
            return actions
            
        for _, row in df.iterrows():
            chrom, pos = clean_chrom(row.get(chr_col, "")), clean_pos(row.get(pos_col, ""))
            if chrom and pos:
                actions[f"{chrom}:{pos}"] = str(row["Action"]).strip().lower()
    except Exception as exc:
        log(f"WARNING: Could not read variant-actions file ({exc}); no manual actions applied.")
    return actions


def read_panel(path: str) -> pd.DataFrame:
    """Reads a PGS Catalog scoring file, automatically bypassing leading comments."""
    open_func = gzip.open if str(path).lower().endswith(".gz") else open
    skip_lines = 0
    with open_func(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                skip_lines += 1
            else:
                break
    try:
        df = pd.read_csv(path, sep="\t", skiprows=skip_lines, dtype=str, compression="infer")
        if len(df.columns) <= 1:
            raise ValueError("tab-delimited parse produced one column")
    except Exception:
        df = pd.read_csv(path, sep=",", skiprows=skip_lines, dtype=str, compression="infer")
        
    df.columns = [str(column).strip().lower().replace(" ", "_") for column in df.columns]
    return df


def pick_col(df: pd.DataFrame, choices: list[str]) -> str | None:
    """Returns the first column name from a list of choices that exists in the dataframe."""
    return next((column for column in choices if column in df.columns), None)


def _read_psam(psam_path: str) -> tuple[pd.DataFrame, str, str]:
    """Reads the PLINK sample information file to extract population/ancestry labels."""
    if not os.path.exists(psam_path):
        raise FileNotFoundError(f"Could not find reference sample file: {psam_path}")
        
    psam = pd.read_csv(psam_path, sep="\t", dtype=str)
    psam.columns = [str(c).lstrip("#") for c in psam.columns]
    
    iid_col = "IID" if "IID" in psam.columns else (psam.columns[0] if len(psam.columns) else None)
    if not iid_col:
        raise ValueError(f"No IID column found in {psam_path}")
        
    pop_col = next((c for c in psam.columns if c.lower() in {"superpop", "super_pop", "superpopulation"}), None)
    if pop_col is None:
        pop_col = next((c for c in psam.columns if c.lower() in {"pop", "population"} or "superpop" in c.lower()), None)
    if pop_col is None:
        raise ValueError(f"Could not find a superpopulation column in {psam_path}; columns: {list(psam.columns)}")
    if psam[iid_col].duplicated().any():
        raise ValueError(f"Duplicate IID values found in {psam_path}")
        
    return psam, iid_col, pop_col


# ==============================================================================
# 3. PLINK DATABASE & CACHING
# ==============================================================================

def load_plink_database(pvar_path: str, afreq_prefix: str, super_pops: list[str]):
    """Loads PVAR variants and their associated population allele frequencies into memory."""
    records_by_pos: dict[str, list[dict[str, Any]]] = {}
    records_by_rsid: dict[str, list[dict[str, Any]]] = {}
    af_rows: dict[str, dict[str, Any]] = {}

    # Load Variants
    if os.path.exists(pvar_path):
        with open(pvar_path, "rt", encoding="utf-8") as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 5:
                    continue
                chrom, pos, native_id = clean_chrom(parts[0]), clean_pos(parts[1]), parts[2].strip()
                ref = clean_allele(parts[3])
                alts = [clean_allele(a) for a in _split_comma_field(parts[4]) if a]
                
                if not chrom or not pos or not ref or not alts:
                    continue
                    
                record = {
                    "native_id": native_id, "chrom": chrom, "pos": pos,
                    "ref": ref, "alts": alts,
                }
                records_by_pos.setdefault(f"{chrom}:{pos}", []).append(record)
                if native_id.lower().startswith("rs"):
                    records_by_rsid.setdefault(native_id, []).append(record)

    # Load Population Frequencies
    for pop in super_pops:
        af_path = f"{afreq_prefix}_{pop}.afreq"
        if not os.path.exists(af_path):
            log(f"WARNING: Missing frequency file {af_path}; AF-based ambiguity resolution may drop variants.")
            continue
            
        try:
            af_df = pd.read_csv(af_path, sep="\t", dtype=str)
        except Exception as exc:
            log(f"WARNING: Could not parse {af_path}: {exc}")
            continue
            
        required = {"ID", "REF", "ALT", "ALT_FREQS"}
        if not required.issubset(af_df.columns):
            log(f"WARNING: {af_path} is missing columns {sorted(required - set(af_df.columns))}.")
            continue

        for _, row in af_df.iterrows():
            native_id = str(row["ID"]).strip()
            ref = clean_allele(row["REF"])
            alts = [clean_allele(a) for a in _split_comma_field(row["ALT"])]
            freq_values = _split_comma_field(row["ALT_FREQS"])
            
            if not native_id or not ref or not alts:
                continue
                
            entry = af_rows.setdefault(native_id, {"ref": ref, "alt_freqs": {}})
            
            if entry["ref"] != ref:
                log(f"WARNING: Conflicting REF alleles for ID {native_id} in {af_path}; ignoring this frequency row.")
                continue
                
            for alt, raw_freq in zip(alts, freq_values):
                try:
                    frequency = float(raw_freq)
                except (TypeError, ValueError):
                    frequency = float("nan")
                if not math.isfinite(frequency) or not 0.0 <= frequency <= 1.0:
                    continue
                entry["alt_freqs"].setdefault(alt, {})[pop] = frequency

    return records_by_pos, records_by_rsid, af_rows


def freq_of(record: dict[str, Any], idx: int, pop: str, af_rows: dict[str, Any]) -> float | None:
    """Calculates the specific population frequency of a target allele."""
    entry = af_rows.get(record["native_id"])
    if not entry:
        return None
        
    if idx > 0:
        alt = record["alts"][idx - 1]
        value = entry["alt_freqs"].get(alt, {}).get(pop)
        return float(value) if value is not None and math.isfinite(float(value)) else None

    # If the target is the Reference allele, calculate it as 1 - (sum of all ALT frequencies)
    alt_freqs = []
    for alt in record["alts"]:
        value = entry["alt_freqs"].get(alt, {}).get(pop)
        if value is None or not math.isfinite(float(value)):
            return None
        alt_freqs.append(float(value))
        
    ref_af = 1.0 - sum(alt_freqs)
    if ref_af < -1e-6 or ref_af > 1.0 + 1e-6:
        return None
    return min(1.0, max(0.0, ref_af))


def _remove_prefix_outputs(prefix: str) -> None:
    """Cleans up PLINK output files matching a specific prefix."""
    for suffix in (".pgen", ".pvar", ".psam", ".log", ".nosex", ".afreq", ".acount"):
        path = prefix + suffix
        if os.path.exists(path):
            os.remove(path)


def _make_coordinate_bed(df: pd.DataFrame, position_column: str, destination: str) -> bool:
    """Creates a BED file of variant coordinates to instruct PLINK which regions to extract."""
    points = set()
    for _, row in df.iterrows():
        chrom = clean_chrom(row.get("chr_clean", ""))
        pos = clean_pos(row.get(position_column, ""))
        if not chrom or not pos:
            continue
        # Support both standard and "chr" prefixed notations
        points.add((chrom, pos, pos))
        points.add((f"chr{chrom}", pos, pos))
        
    if not points:
        return False
    pd.DataFrame(sorted(points), columns=["CHR", "START", "END"]).to_csv(
        destination, sep="\t", index=False, header=False
    )
    return True


def _stat_fingerprint(path: str) -> dict[str, Any]:
    """Generates a lightweight fingerprint of a file using size and modified time."""
    absolute = os.path.abspath(path)
    try:
        stat = os.stat(absolute)
    except OSError:
        return {"path": absolute, "exists": False}
    return {
        "path": absolute,
        "exists": True,
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
    }


def _sha256_file(path: str) -> str:
    """Computes a SHA256 hash of a file for exact caching validation."""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _cache_is_complete(subset_prefix: str, super_pops: list[str]) -> bool:
    """Verifies that all required files for a PLINK subset cache exist."""
    required = [f"{subset_prefix}{suffix}" for suffix in (".pgen", ".pvar", ".psam")]
    required.extend(f"{subset_prefix}_{pop}.afreq" for pop in super_pops)
    return all(os.path.isfile(path) and os.path.getsize(path) > 0 for path in required)


def _build_subset_and_frequencies(
    input_prefix: str,
    subset_prefix: str,
    position_df: pd.DataFrame,
    position_column: str,
    base_cache_dir: str,
    psam_df: pd.DataFrame,
    pop_col: str,
    iid_col: str,
    super_pops: list[str],
    *,
    force_rebuild: bool = False,
) -> bool:
    """Generates or reuses a trimmed PLINK database containing only variants relevant to the panel.
    
    This drastically speeds up downstream matching and scoring by discarding the rest of the genome.
    """
    os.makedirs(base_cache_dir, exist_ok=True)
    bed_path = os.path.join(base_cache_dir, Path(subset_prefix).name + "_coords.bed")
    manifest_path = subset_prefix + ".cache.json"
    
    if not _make_coordinate_bed(position_df, position_column, bed_path):
        log(f"WARNING: No usable coordinates for {subset_prefix}; skipping extraction.")
        try:
            os.remove(manifest_path)
        except FileNotFoundError:
            pass
        return False

    # Extract sample IDs specific to each population for frequency calculation
    keep_paths: dict[str, str] = {}
    for pop in super_pops:
        keep_path = os.path.join(base_cache_dir, f"{pop}.keep")
        pop_rows = (
            psam_df.loc[psam_df[pop_col].astype(str) == str(pop), iid_col]
            .dropna().astype(str)
        )
        pop_rows.to_csv(keep_path, sep="\t", index=False, header=False)
        keep_paths[pop] = keep_path

    # Construct the unique signature to detect if the cache is stale
    signature = {
        "schema_version": CACHE_SCHEMA_VERSION,
        "input_prefix": os.path.abspath(input_prefix),
        "source_files": {
            suffix: _stat_fingerprint(input_prefix + suffix)
            for suffix in (".pgen", ".pvar", ".psam")
        },
        "coordinates_sha256": _sha256_file(bed_path),
        "population_keep_sha256": {
            pop: _sha256_file(keep_paths[pop]) for pop in super_pops
        },
        "super_pops": sorted(str(pop) for pop in super_pops),
        "operation": "plink2 --pfile --extract bed1 --make-pgen; per-pop --keep --freq",
    }

    if not force_rebuild and _cache_is_complete(subset_prefix, super_pops):
        try:
            with open(manifest_path, "rt", encoding="utf-8") as handle:
                cached_manifest = json.load(handle)
            if cached_manifest.get("signature") == signature:
                log(f"Reusing validated reference cache: {subset_prefix}")
                return True
            log(f"Reference cache inputs changed; rebuilding {subset_prefix}.")
        except (OSError, json.JSONDecodeError, AttributeError):
            log(f"No valid cache manifest for {subset_prefix}; rebuilding.")
    elif force_rebuild:
        log(f"Forced cache rebuild requested for {subset_prefix}.")
    else:
        log(f"Reference cache incomplete for {subset_prefix}; rebuilding.")

    # Execute PLINK2 to perform the data extraction
    try:
        os.remove(manifest_path)
    except FileNotFoundError:
        pass
    _remove_prefix_outputs(subset_prefix)
    
    for pop in super_pops:
        for suffix in (".afreq", ".log", ".acount"):
            old = f"{subset_prefix}_{pop}{suffix}"
            if os.path.exists(old):
                os.remove(old)

    cmd = [
        "plink2", "--pfile", input_prefix,
        "--extract", "bed1", bed_path,
        "--make-pgen", "--out", subset_prefix,
        "--threads", "8", "--memory", "16000",
    ]
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
    except subprocess.CalledProcessError as exc:
        if exc.returncode == 13:
            log(f"WARNING: No variants extracted for {input_prefix}; source will contribute no variants.")
            _remove_prefix_outputs(subset_prefix)
            return False
        raise

    # Generate Allele Frequencies per population
    for pop in super_pops:
        freq_prefix = f"{subset_prefix}_{pop}"
        subprocess.run(
            ["plink2", "--pfile", subset_prefix, "--keep", keep_paths[pop],
             "--freq", "--out", freq_prefix, "--threads", "8", "--memory", "16000"],
            check=True, stdout=subprocess.DEVNULL,
        )

    if not _cache_is_complete(subset_prefix, super_pops):
        raise RuntimeError(
            f"PLINK2 completed but cache outputs are incomplete for {subset_prefix}; "
            "the cache manifest was not written."
        )

    # Save the manifest proving the cache is valid
    tmp_manifest = manifest_path + ".tmp"
    with open(tmp_manifest, "wt", encoding="utf-8") as handle:
        json.dump({"signature": signature, "created_at": datetime.now().isoformat()}, handle, indent=2)
        handle.write("\n")
    os.replace(tmp_manifest, manifest_path)
    log(f"Built and validated reference cache: {subset_prefix}")
    
    return True


# ==============================================================================
# 4. ALLELE HARMONIZATION
# ==============================================================================

def _dedup_records(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Removes exact duplicate variants from a list of records."""
    seen = set()
    output = []
    for record in records:
        key = (record["native_id"], record["chrom"], record["pos"], record["ref"], tuple(record["alts"]))
        if key not in seen:
            seen.add(key)
            output.append(record)
    return output


def _record_matches(
    record: dict[str, Any],
    effect: str,
    other: str,
    *,
    assume_forward: bool,
    effect_is_alt: bool,
) -> list[dict[str, Any]]:
    """Evaluates how well a panel variant matches a reference record.
    
    Match Strengths:
      0 = Exact allele-pair match
      1 = Effect-only match (because panel lacked an 'other' allele)
      2 = Effect-only fallback (a provided 'other' allele was wrong, but the effect allele matched)
    """
    allele_set = [record["ref"]] + record["alts"]
    candidates = []

    if effect_is_alt:
        for idx, alt in enumerate(record["alts"], start=1):
            if effect != alt:
                continue
            if not other: strength = 1
            elif other == record["ref"]: strength = 0
            else: strength = 2
            
            candidates.append({
                "native_id": record["native_id"], "rec": record, "idx": idx,
                "flipped": False, "strength": strength,
            })
        return candidates

    matches, _ambiguous, relaxed = match_variants(
        effect, other, allele_set, allow_complement=not assume_forward
    )
    
    if assume_forward:
        matches = [(idx, flipped) for idx, flipped in matches if not flipped]
        
    if not matches:
        return []
        
    if other and not relaxed: strength = 0
    elif not other: strength = 1
    else: strength = 2

    for idx, flipped in matches:
        candidates.append({
            "native_id": record["native_id"], "rec": record, "idx": idx,
            "flipped": flipped, "strength": strength,
        })
    return candidates


def _candidate_key(candidate: dict[str, Any], source: str) -> tuple[Any, ...]:
    """Generates a unique identifier for a specific matched variant."""
    rec = candidate["rec"]
    return (
        source, candidate["native_id"], candidate["idx"], candidate["flipped"],
        rec["chrom"], rec["pos"], rec["ref"], tuple(rec["alts"]),
    )


def _dedup_candidates(candidates: list[dict[str, Any]], source: str) -> list[dict[str, Any]]:
    """Keeps the strongest match when multiple candidates represent the exact same variant."""
    by_key: dict[tuple[Any, ...], dict[str, Any]] = {}
    for candidate in candidates:
        key = _candidate_key(candidate, source)
        if key not in by_key or candidate["strength"] < by_key[key]["strength"]:
            candidate["source"] = source
            by_key[key] = candidate
    return list(by_key.values())


def _candidate_population_afs(candidate: dict[str, Any], af_rows: dict[str, Any], super_pops: list[str]) -> dict[str, float | None]:
    """Fetches population allele frequencies for a specific candidate match."""
    return {
        pop: freq_of(candidate["rec"], candidate["idx"], pop, af_rows)
        for pop in super_pops
    }


def _format_audit_row(
    row: pd.Series,
    *,
    effect: str,
    other: str,
    eaf: float | None,
    palindromic: bool,
    chosen: dict[str, Any] | None,
    final_effect_allele: str,
    status: str,
    note: str,
    requested_build: str,
) -> dict[str, Any]:
    """Formats a single row for the harmonization audit log."""
    rec = chosen["rec"] if chosen else {}
    is_hg38 = str(requested_build).lower() in {"hg38", "grch38", "38"}
    matched_build = ("GRCh38" if is_hg38 else "GRCh37") if chosen and chosen["source"] == "primary" else (("GRCh37" if is_hg38 else "GRCh38") if chosen else "")
    
    return {
        "chr": row.get("chr_clean", ""),
        "pos_hg19": row.get("active_pos_hg19", ""),
        "pos_hg38": row.get("active_pos_hg38", ""),
        "rsID": row.get("rsid_value", ""),
        "effect_allele": effect,
        "other_allele": other,
        "scoring_file_EAF": "" if eaf is None else eaf,
        "palindromic": bool(palindromic),
        "panel_native_id": chosen["native_id"] if chosen else "",
        "harmonized_effect_allele": final_effect_allele,
        "weight_inverted": False,
        "flipped": bool(chosen["flipped"]) if chosen else "",
        "source": chosen["source"] if chosen else "",
        "target_build": "GRCh38" if is_hg38 else "GRCh37",
        "matched_build": matched_build,
        "matched_chrom": rec.get("chrom", "") if chosen else "",
        "matched_pos": rec.get("pos", "") if chosen else "",
        "matched_ref": rec.get("ref", "") if chosen else "",
        "matched_alts": ",".join(rec.get("alts", [])) if chosen else "",
        "match_strength": chosen.get("strength", "") if chosen else "",
        "matched_idx": chosen.get("idx", "") if chosen else "",
        "matched_effect_allele": final_effect_allele,
        "target_effect_allele": effect,
        "target_other_allele": other,
        "is_additive": bool(row.get("is_additive", True)),
        "status": status,
        "note": note,
    }


# ==============================================================================
# 5. MAIN EXECUTION PIPELINE
# ==============================================================================

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build", required=True, help="Build of kg_primary: hg19/GRCh37 or hg38/GRCh38")
    parser.add_argument("--kg_primary", required=True, help="PLINK 2 pfile prefix for the primary build")
    parser.add_argument("--kg_fallback", required=True, help="PLINK 2 pfile prefix for the alternate build")
    parser.add_argument("--prs_file", required=True, help="PGS Catalog scoring file (.txt or .txt.gz)")
    parser.add_argument("--out_dir", required=True)
    parser.add_argument("--base_cache_dir", required=True,
                        help="Shared cache directory reused across harmonization modes")
    parser.add_argument("--rebuild_cache", action="store_true",
                        help="Force rebuilding coordinate subsets/frequencies (otherwise reuse validated cache)")
    parser.add_argument("--variant_actions", default="")
    parser.add_argument("--af_cache_out", required=True)
    parser.add_argument("--assume_forward", type=str, default="False",
                        help="Assert alleles are already in forward/reference orientation")
    parser.add_argument("--effect_is_alt", type=str, default="False",
                        help="Assert effect allele is a literal ALT allele in the forward/reference representation")
    parser.add_argument("--pal_low", type=float, default=0.40)
    parser.add_argument("--pal_high", type=float, default=0.60)
    parser.add_argument("--eaf_pop", default="EUR",
                        help="Population corresponding best to the scoring-file EAF; default EUR")
    parser.add_argument("--af_mode", default="1000G", help="Retained for pipeline compatibility; currently informational")
    args = parser.parse_args()

    if not 0.0 <= args.pal_low < args.pal_high <= 1.0:
        parser.error("Require 0 <= --pal_low < --pal_high <= 1")
        
    assume_forward = str(args.assume_forward).strip().lower() in {"true", "1", "yes"}
    effect_is_alt = str(args.effect_is_alt).strip().lower() in {"true", "1", "yes"}
    if assume_forward and effect_is_alt:
        parser.error("--assume_forward and --effect_is_alt describe different modes; set at most one to True")

    t0 = time.time()
    os.makedirs(args.base_cache_dir, exist_ok=True)
    os.makedirs(args.out_dir, exist_ok=True)
    actions_dict = load_variant_actions(args.variant_actions)

    # -------------------------------------------------------------------------
    # STEP 1: LOAD & NORMALIZE PGS PANEL
    # -------------------------------------------------------------------------
    log("STEP 1: Reading the PGS panel...")
    df = read_panel(args.prs_file)
    if df.empty:
        raise ValueError(f"No variants read from scoring file: {args.prs_file}")

    chr_col = pick_col(df, ["hm_chr", "chr_name", "chr"])
    hg19_pos_col = pick_col(df, ["chr_position", "pos", "hm_pos"])
    hg38_pos_col = pick_col(df, ["hm_pos"])
    requested_is_hg38 = str(args.build).lower() in {"hg38", "grch38", "38"}
    
    if requested_is_hg38:
        effect_col = pick_col(df, ["hm_effect_allele", "effect_allele", "alt"])
        other_col = pick_col(df, ["hm_other_allele", "other_allele"])
    else:
        effect_col = pick_col(df, ["effect_allele", "alt"])
        other_col = pick_col(df, ["other_allele"])
        
    eaf_col = pick_col(df, ["allelefrequency_effect", "hm_effect_allele_frequency", "effect_allele_frequency"])
    rsid_col = pick_col(df, ["hm_rsid", "rsid"])

    if not chr_col: raise ValueError(f"No chromosome column found in panel.")
    if not hg19_pos_col and not hg38_pos_col: raise ValueError(f"No coordinate column found in panel.")
    if not effect_col: raise ValueError(f"No effect_allele/ALT column found in panel.")

    # Calculate additive effect weights natively from OR/HR if necessary
    if "effect_weight" not in df.columns:
        if "or" in df.columns:
            df["effect_weight"] = np.log(pd.to_numeric(df["or"], errors="coerce"))
        elif "hr" in df.columns:
            df["effect_weight"] = np.log(pd.to_numeric(df["hr"], errors="coerce"))
        else:
            raise ValueError("Panel has no effect_weight, OR, or HR column")
            
    df["effect_weight"] = pd.to_numeric(df["effect_weight"], errors="coerce")
    if df["effect_weight"].isna().any():
        bad = int(df["effect_weight"].isna().sum())
        raise ValueError(f"{bad} panel row(s) have missing/non-numeric effect weights")

    chr_candidates = [c for c in ["hm_chr", "chr_name", "chr"] if c in df.columns]
    def row_chrom(row: pd.Series) -> str:
        for col in chr_candidates:
            chrom = clean_chrom(row.get(col, ""))
            if chrom: return chrom
        return ""
    df["chr_clean"] = df.apply(row_chrom, axis=1)

    if "chr_position" in df.columns: df["active_pos_hg19"] = df["chr_position"].map(clean_pos)
    elif "pos" in df.columns: df["active_pos_hg19"] = df["pos"].map(clean_pos)
    elif hg38_pos_col: df["active_pos_hg19"] = df[hg38_pos_col].map(clean_pos)
    else: df["active_pos_hg19"] = ""

    if "hm_pos" in df.columns: df["active_pos_hg38"] = df["hm_pos"].map(clean_pos)
    else: df["active_pos_hg38"] = df["active_pos_hg19"]

    is_hg38 = str(args.build).lower() in {"hg38", "grch38", "38"}
    df["primary_pos"] = df["active_pos_hg38"] if is_hg38 else df["active_pos_hg19"]
    df["fallback_pos"] = df["active_pos_hg19"] if is_hg38 else df["active_pos_hg38"]
    df["effect_clean"] = df[effect_col].map(clean_allele)
    
    if requested_is_hg38 and "effect_allele" in df.columns:
        original_effect = df["effect_allele"].map(clean_allele)
        df["effect_clean"] = df["effect_clean"].where(df["effect_clean"] != "", original_effect)
        
    if other_col:
        other_values = df[other_col].map(clean_allele)
    else:
        other_values = pd.Series([""] * len(df), index=df.index)
        
    if requested_is_hg38 and "other_allele" in df.columns:
        original_other = df["other_allele"].map(clean_allele)
        other_values = other_values.where(other_values != "", original_other)
        
    if requested_is_hg38 and "hm_inferotherallele" in df.columns:
        inferred = df["hm_inferotherallele"].map(clean_allele)
        df["other_clean"] = inferred.where(inferred != "", other_values)
    else:
        df["other_clean"] = other_values
        
    df["other_clean"] = df["other_clean"].str.split("/").str[0].map(clean_allele)
    df["eaf_clean"] = pd.to_numeric(df[eaf_col], errors="coerce") if eaf_col else np.nan
    df["rsid_value"] = df[rsid_col].fillna("").astype(str).str.strip() if rsid_col else ""
    df["is_dominant"] = df.get("is_dominant", pd.Series([False] * len(df), index=df.index)).apply(check_bool)
    df["is_recessive"] = df.get("is_recessive", pd.Series([False] * len(df), index=df.index)).apply(check_bool)
    
    dosage_cols = ["dosage_0_weight", "dosage_1_weight", "dosage_2_weight"]
    if all(column in df.columns for column in dosage_cols):
        df["has_dosage_weights"] = df[dosage_cols].notna().all(axis=1)
    else:
        df["has_dosage_weights"] = False
    df["is_additive"] = ~(df["is_dominant"] | df["is_recessive"] | df["has_dosage_weights"])

    df = df.loc[(df["primary_pos"] != "") | (df["fallback_pos"] != "")].copy()
    if df.empty:
        raise ValueError("No panel variants have a usable coordinate")

    primary_prefix = os.path.join(args.base_cache_dir, "panel_subset_primary")
    fallback_prefix = os.path.join(args.base_cache_dir, "panel_subset_fallback")
    primary_psam, iid_col, pop_col = _read_psam(f"{args.kg_primary}.psam")
    primary_psam = primary_psam.rename(columns={iid_col: "IID", pop_col: "SuperPop"})
    iid_col, pop_col = "IID", "SuperPop"
    
    super_pops = sorted(primary_psam[pop_col].dropna().astype(str).unique().tolist())
    if not super_pops:
        raise ValueError("No superpopulation labels found in primary reference PSAM")
    if args.eaf_pop not in super_pops:
        log(f"WARNING: --eaf_pop {args.eaf_pop} is not in the primary PSAM populations {super_pops}; global AF fallback will be used where possible.")


    # -------------------------------------------------------------------------
    # STEP 2: BUILD PLINK SUBSETS
    # -------------------------------------------------------------------------
    log("STEP 2: Checking shared coordinate-subset caches and population allele frequencies...")
    primary_ok = _build_subset_and_frequencies(
        args.kg_primary, primary_prefix, df, "primary_pos", args.base_cache_dir,
        primary_psam, pop_col, iid_col, super_pops,
        force_rebuild=args.rebuild_cache,
    )
    fallback_ok = _build_subset_and_frequencies(
        args.kg_fallback, fallback_prefix, df, "fallback_pos", args.base_cache_dir,
        primary_psam, pop_col, iid_col, super_pops,
        force_rebuild=args.rebuild_cache,
    )

    log("STEP 3: Loading reference records and allele frequencies...")
    p_pos, p_rsid, p_af = load_plink_database(f"{primary_prefix}.pvar", primary_prefix, super_pops) if primary_ok else ({}, {}, {})
    f_pos, f_rsid, f_af = load_plink_database(f"{fallback_prefix}.pvar", fallback_prefix, super_pops) if fallback_ok else ({}, {}, {})

    dbs = [
        {"source": "primary", "chrom": None, "by_pos": p_pos, "by_rsid": p_rsid, "af": p_af},
        {"source": "fallback", "chrom": None, "by_pos": f_pos, "by_rsid": f_rsid, "af": f_af},
    ]


    # -------------------------------------------------------------------------
    # STEP 4: HARMONIZE VARIANTS
    # -------------------------------------------------------------------------
    log(f"STEP 4: Harmonizing {len(df)} panel rows...")
    audit_lines: list[dict[str, Any]] = []
    mapped_rows: list[dict[str, Any]] = []
    drop_reasons: dict[str, int] = {}
    kept = 0

    for _, row in df.iterrows():
        chrom = row["chr_clean"]
        primary_pos = clean_pos(row["primary_pos"])
        fallback_pos = clean_pos(row["fallback_pos"])
        primary_key = f"{chrom}:{primary_pos}" if primary_pos else ""
        fallback_key = f"{chrom}:{fallback_pos}" if fallback_pos else ""
        
        effect, other = row["effect_clean"], row["other_clean"]
        try:
            eaf = float(row["eaf_clean"]) if pd.notna(row["eaf_clean"]) else None
        except (TypeError, ValueError):
            eaf = None
        rsid = str(row["rsid_value"]).strip()

        manual_action = actions_dict.get(primary_key, actions_dict.get(fallback_key, "keep"))
        if manual_action == "drop":
            note = "Manual_Drop"
            audit_lines.append(_format_audit_row(
                row, effect=effect, other=other, eaf=eaf,
                palindromic=is_palindromic(effect, other), chosen=None,
                final_effect_allele="", status="Dropped", note=note,
                requested_build=args.build,
            ))
            drop_reasons[note] = drop_reasons.get(note, 0) + 1
            continue

        # Look for the variant in both the primary and fallback reference databases
        source_candidates: dict[str, list[dict[str, Any]]] = {"primary": [], "fallback": []}
        for db in dbs:
            source = db["source"]
            pos = primary_pos if source == "primary" else fallback_pos
            pos_key = primary_key if source == "primary" else fallback_key
            records = []
            
            if pos_key:
                records.extend(db["by_pos"].get(pos_key, []))
            if rsid and rsid.lower().startswith("rs"):
                records.extend(db["by_rsid"].get(rsid, []))
                
            records = _dedup_records(records)
            candidates = []
            for record in records:
                candidates.extend(_record_matches(
                    record, effect, other,
                    assume_forward=assume_forward,
                    effect_is_alt=effect_is_alt,
                ))
            source_candidates[source] = _dedup_candidates(candidates, source)

        # Pick the strongest match overall (Preferring exact matches over fallback builds)
        pool: list[dict[str, Any]] = []
        selected_strength = None
        for strength in (0, 1, 2):
            for source in ("primary", "fallback"):
                same_strength = [c for c in source_candidates[source] if c["strength"] == strength]
                if same_strength:
                    pool = same_strength
                    selected_strength = strength
                    break
            if pool:
                break

        final_pal = is_palindromic(effect, other)
        for candidate in pool:
            rec = candidate["rec"]
            if len(rec["alts"]) == 1 and is_palindromic(rec["ref"], rec["alts"][0]):
                final_pal = True

        chosen = None
        status = "Dropped"
        note = "NotFoundInReference" if not (p_pos or f_pos) else "UnmatchedAllele"
        
        if not pool:
            any_record = bool(
                (primary_key and p_pos.get(primary_key))
                or (fallback_key and f_pos.get(fallback_key))
                or (rsid and (p_rsid.get(rsid) or f_rsid.get(rsid)))
            )
            note = "OtherAllelePairMismatch" if any_record and other else ("UnmatchedAllele" if any_record else "NotFoundInReference")
        else:
            known_orientation = assume_forward or effect_is_alt
            # Address ambiguous Palindromic SNPs that need population frequency to resolve
            if final_pal and not known_orientation and (eaf is None or args.pal_low <= eaf <= args.pal_high):
                note = "AmbiguousPalindromic_NoEAF" if eaf is None else "AmbiguousPalindromic_5050"
                pool = []
            elif len({_candidate_key(c, c["source"]) for c in pool}) > 1:
                candidate_records = []
                for candidate in pool:
                    af_map = p_af if candidate["source"] == "primary" else f_af
                    candidate_records.append({
                        "key": _candidate_key(candidate, candidate["source"]),
                        "pop_afs": _candidate_population_afs(candidate, af_map, super_pops),
                    })
                label = "AmbiguousPalindromic" if final_pal else "AmbiguousRecord"
                key, note = resolve_candidates_global(
                    eaf, candidate_records, label=label,
                    af_low=args.pal_low, af_high=args.pal_high,
                    preferred_pop=args.eaf_pop,
                )
                if key is not None:
                    chosen = next(c for c in pool if _candidate_key(c, c["source"]) == key)
                    status = "Kept"
            elif len(pool) == 1:
                chosen, status = pool[0], "Kept"
                if selected_strength == 0:
                    note = "StrictPairMatch"
                elif effect_is_alt:
                    note = "EffectAlleleALTMatch" if selected_strength == 1 else "EffectAlleleALT_OtherMismatch"
                elif selected_strength == 1:
                    note = "EffectAlleleOnlyMatch"
                else:
                    note = "DirectMatch_OtherMismatch"
                if chosen["flipped"]:
                    note = "StrandFlipped"
            else:
                note = "AmbiguousRecord_Unresolved"

        # Apply final match and construct mapped rows
        final_effect_allele = ""
        if chosen is not None:
            record = chosen["rec"]
            allele_set = [record["ref"]] + record["alts"]
            final_effect_allele = allele_set[chosen["idx"]]
            kept += 1
            
            mapped = row.to_dict()
            mapped.update({
                "native_id": chosen["native_id"],
                "harmonized_eff": final_effect_allele,
                "weight_multiplier": 1.0,
                "source": chosen["source"],
                "matched_ref": record["ref"],
                "matched_alts": ",".join(record["alts"]),
                "matched_pos": record["pos"],
                "matched_chrom": record["chrom"],
                "match_strength": selected_strength,
                "flipped": bool(chosen["flipped"]),
                "matched_idx": chosen["idx"],
            })
            mapped_rows.append(mapped)
        else:
            drop_reasons[note] = drop_reasons.get(note, 0) + 1

        audit_lines.append(_format_audit_row(
            row, effect=effect, other=other, eaf=eaf,
            palindromic=final_pal, chosen=chosen,
            final_effect_allele=final_effect_allele,
            status=status, note=note, requested_build=args.build,
        ))

    mapped_df = pd.DataFrame(mapped_rows)
    audit_path = os.path.join(args.out_dir, "panel_harmonisation_audit.tsv")
    pd.DataFrame(audit_lines).to_csv(audit_path, sep="\t", index=False)
    
    log(f"Harmonization result: kept {kept}/{len(df)}; dropped {len(df) - kept}.")
    log(f"Drop reasons: {drop_reasons}")

    # Output formatted weights for downstream scoring
    weight_col = "effect_weight"
    if not mapped_df.empty:
        mapped_df[weight_col] = pd.to_numeric(mapped_df[weight_col], errors="raise") * mapped_df["weight_multiplier"]
        mapped_df[["native_id", "harmonized_eff", weight_col]].to_csv(
            os.path.join(args.out_dir, "formatted_panel.txt"), sep="\t", index=False, header=False
        )
    else:
        open(os.path.join(args.out_dir, "formatted_panel.txt"), "w", encoding="utf-8").close()


    # -------------------------------------------------------------------------
    # STEP 5: OUTPUT ALLELE FREQUENCY CACHE
    # -------------------------------------------------------------------------
    log("STEP 5: Writing unified allele-frequency cache...")
    if args.af_cache_out != "NONE" and not mapped_df.empty:
        af_parent = os.path.dirname(os.path.abspath(args.af_cache_out))
        os.makedirs(af_parent, exist_ok=True)
        af_records = []
        
        for _, row in mapped_df.iterrows():
            source_af = p_af if row["source"] == "primary" else f_af
            rec = {
                "chr": row["chr_clean"],
                "pos_hg19": row["active_pos_hg19"],
                "pos": row["active_pos_hg38"],
                "effect_clean": row["effect_clean"],
                "other_clean": row["other_clean"],
                "harmonized_effect_allele": row["harmonized_eff"],
                "REF": row["matched_ref"],
                "ALT": row["matched_alts"],
                "native_id": row["native_id"],
                "source": row["source"],
            }
            record_for_af = {
                "native_id": row["native_id"],
                "ref": row["matched_ref"],
                "alts": str(row["matched_alts"]).split(","),
            }
            af_entry = source_af.get(row["native_id"])
            if af_entry:
                allele_set = [record_for_af["ref"]] + record_for_af["alts"]
                if row["harmonized_eff"] in allele_set:
                    idx = allele_set.index(row["harmonized_eff"])
                    for pop in super_pops:
                        rec[f"{pop}_AF"] = freq_of(record_for_af, idx, pop, source_af)
                else:
                    for pop in super_pops:
                        rec[f"{pop}_AF"] = None
            else:
                for pop in super_pops:
                    rec[f"{pop}_AF"] = None
                    
            rec["query_id_hg19"] = f"{row['chr_clean']}:{row['active_pos_hg19']}_{row['other_clean']}>{row['effect_clean']}"
            rec["query_id_hg38"] = f"{row['chr_clean']}:{row['active_pos_hg38']}_{row['other_clean']}>{row['effect_clean']}"
            af_records.append(rec)
            
        pd.DataFrame(af_records).to_csv(args.af_cache_out, index=False)


    # -------------------------------------------------------------------------
    # STEP 6: SCORE REFERENCE POPULATION FOR Z-SCORE CALIBRATION
    # -------------------------------------------------------------------------
    log("STEP 6: Scoring the primary and fallback reference subsets...")
    if mapped_df.empty:
        raise RuntimeError("No variants were harmonized; refusing to write fake ancestry stats")

    score_dfs = []
    primary_additive = mapped_df[(mapped_df["source"] == "primary") & (mapped_df["is_additive"] == True)]
    fallback_additive = mapped_df[(mapped_df["source"] == "fallback") & (mapped_df["is_additive"] == True)]

    # Score primary build variants
    if not primary_additive.empty:
        if not os.path.exists(f"{primary_prefix}.pgen"):
            raise FileNotFoundError("Primary variants were harmonized but primary subset PGEN is missing")
            
        score_file = os.path.join(args.out_dir, "formatted_panel_primary.txt")
        primary_additive[["native_id", "harmonized_eff", weight_col]].to_csv(
            score_file, sep="\t", index=False, header=False
        )
        score_prefix = os.path.join(args.out_dir, "ref_scores_primary")
        
        subprocess.run([
            "plink2", "--pfile", primary_prefix, "--score", score_file,
            "1", "2", "3", "cols=+scoresums", "--out", score_prefix,
            "--threads", "8", "--memory", "16000",
        ], check=True, stdout=subprocess.DEVNULL)
        
        score_df = pd.read_csv(f"{score_prefix}.sscore", sep="\t")
        score_df.columns = [str(c).lstrip("#") for c in score_df.columns]
        
        if not {"IID", "SCORE1_SUM"}.issubset(score_df.columns):
            raise ValueError(f"Unexpected PLINK primary score columns: {list(score_df.columns)}")
            
        score_dfs.append(score_df[["IID", "SCORE1_SUM"]].rename(columns={"SCORE1_SUM": "SCORE_PRIMARY"}))
        log(f"Primary score: {len(score_df)} IIDs; {len(primary_additive)} panel rows.")

    # Score fallback build variants
    if not fallback_additive.empty:
        if not os.path.exists(f"{fallback_prefix}.pgen"):
            raise FileNotFoundError("Fallback variants were harmonized but fallback subset PGEN is missing")
            
        score_file = os.path.join(args.out_dir, "formatted_panel_fallback.txt")
        fallback_additive[["native_id", "harmonized_eff", weight_col]].to_csv(
            score_file, sep="\t", index=False, header=False
        )
        score_prefix = os.path.join(args.out_dir, "ref_scores_fallback")
        
        subprocess.run([
            "plink2", "--pfile", fallback_prefix, "--score", score_file,
            "1", "2", "3", "cols=+scoresums", "--out", score_prefix,
            "--threads", "8", "--memory", "16000",
        ], check=True, stdout=subprocess.DEVNULL)
        
        score_df = pd.read_csv(f"{score_prefix}.sscore", sep="\t")
        score_df.columns = [str(c).lstrip("#") for c in score_df.columns]
        
        if not {"IID", "SCORE1_SUM"}.issubset(score_df.columns):
            raise ValueError(f"Unexpected PLINK fallback score columns: {list(score_df.columns)}")
            
        score_dfs.append(score_df[["IID", "SCORE1_SUM"]].rename(columns={"SCORE1_SUM": "SCORE_FALLBACK"}))
        log(f"Fallback score: {len(score_df)} IIDs; {len(fallback_additive)} panel rows.")

    if not score_dfs:
        raise RuntimeError("No additive score files were generated; refusing to write ancestry_stats.json")

    # Combine scores across builds for individuals present in both datasets
    if len(score_dfs) == 1:
        combined_scores = score_dfs[0].copy()
        score_column = next(c for c in combined_scores.columns if c != "IID")
        combined_scores = combined_scores.rename(columns={score_column: "SCORE1_SUM"})
        score_cohort_note = "single reference subset"
    else:
        combined_scores = score_dfs[0].merge(
            score_dfs[1], on="IID", how="inner", validate="one_to_one"
        )
        if combined_scores.empty:
            raise RuntimeError("Primary and fallback score files have no overlapping IIDs")
            
        combined_scores["SCORE1_SUM"] = combined_scores["SCORE_PRIMARY"] + combined_scores["SCORE_FALLBACK"]
        combined_scores = combined_scores[["IID", "SCORE1_SUM"]]
        score_cohort_note = "intersection of primary and fallback IID sets (no imputation)"

    # Generate Mean and SD Statistics per SuperPopulation
    label_df = primary_psam[["IID", "SuperPop"]].copy()
    label_df["IID"] = label_df["IID"].astype(str)
    combined_scores["IID"] = combined_scores["IID"].astype(str)
    combined_scores = combined_scores.merge(label_df, on="IID", how="inner", validate="one_to_one")
    
    if combined_scores.empty:
        raise RuntimeError("No scored IIDs matched the primary PSAM population labels")
    log(f"Statistics cohort: {len(combined_scores)} IIDs; {score_cohort_note}.")

    stats_dict = {}
    for pop in super_pops:
        pop_scores = pd.to_numeric(
            combined_scores.loc[combined_scores["SuperPop"].astype(str) == str(pop), "SCORE1_SUM"],
            errors="coerce",
        ).dropna()
        if len(pop_scores) < 2:
            raise RuntimeError(f"Fewer than two scored individuals for population {pop}; cannot estimate SD")
            
        stats_dict[pop] = {
            "mean": float(pop_scores.mean()),
            "sd": float(pop_scores.std(ddof=1)),
        }

    stats_path = os.path.join(args.out_dir, "ancestry_stats.json")
    with open(stats_path, "w", encoding="utf-8") as handle:
        json.dump(stats_dict, handle, indent=4, allow_nan=False)

    log(f"Wrote {stats_path}.")
    log("Allele-match mode: " + ("effect-is-ALT" if effect_is_alt else "known-forward" if assume_forward else "full harmonization"))
    log(f"EAF resolution population: {args.eaf_pop}.")
    log(f"ALL FINISHED in {time.time() - t0:.1f}s.")


if __name__ == "__main__":
    main()