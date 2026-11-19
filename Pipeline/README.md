# PRS Calculation Snakemake Pipeline

This directory contains the Snakemake workflow (`Snakefile`), helper scripts (`scripts/`), reference configuration (`reference_config.yaml`), and per-rule Conda environments (`envs/`) used to process sequencing data (VCF, gVCF, or BAM) into ancestry-calibrated Polygenic Risk Scores (PRS).

## Workflow Overview

```text
                       ┌───────────────────────────────┐
                       │ PGS Panel + 1000G References  │
                       └───────────────┬───────────────┘
                                       ▼
                      Calibrate_Reference_Distributions
                             (panel cache + stats)
                                       │
Input (BAM) ─► CallVariantsGATK_from_BAM ─► compress_index_vcf ─┐
                                                                ▼
Input (VCF / gVCF) ────────────────────────────────────────► LEFTALIGN
                                                        ┌───────┴───────┐
                                                        ▼               │
                                                Somalier_Extract        │
                                                        ▼               │
                                                Somalier_Ancestry       │
                                                        ▼               │
                                               Build_Ancestry_File      │
                                                        │               │
                                                        ▼               ▼
                                             Prepare_PRS_and_Targets ◄──┘
                                                │               │
                                          (_prep.vcf)    (_targets.txt)
                                                │               ▼
                                                │   BeagleImputation_Targeted
                                                │       (_imputed.vcf.gz)
                                                │               │
                                                └──────► Merge_Beagle ◄─┘
                                                        (_merged.vcf)
                                                                ▼
                                                        Score_PRS_Imputed
                                                       (_prs_results.csv)
```

---

## Pipeline Rules (Step-by-Step)

### 1. Reference Calibration (`Calibrate_Reference_Distributions`)
* **Script / Env:** `scripts/build_reference_cache.py` (`envs/plink2.yaml`)
* **What it does:** Harmonizes the selected PGS scoring file against primary and fallback 1000 Genomes PLINK reference panels (matching `GRCh38`/`hg38` or `GRCh37`/`hg19`). Generates allele-frequency caches, a harmonisation audit table (`panel_harmonisation_audit.tsv`), and reference population PRS distributions (`ancestry_stats.json`) used for Z-score calibration. Writes `calibration.done` once validated.

### 2. Variant Calling from BAM (`CallVariantsGATK_from_BAM` & `compress_index_vcf`)
* **Tools / Env:** `gatk HaplotypeCaller` (`envs/gatk.yaml`), `bgzip` + `tabix` (`envs/bcftools.yaml`)
* **What it does:** Runs **only** when `input_type=bam` (skipped for VCF/gVCF inputs). Calls variants across the full genome for WGS or restricted to `panel_coords_primary.bed` when `seq_type` contains `"Targeted"`, then compresses and indexes the output VCF.

### 3. Normalization & Left-Alignment (`LEFTALIGN`)
* **Tools / Env:** `bcftools` (`>=1.16`), `tabix` (`envs/leftalign.yaml`)
* **What it does:**
  1. Fixes VCF header type definitions (`QR` and `QA` `Integer` → `Float`) to prevent parser errors on FreeBayes VCFs.
  2. Runs `bcftools norm -a -c w -f <ref> -m -any` to atomize complex variants, left-align indels against the reference FASTA, and split multiallelic sites into biallelic rows.
  3. Renames the VCF sample column to `base_sample_name` and creates a `.tbi` index.

### 4. Genetic Ancestry Inference (`Somalier_Extract`, `Somalier_Ancestry`, `Build_Ancestry_File`)
* **Tools / Env:** `somalier` (`envs/somalier.yaml`) + local Python rule
* **What it does:** Extracts ancestry-informative sites from the normalized VCF and projects the sample onto 1000 Genomes reference sketches to predict genetic ancestry and PCs.
  * If `known_ancestry` is provided, it overrides the predicted population label while keeping Somalier PCs if available.
  * If Somalier fails (e.g., too few overlapping sites on a small targeted panel) and no override is set, it gracefully falls back to `EUR` and records the fallback in `ancestry_notes`.

### 5. Variant QC & Imputation Target Selection (`Prepare_PRS_and_Targets`)
* **Script / Env:** `scripts/prep_vcf.py` (`envs/python.yaml`)
* **What it does:** Filters variants using configured depth (`min_dp`), genotype quality (`min_gq`), allele balance (`min_af`), gnomAD low-complexity regions (LCRs), palindromic frequency bounds (`pal_low`/`pal_high`), strand/effect flags, and any manual `variant_actions`. Outputs:
  * `{sample}_prep.vcf`: QC-filtered observed genotypes.
  * `{sample}_targets.txt`: Chromosomes containing missing PGS variants that require Beagle imputation (empty if `skip_beagle=True` or no targets are missing).

### 6. Targeted Imputation (`BeagleImputation_Targeted`)
* **Tools / Env:** `bcftools`, `beagle`, `tabix` (`envs/bcftools.yaml`, 8 threads, `-Xmx16g`)
* **What it does:** Filters the normalized VCF to biallelic SNPs/indels and runs Beagle chromosome-by-chromosome only for chromosomes listed in `{sample}_targets.txt` using the corresponding `.bref3` reference panels.
  * Concatenates imputed chromosomes into `{sample}_imputed.vcf.gz` (or writes a valid empty VCF if no chromosomes need imputation).
  * Fails immediately on filesystem write errors (`No space left on device`); for algorithm/reference failures on individual chromosomes, logs a warning and allows downstream AF fallback.

### 7. Merge & Final PRS Scoring (`Merge_Beagle` & `Score_PRS_Imputed`)
* **Scripts / Env:** `scripts/merge_post_beagle.py`, `scripts/compute_prs_zscore.py` (`envs/python.yaml`)
* **What it does:**
  1. Merges observed genotypes from `_prep.vcf` with imputed genotypes from `_imputed.vcf.gz` that meet the dosage $R^2$ threshold (`beagle_dr2`).
  2. Calculates the sample's PRS and ancestry-calibrated Z-score against `ancestry_stats.json` (bypassed when `af_mode=CUSTOM`), writing the final table to `{sample}_prs_results.csv`.

---

## Key Configuration Parameters

Parameters are passed via Snakemake `--config key=value` (typically generated automatically by `app/main.py`):

| Parameter | Default | Description |
|---|---|---|
| `input_file` | `""` | Path to input `.vcf`, `.vcf.gz`, `.g.vcf(.gz)`, or `.bam` |
| `input_type` | `vcf` | Input format: `vcf`, `gvcf`, or `bam` (`.g.vcf` files are also auto-detected) |
| `seq_type` | `WGS` | Sequencing type (e.g., `Whole Genome Sequencing (WGS)` or `Targeted Sequencing (Panel)`) |
| `genome_build` | `GRCh38` | Genome assembly: `GRCh38` (`hg38`) or `GRCh37` (`hg19`) |
| `pgs_panel` / `pgs_file` | — | PGS Catalog panel ID and path to harmonized scoring file |
| `af_mode` | `1000G` | Allele frequency source: `1000G`, `PGS_FILE`, or `CUSTOM` |
| `known_ancestry` | `RUN_SOMALIER` | Leave as `RUN_SOMALIER` for automatic prediction, or supply a manual population override |
| `assume_forward` / `effect_is_alt` | `False` | Strand/allele harmonization toggles |
| `skip_beagle` | `False` | Skip Beagle imputation and rely solely on observed calls + AF fallback |
| `min_dp` / `min_gq` / `min_af` | `15` / `50` / `0.20` | Variant QC thresholds for depth, genotype quality, and allele fraction |
| `beagle_dr2` | `0.90` | Minimum Beagle imputation $DR^2$ required to retain an imputed site |
| `pal_low` / `pal_high` | `0.40` / `0.60` | Ambiguous allele-frequency window for palindromic SNPs |
| `variant_actions` | `""` | Optional TSV of manual per-variant overrides (include/exclude/flip) |
| `reference_config` | `/prs-app/Pipeline/reference_config.yaml` | YAML mapping reference FASTAs, 1000G PLINK prefixes, LCR BEDs, Somalier files, and Beagle panels |
| `tmp_dir` | `/prs-app/tmp` | Root directory for intermediate files, caches, and final outputs |

---

## Caching & Fingerprinting

The workflow uses deterministic cache keys so expensive steps are not repeated unnecessarily, while automatically invalidating caches if inputs or parameters change:

* **Reference Calibration Cache (`Pipeline/inputs/reference_cache/<pgs_panel>/`):**
  * Shared base coordinate subsets (`<panel>_<build>_base`) are reused across parameter variants.
  * Mode-specific caches encode parameter suffixes (`_fwd`, `_alt`, `_pgsaf`, `_act_<sha256>`) and a 12-character SHA-256 fingerprint of the PGS file and 1000G reference file sizes/timestamps.
  * The same suffix string is appended to output sample filenames (`SAMPLE`) so runs with different harmonization settings never overwrite each other.
* **Shared WGS Caches (enabled when `seq_type="Whole Genome Sequencing (WGS)"`):**
  * Keyed by `wgs_cache_key` and validated via `.manifest` files storing input SHA-256 checksums and reference signatures.
  * Reuses **left-aligned VCFs** (`leftalign_wgs_cache/`), **Somalier ancestry results** (`somalier_wgs_cache/`), and **per-chromosome Beagle outputs** (`beagle_wgs_cache/`) across different PGS panels run on the same whole genome.

---

## Final Outputs (`rule all`)

For each run, Snakemake produces two primary targets under `{tmp_dir}`:

* `{tmp_dir}/prs_imputed_full/{SAMPLE}_prs_results.csv` — Final PRS results and ancestry-calibrated Z-scores
* `{tmp_dir}/beagle_full/{SAMPLE}_imputed.vcf.gz` (plus `.tbi`) — Consolidated Beagle-imputed VCF

Intermediate files are organized under `{tmp_dir}` in `variant_calling_full/`, `leftalign_raw_full/`, `ancestry_raw_full/`, `prs_prep_full/`, and `prs_imputed_full/`.