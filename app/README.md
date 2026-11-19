# PRS Streamlit Web Application (`app/`)

This directory contains the Streamlit web interface (`main.py`) that configures, launches, monitors, and visualizes results from the Snakemake workflow in `../Pipeline/Snakefile`.

## Files

| File | Role |
|---|---|
| `main.py` | Main Streamlit application (UI configuration, Snakemake runner, live log parser, and results dashboard) |
| `pgs_fetcher.py` | Fetches PGS Catalog metadata and harmonized scoring files by PGS ID |
| `environment.yml` | Conda environment specification for running the web app |

---

## Application Modes

The sidebar top toggle switches between two modes:

### 1. Run New Pipeline
Guides the user through a 7-step configuration workflow before launching Snakemake:

1. **Panel Selection & Editor**
   * Scans `../Pipeline/inputs/*.txt.gz` and displays cache status for each panel (`Ready` when `reference_cache/<panel>/<mode>/calibration.done` exists, `Partially Ready` when an AF cache CSV exists, or `Cache Missing`).
   * Supports selecting a local panel, fetching a new score from the **PGS Catalog** by ID (with Discovery, Training, and Evaluation ancestry metadata), or uploading a **custom panel** (`.csv`, `.tsv`, or `.txt` containing `chr`, `pos`, `ref`, `alt`, and `effect_weight`, `OR`, or `HR`, saved as `CUSTOM_<name>.txt.gz`).
   * **Advanced Panel & AF Editor:** Provides harmonization toggles (`assume_forward`, `effect_is_alt`, `skip_beagle`) and an interactive table to edit values or assign per-variant actions (`Keep`, `Drop`, `Force Impute`). Edited panels can be saved as reusable `CUSTOM_*` files.
2. **Global Configuration:** Selects genome assembly (`GRCh38 (hg38)` or `GRCh37 (hg19)`).
3. **Population Allele Frequencies:** Selects `PGS_FILE` (shown automatically when `allelefrequency_effect` is present in the panel), `1000G` (uses or builds the local 1000 Genomes cache), `CUSTOM` (upload `.csv`/`.tsv`), or `NONE` (Raw PRS only).
4. **Patient Ancestry:** Runs automatic ancestry inference via Somalier (`RUN_SOMALIER`) or applies a manual ancestry label (`EUR`, `AFR`, `AMR`, `EAS`, `SAS`, or `Other / Unknown`).
5. **Patient Data Input:** Selects sequencing type (`Whole Genome Sequencing (WGS)` or `Targeted Sequencing (Panel)`), optional custom sample name, and input source—either browser upload (saved to `../Pipeline/data/inputs/`) or an absolute server file path for large `.vcf`, `.vcf.gz`, `.g.vcf`, or `.bam` files.
6. **QC Thresholds & Run Config:** Configures `min_dp` (default `15`), `min_gq` (default `50`), `min_af` (default `0.20`), `beagle_dr2` (default `0.90`), CPU cores (default `8`), and palindromic EAF discard bounds (`0.40–0.60`).
7. **Execute:** Writes any non-`Keep` variant actions to `../../tmp-streamlit/configs/<sample>_actions.tsv` and launches `snakemake` inside `../Pipeline` via `subprocess.Popen`.

During execution, the main panel parses Snakemake's stdout in real time to display:
* **Pipeline Steps:** Scheduled DAG rules and job counts.
* **Currently Executing:** Active rules as jobs start and finish.
* **Detailed Terminal Log:** Expandable live log (rolling buffer of the last 1,500 lines).

### 2. Visualize Existing Results
Renders the full results dashboard from previously generated outputs without running Snakemake. Accepts drag-and-drop uploads for:
* **Required:** PRS results file (`*_prs_results.csv`)
* **Optional:** Somalier ancestry probabilities (`*.somalier-ancestry.tsv`) and interactive PCA report (`*.somalier-ancestry.html`)

---

## Results Dashboard

When a pipeline finishes (or when existing results are loaded), `render_dashboard()` displays up to four sections:

1. **PRS Summary Table & Export:** Displays formatted PRS metrics (rounded to 3 decimal places) with a download button for the complete unmodified `*_prs_results.csv` (including hidden internal distribution metrics).
2. **Contextualized Population Mapping:**
   * Metric cards for **Deterministic Min**, **95% UI Lower**, **Categorical Z-Score (`Z1`)**, **95% UI Upper**, **Deterministic Max**, and **PC-Adjusted Z-Score (`Z2`)** when available for WGS admixture adjustment.
   * Interactive Plotly standard normal distribution showing the population mean, patient `Z1`/`Z2` lines, shaded 95% uncertainty interval, and deterministic bounds.
3. **Somalier Ancestry Prediction:** Population confidence bar chart and probability table parsed from `*_prob` columns (shown when Somalier output is available).
4. **Interactive Ancestry PCA Map:** Embeds the Somalier HTML PCA visualization with the legend enabled and the patient sample highlighted as a star marker.

---

## Path & Working Directory Requirements

`main.py` resolves pipeline and output paths relative to the `app/` directory:

| Relative Path Used in `main.py` | Resolved Location | Purpose |
|---|---|---|
| `../Pipeline` | `<repo>/Pipeline` | Working directory for Snakemake execution |
| `../Pipeline/inputs` | `<repo>/Pipeline/inputs` | PGS panels, custom panels, and `reference_cache/` |
| `../Pipeline/data/inputs` | `<repo>/Pipeline/data/inputs` | Destination for browser-uploaded patient files |
| `../../tmp-streamlit` | `/tmp-streamlit` (in container) | Pipeline outputs, intermediate files, and `configs/` |

Because of these relative paths, always launch Streamlit with `app/` as the working directory (e.g., `--pwd /prs-app/app` in Singularity or `WORKDIR /prs-app/app` in Docker):

```bash
python -m streamlit run main.py --server.port=8501 --server.address=0.0.0.0
```