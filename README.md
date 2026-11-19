# PRS Pipeline & Streamlit App

A containerized Polygenic Risk Score (PRS) calculation pipeline powered by Snakemake with an interactive Streamlit web interface. It supports running via **Docker** (local workstations) or **Singularity / Apptainer** (shared servers and HPC clusters without root access).

## Repository Structure

```text
.
├── app/                 # Streamlit web application (main.py)
├── Pipeline/            # Snakemake workflow, scripts, and Conda environments
├── Dockerfile           # Container image definition
├── docker-compose.yml   # Docker Compose configuration
└── setup.sh             # Downloads required genomic reference files
```

---

## Storage & General Requirements

Plan for **at least 150 GB of free disk space** on a writable drive before setting up the repository:

| Component | Disk Size | Notes |
|---|---|---|
| **Reference datasets** (`references/`) | **~127 GB** | Downloaded by `setup.sh` (1000G GRCh38 PLINK ~85 GB, 1000G GRCh37 PLINK ~22 GB, Beagle GRCh38/GRCh37 `bref3` ~14 GB, Somalier ~530 MB, plus reference FASTAs) |
| **Container image** (`app.sif` / Docker image) | **~740 MB** | Final `.sif` image size (building the image temporarily requires **~15 GB** of scratch space, which can be deleted afterward) |
| **Snakemake environments & caches** (`Pipeline/`, `conda_pkgs/`) | **~7 GB** | Created automatically on the first pipeline run and reused afterward |

* **Software:** Either **Docker + Docker Compose** *or* **Singularity / Apptainer**
* **CPU / RAM:** Multi-core CPU (8 cores recommended; used by Snakemake/Beagle) and at least 16 GB RAM for Beagle imputation

---

## Step 1: Download Reference Files (Required First)

Before launching the app with either Docker or Singularity, download and configure the reference genomes and annotation datasets (~127 GB) from the root of this repository:

```bash
chmod +x setup.sh
./setup.sh
```

Wait until `setup.sh` completes before building or running the container.

---

## Option A: Run with Docker

From the repository root directory, build and start the container:

```bash
docker compose up --build
```

> **Saving pipeline outputs directly into the repo folder (optional):** By default, pipeline runs write to `/tmp-streamlit` inside the container. To persist all outputs directly inside `./tmp-streamlit` on your host machine, add `- ./tmp-streamlit:/tmp-streamlit` under `volumes:` in `docker-compose.yml` before running `docker compose up --build`.

Once started, open **http://localhost:8501** in your browser. Press `Ctrl + C` in the terminal to stop the app.

---

## Option B: Run with Singularity / Apptainer (HPC & Shared Servers)

Use this method if Docker is unavailable on your server. Singularity builds a standalone `.sif` image (~740 MB) and runs without requiring root privileges at runtime.

### 1. Convert the `Dockerfile` to a Singularity Definition File
From the repository root directory, install `spython` and generate `Singularity.def`:

```bash
pip install spython
spython recipe Dockerfile > Singularity.def
```

### 2. Build `app.sif`

Choose the build method that matches your server permissions:

#### Method 1: With `sudo` or `--fakeroot` (Standard)
If you have admin rights or `--fakeroot` is configured on your server/compute node:

```bash
# With sudo (personal Linux machine)
sudo singularity build app.sif Singularity.def

# OR with fakeroot (HPC cluster with subuid/subgid configured)
singularity build --fakeroot app.sif Singularity.def
```

#### Method 2: Without `sudo` or `--fakeroot` (Unprivileged Build with Static `proot`)
If neither `sudo` nor `--fakeroot` is available, use the statically linked `proot` binary and redirect Singularity's cache and temporary build directory (~15 GB needed during the build) to drives with sufficient space:

```bash
# 1. Download the static proot binary (avoids missing shared-library errors inside the container)
mkdir -p ~/bin
curl -L [https://proot.gitlab.io/proot/bin/proot](https://proot.gitlab.io/proot/bin/proot) -o ~/bin/proot
chmod +x ~/bin/proot
export PATH="$HOME/bin:$PATH"

# 2. Point Singularity's cache and temporary build folders to drives with sufficient space
#    - Cache (~1-2 GB): use $HOME (must allow chmod 700)
#    - Temp dir (~15 GB): use the current storage mount
mkdir -p "$HOME/singularity_build/cache" "$(pwd)/singularity_tmp"
chmod 700 "$HOME/singularity_build/cache"

export SINGULARITY_CACHEDIR="$HOME/singularity_build/cache"
export SINGULARITY_TMPDIR="$(pwd)/singularity_tmp"
export APPTAINER_CACHEDIR="$SINGULARITY_CACHEDIR"
export APPTAINER_TMPDIR="$SINGULARITY_TMPDIR"
export TMPDIR="$SINGULARITY_TMPDIR"

# 3. Build the image
singularity build app.sif Singularity.def
```

Once `Build complete: app.sif` appears, clean up the temporary build files to reclaim disk space:

```bash
rm -rf "$(pwd)/singularity_tmp" "$HOME/singularity_build"
```

#### Method 3: Build on a Local Machine with Docker and Transfer
If building directly on the server is restricted, build the image on a computer with Docker, transfer the archive, and convert it on the server (converting a `.tar` archive requires no root or `proot`):

```bash
# On your local machine (inside the repo folder)
docker build -t prs-app .
docker save prs-app -o prs-app.tar
scp prs-app.tar <user>@<server>:/path/to/repo/

# On the server (inside the repo folder)
singularity build app.sif docker-archive://prs-app.tar
rm -f prs-app.tar
```

### 3. Prepare Persistent Folders and Launch the App

Create local folders inside the repository so Snakemake environments, Conda package caches, temporary runtime files, and **all pipeline results** are stored permanently on disk instead of inside temporary container memory:

```bash
mkdir -p snakemake_envs conda_pkgs tmp_runtime tmp-streamlit
```

Launch the Streamlit app from the repository root directory:

```bash
singularity exec --cleanenv --writable-tmpfs \
  --env "PATH=/opt/conda/envs/prs-app-env/bin:/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin" \
  --env "CONDA_PKGS_DIRS=/prs-app/conda_pkgs" \
  --env "TMPDIR=/prs-app/tmp_runtime" \
  --env "CONDA_REPODATA_USE_SHARDS=false" \
  --bind "$(pwd):/prs-app" \
  --bind "$(pwd)/snakemake_envs:/prs-app/.snakemake" \
  --bind "$(pwd)/tmp-streamlit:/tmp-streamlit" \
  --pwd /prs-app/app \
  app.sif \
  /opt/conda/envs/prs-app-env/bin/python -m streamlit run main.py --server.port=8501 --server.address=0.0.0.0
```

**Why these flags matter:**
* `--bind "$(pwd)/tmp-streamlit:/tmp-streamlit"` saves all pipeline results directly to `./tmp-streamlit/` inside your repository folder instead of losing them in temporary container memory when the container stops.
* `--pwd /prs-app/app` ensures relative paths to `../Pipeline/Snakefile` resolve correctly.
* `--env "CONDA_PKGS_DIRS=..."` and `--env "TMPDIR=..."` keep package downloads off limited `/tmp` or tmpfs memory.
* `--env "CONDA_REPODATA_USE_SHARDS=false"` prevents SQLite `database is locked` errors on network/HPC filesystems.

Leave this terminal open while using the app. Press `Ctrl + C` to stop it.

---

## Accessing the Web Interface

* **Local machine:** Open `http://localhost:8501` in your web browser.
* **Remote server / HPC cluster:** Open a **new terminal on your local computer** to forward port `8501`:
  ```bash
  ssh -L 8501:localhost:8501 <your_username>@<server_address>
  ```
  Then open `http://localhost:8501` in your local browser.

---

## Outputs & Caches

When launched with the Singularity command above (or with `./tmp-streamlit:/tmp-streamlit` mounted in Docker), all results and caches are saved inside the repository folder:

| Path | Contents |
|---|---|
| `tmp-streamlit/prs_imputed_full/` | Final PRS score tables (`*_prs_results.csv`) and merged VCFs |
| `tmp-streamlit/beagle_full/` | Imputed VCFs (`*_imputed.vcf.gz`) and indexes |
| `tmp-streamlit/ancestry_raw_full/` | Somalier ancestry predictions |
| `tmp-streamlit/leftalign_raw_full/` | Normalized and left-aligned VCFs |
| `tmp-streamlit/configs/` | Per-run variant action configs |
| `Pipeline/inputs/reference_cache/` | Cached reference calibration distributions per PGS panel |
| `Pipeline/.snakemake/conda/` | Snakemake rule Conda environments (created once on first run) |
| `conda_pkgs/` | Cached Conda package archives (Singularity) |

> **First-run note:** The first time you run a pipeline for a new PGS panel, Snakemake will create its rule-specific Conda environments and build the reference calibration cache (which can take ~10 minutes). Subsequent runs reuse both caches automatically and finish much faster.