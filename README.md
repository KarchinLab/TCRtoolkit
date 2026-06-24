# TCR-Toolkit-SCRATCH

A unified Nextflow DSL2 pipeline for T-cell receptor (TCR) analysis that handles **bulk**, **single-cell**, and **combined** TCR data under a single entry point.

Internally, the pipeline integrates two complementary analysis engines:

- **TCRtoolkit** — bulk and pseudo-bulk repertoire analysis (clonotype statistics, generation probabilities, TCR sharing, convergence, antigen specificity)
- **SCRATCH-TCR** — single-cell GEX + TCR co-analysis (VDJ QC, T-cell integration, clonotype clustering, repertoire profiling)

Users interact with one pipeline and one set of outputs. The choice of which engines run is determined automatically from the inputs provided, or set explicitly with `--mode`.

---

## Table of Contents

- [Overview](#overview)
- [Requirements](#requirements)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Pipeline Modes](#pipeline-modes)
  - [Mode 1: Bulk](#mode-1-bulk)
  - [Mode 2: Single-cell](#mode-2-single-cell)
    - [Full SC sub-mode](#full-sc-sub-mode)
    - [VDJ-only sub-mode](#vdj-only-sub-mode)
- [Parameters](#parameters)
- [Outputs](#outputs)
- [Architecture](#architecture)
- [Bridge Modules](#bridge-modules)
- [Configuration](#configuration)
- [Contact](#contact)

---

## Overview

| Mode | Sub-mode | Input | What runs |
|---|---|---|---|
| `bulk` | — | AIRR / Adaptive / CellRanger TSV files | Full TCRtoolkit repertoire analysis |
| `singlecell` | **Full SC** | Seurat RDS + Cell Ranger VDJ outs + sample sheet | VDJ QC → T-cell integration → TCRtoolkit bulk analysis → SC clustering (CONGA, GLIPH2, TCRdist3, GIANA) → Consensus → Master Summary |
| `singlecell` | **VDJ-only** | Cell Ranger VDJ outs + sample sheet *(no GEX object)* | VDJ QC → TCRtoolkit bulk analysis only; SC clustering modules are skipped |

The mode is auto-detected from the inputs you provide. Within `singlecell` mode, the sub-mode is also auto-detected: if `--input_annotated_object` is absent, the pipeline switches to VDJ-only automatically. You can also set the mode explicitly with `--mode bulk | singlecell | combined`.

---

## Requirements

- [Nextflow](https://www.nextflow.io/) ≥ 22.10
- Java 11 or later
- [Docker](https://docs.docker.com/engine/install/) **or** [Singularity](https://sylabs.io/singularity/)
- Git

Two containers are used automatically — no manual installation of tools required:

| Container | Used by |
|---|---|
| `ghcr.io/karchinlab/tcrtoolkit:main` | Bulk analysis processes |
| `syedsazaidi/scratch-tcr:latest` | Single-cell analysis processes |

---



## Installation

```bash
git clone https://github.com/KarchinLab/TCRtoolkit.git
cd TCR-Toolkit
```
---

## Quick Start

### Mode 1 — Bulk only

```bash
nextflow run main.nf \
  --samplesheet samplesheet.csv \
  --input_format airr \
  --workflow_level "sample,compare" \
  --outdir results/bulk \
  --project_name my_bulk_run
```

### Mode 2a — Single-cell (Full SC, with GEX object)

```bash
nextflow run main.nf \
  --input_annotated_object annotated_seurat.RDS \
  --input_vdj_contigs "data/VDJ/*/outs" \
  --sample_sheet sc_samplesheet.csv \
  --outdir results/singlecell \
  --project_name my_sc_run
```

### Mode 2b — Single-cell (VDJ-only, no GEX object)

```bash
nextflow run main.nf \
  --input_vdj_contigs "data/VDJ/*/outs" \
  --sample_sheet sc_samplesheet.csv \
  --outdir results/singlecell \
  --project_name my_sc_run
```

Omitting `--input_annotated_object` triggers VDJ-only mode automatically.


### Profile

Add `-profile singularity or docker` to any of the commands above:

```bash
nextflow run main.nf -profile singularity/docker \
  --input_annotated_object annotated_seurat.RDS \
  --input_vdj_contigs "data/VDJ/*/outs" \
  --sample_sheet sc_samplesheet.csv \
  --outdir results
```


---

## Pipeline Modes

### Mode 1: Bulk

For bulk or pseudo-bulk TCR sequencing data. Runs the full TCRtoolkit repertoire analysis pipeline.

```
Input (AIRR / Adaptive / CellRanger TSVs)
  │
  ├── INPUT_CHECK          — validate samplesheet
  ├── CONVERT (optional)   — Adaptive or CellRanger → AIRR format
  ├── ANNOTATE             — CDR3 deduplication + OLGA generation probabilities
  │
  ├── SAMPLE               — per-sample statistics
  │     ├── V/D/J gene usage
  │     ├── TCRdist3 distance matrices + histograms
  │     ├── OLGA pGen histograms
  │     ├── Convergence analysis
  │     ├── TCR phenotyping
  │     └── VDJdb antigen-specificity matching
  │
  ├── PATIENT (optional)   — patient-level analysis
  │     ├── GIANA clonotype clustering
  │     └── GLIPH2 motif clustering (optional)
  │
  └── COMPARE (optional)   — cross-sample TCR sharing
        ├── Sharing histogram
        └── Sharing vs log10(pGen) scatter plot
```

**Samplesheet format (bulk):**

```csv
sample,file,timepoint,origin,patient,subject_id
sample_A,/path/to/sample_A_airr.tsv,pre,blood,patient_1,P1
sample_B,/path/to/sample_B_airr.tsv,post,blood,patient_1,P1
```

---

### Mode 2: Single-cell

The single-cell mode has two sub-modes selected automatically based on whether a GEX-annotated Seurat object is provided.

---

#### Full SC sub-mode

Requires Cell Ranger VDJ output **and** a pre-processed Seurat RDS with GEX annotations. Runs the complete SCRATCH-TCR analysis alongside TCRtoolkit.

```
Input: Seurat RDS + Cell Ranger VDJ outs + sample sheet
  │
  ├── VDJ QC               — filter contigs (productive, high-confidence, UMI/CDR3 thresholds)
  ├── T-cell Integration   — merge QC'd contigs into Seurat; assign clonotypes; subset T cells
  │
  ├── [Bridge 1: SC → Bulk]
  │     Converts per-cell export table → per-sample AIRR TSVs
  │     ↓
  │   TCRtoolkit bulk analysis on SC-derived data
  │     ├── ANNOTATE → OLGA pGen
  │     ├── SAMPLE   → diversity, convergence, VDJdb
  │     ├── PATIENT  → GIANA / GLIPH2 clustering
  │     └── COMPARE  → TCR sharing across samples
  │
  ├── SC downstream modules (parallel with TCRtoolkit):
  │     ├── TCRi       — TCR-based immune phenotyping scores
  │     ├── CoNGA      — GEX + TCR graph co-analysis
  │     ├── GLIPH2     — motif-based antigen specificity clustering
  │     ├── TCRdist3   — distance-based clonotype clustering
  │     ├── GIANA      — sequence-similarity clonotype clustering
  │     └── Repertoire — clonal diversity, sharing, flux
  │
  ├── Consensus Clustering — majority-vote across GLIPH2, TCRdist3, GIANA
  └── Master Summary       — aggregated HTML report
```

**Selectively disable SC modules:**

```bash
nextflow run main.nf \
  --input_annotated_object seurat.RDS \
  --input_vdj_contigs "VDJ/*/outs" \
  --sample_sheet samplesheet.csv \
  --run_tcri false \
  --run_conga false \
  --run_gliph2 true \
  --run_tcrdist3 true \
  --run_giana true \
  --run_consensus true \
  --run_repertoire true \
  --run_master_summary true
```

---

#### VDJ-only sub-mode

For when you have Cell Ranger VDJ output but **no GEX Seurat object**. T-cell integration and all SC clustering modules are skipped. TCRtoolkit bulk analysis is still performed on the VDJ contigs.

Triggered automatically when `--input_annotated_object` is not provided.

```
Input: Cell Ranger VDJ outs + sample sheet  (no Seurat RDS)
  │
  ├── VDJ QC               — filter contigs (productive, high-confidence, UMI/CDR3 thresholds)
  │
  ├── [Bridge 3: VDJ → Bulk]
  │     Converts contigs_after_qc.tsv (TRB rows) → per-sample AIRR TSVs
  │     ↓
  │   TCRtoolkit bulk analysis on VDJ-derived data
  │     ├── ANNOTATE → OLGA pGen
  │     ├── SAMPLE   → diversity, convergence, VDJdb
  │     ├── PATIENT  → GIANA / GLIPH2 clustering
  │     └── COMPARE  → TCR sharing across samples
  │
  └── (SC modules not run — no GEX object available)
```

> **Note:** In VDJ-only mode, clonotype counts are derived by grouping cells with the same CDR3β + V + J gene. This gives repertoire-level TCR analysis equivalent to bulk sequencing, but without any GEX-informed cell type annotations.

---

**Sample sheet format (both single-cell sub-modes):**

```csv
sample,path
HRS371754,/data/SCRATCH_ALIGN-CELLRANGER_VDJ/HRS371754/outs
HRS371755,/data/SCRATCH_ALIGN-CELLRANGER_VDJ/HRS371755/outs
```



## Parameters

### Universal

| Parameter | Default | Description |
|---|---|---|
| `--mode` | auto-detected | `bulk` \| `singlecell` \| `combined` |
| `--project_name` | timestamped | Name used in report titles and output directories |
| `--outdir` | `results` | Output directory |
| `--max_cpus` | `20` | Max CPUs per process |
| `--max_memory` | `200.GB` | Max memory per process |
| `--max_time` | `240.h` | Max walltime per process |

### Shared metadata columns

These column names are used in both the bulk samplesheet and the Seurat object metadata:

| Parameter | Default | Description |
|---|---|---|
| `--sample_col` | `orig.ident` | Sample identifier column |
| `--patient_col` | `patient_id` | Patient identifier column |
| `--condition_col` | `condition` | Condition/group column |
| `--timepoint_col` | `timepoint` | Timepoint column |
| `--batch_col` | `batch` | Batch column |

### Bulk (TCRtoolkit) parameters

| Parameter | Default | Description |
|---|---|---|
| `--samplesheet` | — | Path to bulk samplesheet CSV |
| `--input_format` | `airr` | `airr` \| `adaptive` \| `cellranger` |
| `--workflow_level` | `sample,compare` | Comma-separated: `sample`, `patient`, `compare`, `convert` |
| `--sobject_gex` | — | Seurat RDS for pseudo-bulk-by-phenotype (CellRanger mode only) |
| `--olga_chunk_length` | `100000` | OLGA parallelization chunk size |
| `--matrix_sparsity` | `sparse` | TCRdist3 matrix sparsity |
| `--use_gliph2` | `false` | Run GLIPH2 at patient level |
| `--threshold` | `7.0` | GIANA similarity threshold |

### Single-cell (SCRATCH-TCR) parameters

| Parameter | Default | Description |
|---|---|---|
| `--input_annotated_object` | — | Path to annotated Seurat `.RDS` *(optional — omit for VDJ-only mode)* |
| `--input_vdj_contigs` | — | Glob path to Cell Ranger VDJ `outs/` directories |
| `--sample_sheet` | — | Path to SC sample sheet CSV |
| `--run_tcri` | `true` | Run TCRi module |
| `--run_conga` | `true` | Run CoNGA module |
| `--run_gliph2` | `true` | Run GLIPH2 SC module |
| `--run_tcrdist3` | `true` | Run TCRdist3 SC module |
| `--run_giana` | `true` | Run GIANA SC module |
| `--run_consensus` | `true` | Run Consensus clustering |
| `--run_repertoire` | `true` | Run Repertoire module |
| `--run_master_summary` | `true` | Run Master Summary |
| `--vdj_require_productive` | `true` | Keep productive contigs only |
| `--vdj_keep_paired_only` | `false` | Keep only paired alpha-beta contigs |
| `--clone_call_preference` | `aa` | Clonotype definition: `aa` or `nt` |
| `--gliph_reference_bundle` | built-in | Path to GLIPH2 reference `.RData` |
| `--tcrdist_radius` | `24` | TCRdist3 neighbor radius |
| `--consensus_min_methods` | `2` | Minimum methods agreeing for consensus label |

Full parameter documentation is in `nextflow.config`.

---

## Outputs

All results are written to `--outdir`. The directory structure depends on the mode:

### Bulk mode outputs

```
results/
├── sample_stats/         — per-sample V/D/J usage, clone statistics
├── tcrdist3/             — pairwise distance matrices, histograms
├── olga/                 — generation probability histograms per sample
├── convergence/          — TCR convergence analysis per sample
├── tcrpheno/             — TCR phenotype annotations
├── vdjdb/                — antigen-specificity matches
├── patient/              — patient-level concatenated CDR3, GIANA/GLIPH2 clusters
└── compare/              — cross-sample TCR sharing table, histogram, scatter plot
```

### Single-cell mode outputs

**Full SC sub-mode:**

```
results/
├── bridge/sc_to_bulk/              — per-sample AIRR TSVs derived from SC data
├── [bulk analysis outputs]         — same as bulk mode above, derived from SC data
├── VDJ_QC/                         — VDJ_QC_analysis.html + QC tables/figures
├── TCell_Integration/              — TCell_Integration_Report.html + Seurat RDS
├── TCRi/                           — TCRi_Report.html
├── CoNGA/                          — CoNGA_Report.html
├── GLIPH2/                         — GLIPH2_Report.html
├── TCRdist3/                       — TCRdist3_Report.html
├── GIANA/                          — GIANA_Report.html
├── Repertoire/                     — Repertoire_Report.html
├── Consensus_Clustering/           — Clonotype_Clustering_Consensus_Report.html
└── Master_Summary/                 — Master_Summary_Report.html
```

**VDJ-only sub-mode:**

```
results/
├── bridge/vdj_to_bulk/             — per-sample AIRR TSVs derived from VDJ contigs
├── [bulk analysis outputs]         — same as bulk mode above, derived from VDJ data
└── VDJ_QC/                         — VDJ_QC_analysis.html + QC tables/figures
```


---

## Architecture

The pipeline is organized into two layers:

```
main.nf                      ← detects mode + sub-mode, routes to one of three workflows
  │
  ├── workflows/bulk.nf       ← Scenario 1: TCRtoolkit bulk analysis
  ├── workflows/singlecell.nf ← Scenario 2: branches on GEX object presence
        ├── Full SC path: VDJ_QC → TCELL_INTEGRATION → Bridge 1 → TCRtoolkit + SC modules
        └── VDJ-only path: VDJ_QC → Bridge 3 → TCRtoolkit bulk analysis only
 
```

---

## Bridge Modules

The bridge modules are the only genuinely new code in this pipeline. They connect the two analysis engines.

### Bridge 1 — `SC_TO_BULK`

**Purpose:** Enables TCRtoolkit's bulk analysis on single-cell data.

**How it works:**

1. Reads the per-cell `export_cells.tsv` produced by SCRATCH-TCR's T-cell integration step
2. Extracts beta chain CDR3 (`junction_aa`), TRBV (`v_call`), TRBJ (`j_call`) from the `CTaa` and `CTgene` columns
3. Groups cells by sample, counts cells per clonotype → `duplicate_count`
4. Writes one AIRR-format TSV per sample to `bulk_samples/`
5. Writes a `synthetic_samplesheet.csv` that TCRtoolkit's input modules consume

**Input:** `TCell_Integration_Report/tables/tcr_export_cells_with_embedding.tsv`  
**Output:** `bulk_samples/<sample>_bulk.tsv` (one per sample) + `synthetic_samplesheet.csv`

### Bridge 3 — `VDJ_TO_BULK`

**Purpose:** Enables TCRtoolkit bulk analysis when no GEX Seurat object is available (VDJ-only mode).

**How it works:**

1. Reads `contigs_after_qc.tsv` produced by the VDJ QC step
2. Filters to TRB (beta chain) rows only
3. Maps columns: `cdr3` → `junction_aa`, `v_gene` → `v_call`, `j_gene` → `j_call`
4. Groups cells by sample and clonotype, counting unique barcodes → `duplicate_count`
5. Writes one AIRR-format TSV per sample to `bulk_samples/`
6. Writes a `synthetic_samplesheet.csv` for downstream TCRtoolkit modules

**Input:** `VDJ_QC/tables/contigs_after_qc.tsv`  
**Output:** `bulk_samples/<sample>_bulk.tsv` (one per sample) + `synthetic_samplesheet.csv`

**Used automatically** when `--input_annotated_object` is not provided.



---






### Bridge 2 — `CLUSTER_TO_SC` (optional)

**Purpose:** Maps TCRtoolkit's bulk GIANA/GLIPH2 cluster assignments back onto single cells for cross-pipeline annotation.

**How it works:**

1. Reads bulk cluster files (CDR3b → cluster_id) from TCRtoolkit's patient-level GIANA and GLIPH2 outputs
2. Joins onto `export_cells.tsv` using CDR3b as the key
3. Outputs per-cell TSVs with bulk cluster columns added

**When to use:** When you want to visualize or compare bulk clustering results overlaid on the single-cell UMAP. Not required for the core Consensus clustering step (which uses SCRATCH-TCR's own SC-level GLIPH2/TCRdist3/GIANA).

**Invoke manually** after a singlecell or combined run:

```bash
nextflow run main.nf \
  ... \
  --giana_clusters_file results/patient/PATIENT_giana.txt \
  --gliph_clusters_file results/patient/PATIENT_gliph2.txt
```

---

## Configuration

### Resource limits

Set global resource caps in `nextflow.config` or on the command line:

```bash
nextflow run main.nf ... --max_cpus 32 --max_memory 256.GB --max_time 48.h
```

### Institutional profile (HPC)

Create `conf/lsf.config` following the [nf-core institutional profile guide](https://nf-co.re/docs/tutorials/use_nf-core_pipelines/config_institutional_profile), then run:

```bash
nextflow run main.nf -profile lsf,singularity ...
```

### Container override

The SCRATCH-TCR container is set via `params.container` in `nextflow.config`. To use a local Singularity image:

```bash
nextflow run main.nf -profile singularity \
  --container /path/to/scratch_tcr.sif \
  ...
```

---

## Contact

For issues and contributions, please open a GitHub issue or pull request.
