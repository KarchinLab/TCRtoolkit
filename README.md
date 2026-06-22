# TCRtoolkit: A T Cell Repertoire Analysis Platform

![TCRtoolkit-banner](./assets/images/TCRtoolkit-banner.png)

Thanks for checking out `TCRtoolkit`, the platform for T Cell Repertoire analysis! `TCRtoolkit` is wrapped in Nextflow, written in python, and uses Docker to manage dependencies.

<p align="center">
    <img src=assets/TCR_profiling.png/>
</p>

We support bulk and single cell pseudo-bulk TCR sequencing data in either AIRR or Adaptive Biotechnologies formats.

##  Requirements

1. Nextflow

Nextflow can be used on any POSIX-compatible system (Linux, OS X, WSL). It requires Bash 3.2 (or later) and Java 11 (or later, up to 18) to be installed.

```bash
wget -qO- https://get.nextflow.io | bash
chmod +x nextflow
```

The nextflow executable is now available to run on the command line. The executable can be moved to a directory in your $PATH variable so you can run it from any directory.

2. Docker

`TCRtoolkit` runs in Docker containers available at [GHCR](https://github.com/KarchinLab/TCRtoolkit/pkgs/container/tcrtoolkit). Depending on what operating system you are running (Linux, MacOS, WSL), please refer to the [Docker documentation](https://docs.docker.com/engine/install/) for installation instructions.

## Quick Start

Below is a minimal example of how to run the pipeline. The `minimal-example` dataset provided is a small subset of the dataset from this manuscript by [Tumeh and Ribas et al. (2014)](https://www.nature.com/articles/nature13954). Note that the results are simply for demonstration purposes and are not intended for biological interpretation.

With the update to [Nextflow strict syntax](https://docs.seqera.io/nextflow/strict-syntax#using-legacy-parameter-declarations), non-default parameters should be supplied in `params.yml` rather than the command line, to ensure that number- and boolean-type parameters do not get cast as strings.

```bash
nextflow run KarchinLab/TCRtoolkit \
    -params-file params.yml
```

## Input Formats
 
`TCRtoolkit` accepts three input formats, specified via `--input_format`:
 
| Format | Description |
|---|---|
| `adaptive` | Adaptive Biotechnologies output files |
| `cellranger` | 10x Genomics CellRanger 'airr_rearrangement.tsv' output files (single-cell pseudo-bulk) |
| `airr` | AIRR-compliant tab-separated files |
 
## Workflow Levels
 
The pipeline supports multiple levels of analysis, controlled by `--workflow_level`:
 
| Level | Description |
|---|---|
| `sample` | Per-sample QC and repertoire statistics |
| `patient` | Patient-level clonotype aggregation and comparison |
| `compare` | Cross-cohort repertoire comparison and overlap |
 
Levels can be combined: `--workflow_level sample,patient,compare`
 
## HTML Reports
 
After the pipeline finishes, `TCRtoolkit` generates interactive HTML reports using [Quarto](https://quarto.org/). Four main report notebooks are rendered automatically:
 
| Notebook | Description |
|---|---|
| `template_qc.qmd` | Quality control metrics and filtering summary |
| `template_discovery_brief.qmd` | Repertoire discovery most relevant information  |
| `template_details_part1.qmd` | Detailed repertoire analysis, part 1 |
| `template_details_part2.qmd` | Detailed repertoire analysis, part 2  |
 
### Conditional Report Sections
 
Certain sub-reports are automatically appended based on input and workflow options:
 
- `--input_format cellranger` → includes single-cell phenotype report
- `--input_format adaptive` → includes bulk phenotype report
- `--workflow_level sample,patient,compare` (Patient workflow enabled) → includes patient-level clonotype analysis
- `--use_gliph2` → additionally includes GLIPH2 clustering report

## Key Parameters
 
| Parameter | Default | Description |
|---|---|---|
| `--samplesheet` | — | Path or URL to sample sheet CSV |
| `--outdir` | `out` | Output directory |
| `--input_format` | `airr` | Input format: `airr`, `adaptive`, or `cellranger` |
| `--workflow_level` | `sample,compare` | Analysis level(s): `sample`, `patient`, `compare` |
| `--use_gliph2` | `false` | Enable GLIPH2 CDR3 motif clustering |
| `--sobject_gex` | — | Path to TSV file containing cell-barcode phenotypes for pseudo-bulk phenotyping |
| `--max_memory` | `768.GB` | Maximum memory allocation |
| `--max_cpus` | `192` | Maximum CPU allocation |
 