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
> [!IMPORTANT]
> If having an error similar to `qemu: uncaught target signal 11 (Segmentation fault) - core dumped` (observed on Mac OS Tahoe i.e. `26.5.1`), rebuild the Docker container locally, using the provided Dockerfile, and use that in the params, e.g. `--container <my image>`. It features arch selectors to pick the right `quarto` installation.

## Input Formats
 
`TCRtoolkit` accepts two input formats, specified via `--input_format`:
 
| Format | Description |
|---|---|
| `adaptive` | Adaptive Biotechnologies output files |
| `airr` | AIRR-compliant tab-separated files |
 
## Workflow Levels
 
The pipeline supports multiple levels of analysis, controlled by `--workflow_level`:
 
| Level | Description |
|---|---|
| `sample` | Per-sample QC and repertoire statistics |
| `patient` | Patient-level clonotype clustering: GIANA, GLIPH2 (`--use_gliph2`) and TCRdist3, each over all of a patient's samples |
| `compare` | Cross-cohort repertoire comparison and overlap |
 
Levels can be combined: `--workflow_level sample,patient,compare`
 
## HTML Reports
 
After the pipeline finishes, `TCRtoolkit` generates interactive HTML reports using [Quarto](https://quarto.org/). Report notebooks are rendered automatically depending on which `--workflow_level` stages ran:
 
| Notebook | Rendered when | Description |
|---|---|---|
| `template_qc.qmd` | `sample` | Quality control metrics and filtering summary |
| `template_discovery_brief.qmd` | `sample` + `compare` | Highest-impact repertoire findings |
| `template_details_sample.qmd` | `sample` | Detailed per-sample repertoire analysis |
| `template_details_compare.qmd` | `sample` + `compare` | Detailed cross-sample / longitudinal comparison analysis |
| `template_details_patient.qmd` | `sample` + `patient` | Patient-level clonotype clustering (GIANA, and GLIPH2 when `--use_gliph2`) |
 
### Conditional Report Sections
 
Certain sub-reports are automatically appended based on input and workflow options:
 
- TCRpheno (phenotype predicted from the TCR β sequence) runs on bulk data and the single-cell VDJ-only route, and its section is added to the discovery brief. It is skipped on the single-cell VDJ+GEX route, where the GEX object gives each cell's measured phenotype.
- The per-sample TCRdist3 sections read each sample's block of its patient's TCRdist3 matrix, so they need `patient` in `--workflow_level`.
- `--workflow_level` includes `sample` and `patient` → renders `template_details_patient.qmd` (patient-level clonotype analysis)
- `--use_gliph2` → additionally includes GLIPH2 clustering in the patient details report
- `--mode singlecell` → the same five notebooks are rendered, each with a single-cell section appended (`notebooks/_sc_*.qmd`): VDJ / GEX-matching QC, cell states of expanded and dynamically changing clones, CoNGA and TCRi, blood vs tissue, clone tracking over time, and cell-level consensus clusters. The report samplesheet is built from `--sample_sheet` (`patient_id`, `timepoint`, `origin`/`tissue` columns).

### Repertoire metrics (both modalities)

Diversity, clone-size bins, overlap and clone sharing are computed once, by `REPERTOIRE_METRICS` (`bin/repertoire_metrics.py`) in the `sample` stage, and every report notebook and the Master Summary reads its tables instead of recomputing them:

| Output | Contents |
|---|---|
| `sample/sample_stats.csv` | Per-sample stats with diversity: `num_clones`, `num_TCRs`, `clonality`, `simpson_index[_corrected]`, `shannon_entropy`, `inverse_simpson`, `gini_coefficient`, `hill_q0`–`hill_q2`, `clone_unit` |
| `repertoire/diversity_by_group.tsv` | The same metrics on clones pooled per patient / condition / timepoint / origin |
| `repertoire/clone_size_bins.tsv` | Clones and TCRs per frequency bin (`--clone_bin_breaks`, `--clone_bin_labels`) |
| `repertoire/overlap_samples.tsv` | Every sample pair: shared clones, Jaccard, Sørensen, Morisita-Horn (within-patient pairs included) |
| `repertoire/overlap_patients.tsv` | Every patient pair (clones pooled per patient): shared clones, Jaccard |
| `repertoire/clone_sharing.tsv` | Per clone: samples, patients, conditions and timepoints it occurs in |
| `repertoire/repertoire_summary.tsv` | Cohort rollup |

The clone unit depends on the modality: bulk counts TRB clonotypes (reads/templates), and single-cell counts paired clones (`A:<CDR3a>|B:<CDR3b>`, cells), taken from `paired_clones.tsv` written by the route's bridge (`SC_TO_CDR3` or `VDJ_TO_BULK`), which define the clone identically. The public-clone analysis in `tcrsharing` stays on TRB CDR3 because generation probability (OLGA) is defined per CDR3.

### Clonotype clustering and Master Summary (both modalities)

After the cohort comparison, both modalities run the same chain: a per-cell (single-cell) or per-clonotype (bulk) export → `CLUSTER_TO_SC` (GIANA / GLIPH2 / TCRdist3 clusters mapped onto those rows) → `CONSENSUS_CLUSTERING` → `CLUSTER_ROLLUP` → `MASTER_SUMMARY`. Without a GEX object these steps run on the export table alone. Bulk runs need `compare` plus `sample` and/or `patient`; disable with `--run_clonotype_clustering false`.

All three methods cluster per patient, over all of that patient's samples, and their clusters are mapped onto the rows of that patient only (matched on patient and CDR3β). TCRdist3 clusters are the connected components of the patient's distance matrix at `--tcrdist_radius`, with at least two distinct CDR3β. Because the three methods label clusters independently, the consensus does not compare labels: two clonotypes are linked when at least `--consensus_min_methods` methods (default 2) place them in the same cluster, and each connected group of linked clonotypes is one consensus cluster. `consensus_cluster_map.tsv` lists, per consensus cluster, its clonotypes, cells, the methods whose agreement built it (`consensus_methods`) and the member clusters of each method.

With a GEX object (`--input_annotated_object`), `MASTER_SUMMARY` also renders `SC_GEX_Report.html` from `modules/scratch/MASTER_SUMMARY/_sc_gex.qmd` (cell-state map, composition by sample and condition, expansion and clone sharing by cell state, marker-gene expression of expanded cells, CoNGA, TCRi and consensus clusters by cell state), and the same section is appended to the Master Summary. Other routes get a short stub in its place.

## Key Parameters
 
| Parameter | Default | Description |
|---|---|---|
| `--samplesheet` | — | Path or URL to sample sheet CSV |
| `--outdir` | `out` | Output directory |
| `--input_format` | `airr` | Input format: `airr` or `adaptive` |
| `--workflow_level` | `sample,compare` | Analysis level(s): `sample`, `patient`, `compare` |
| `--use_gliph2` | `false` | Enable GLIPH2 CDR3 motif clustering |
| `--max_memory` | `768.GB` | Maximum memory allocation |
| `--max_cpus` | `192` | Maximum CPU allocation |
 
