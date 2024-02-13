# Bulk TCR repertoire analysis

Thanks for checking out the `tcr-toolkit`!

This repo contains the minimum viable pipeline for bulk TCR repertoire analysis.
The primary input are T cell clone counts from Bulk DNA TCR sequencing data from
Adaptive Biotechnologies. Future iterations of this pipeline will include 
options to process and analyze raw sequencing data from DNA and RNA TCRseq data. 

##  Requirements
* Nextflow
* Docker
* bash command line

## Installation

```{python}
## code here
```

## Running the pipeline

```
nextflow run main.nf \
    --project_name=ribas_pd1 \
    --sample_table=/lab/projects1/btc/bulk-tcrseq/assets/ribas_pd1_sample_table.csv \
    --patient_table=/lab/projects1/btc/bulk-tcrseq/assets/ribas_pd1_patient_table.csv \
    --output_dir=<outdir>