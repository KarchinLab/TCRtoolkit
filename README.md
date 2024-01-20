# Bulk TCR repertoire analysis

This repo contains the minimum viable pipeline for bulk TCR repertoire analysis.
The primary input are T cell clone counts from Bulk DNA TCR sequencing data from
Adaptive Biotechnologies. Future iterations of this pipeline will include 
options to process and analyze raw sequencing data from DNA and RNA TCRseq data. 

## Requirements

### Nextflow

This pipeline is written in Nextflow, a workflow management system. To install Nextflow, follow the instructions [here](https://www.nextflow.io/docs/latest/getstarted.html#installation).

### Docker

This pipeline uses Docker containers to run the analysis. To install Docker, follow the instructions [here](https://docs.docker.com/get-docker/).

## Running the pipeline

```
nextflow run main.nf \
    --project_name=ribas_pd1 \
    --sample_table=assets/ribas_pd1_sample_table.csv \
    --patient_table=assets/ribas_pd1_patient_table.csv \
    --output_dir=<outdir>

nextflow run main.nf \
    --project_name=Neutrophils_dominate_NSCL \
    --sample_table=assets/neutrophils_dominate_sample_table.csv \
    --patient_table=assets/neutrophils_dominate_patient_table.csv \
    --output_dir=results/out51

nextflow run main.nf \
    --project_name=Durable_PD1_Response \
    --sample_table=assets/durable_pd1_response_sample_table.csv \
    --patient_table=assets/durable_pd1_response_patient_table.csv \
    --output_dir=results/out52

nextflow run main.nf \
    --project_name=Bladder_Cancer \
    --sample_table=assets/bladder_cancer_sample_table.csv \
    --patient_table=assets/bladder_cancer_patient_table.csv \
    --output_dir=results/out53
```

Running pipeline from specific entrypoint

```
nextflow run main.nf -entry COMPARE \
    --project_name=ribas_pd1 \
    --sample_table=assets/ribas_pd1_sample_table.csv \
    --patient_table=assets/ribas_pd1_patient_table.csv \
    --output_dir=
```
```