# TCRtoolkit: A T Cell Repertoire Analysis Platform

![TCRtoolkit-banner](./assets/TCRtoolkit-banner.png)

Thanks for checking out `TCRtoolkit`, the platform for T Cell Repertoire analysis! `TCRtoolkit` is wrapped in NextFlow, written in python, and uses Docker to manage dependencies.

This platform is designed to be a flexible, scalable, and easy-to-use tool for analyzing TCR sequencing data. Because `TCRtoolkit` is built on NextFlow, it can be run via the command line locally, on HPC environments, and even cloud-based systems like AWS, GCP, and Azure. In addition, we plan to provide a web-based interface for running the pipeline with no coding required.

We currently support bulk TCRseq data from Adaptive Biotechnologies, but plan to add single cell and spatial TCRseq datatypes in the near future.

##  Requirements

### 1. Nextflow

Nextflow can be used on any POSIX-compatible system (Linux, OS X, WSL). It requires Bash 3.2 (or later) and Java 11 (or later, up to 18) to be installed.

```{bash}
wget -qO- https://get.nextflow.io | bash
chmod +x nextflow
```

The nextflow executable is now available to run on the command line. The executable can be moved to a directory in your $PATH variable so you can run it from any directory, for example: `mv nextflow /usr/local/bin`.

### A note about Docker

Docker is a platform to enable applications to run in a consistent environment across different computing environments. Docker images will be automatically downloaded when the pipeline is run, no need for user installation. We plan to add support for Singularity containers in the near future.

## Running the pipeline

An example execution of `TCRtoolkit` below. 

```{bash}
## clone the repository
git clone https://github.com/KarchinLab/TCRtoolkit.git
cd TCRtoolkit

nextflow run main.nf \
    --project_name=Bladder_Cancer \
    --sample_table=/lab/projects1/btc/tcr-toolkit/assets/bladder_cancer_sample_table.csv \
    --patient_table=/lab/projects1/btc/tcr-toolkit/assets/bladder_cancer_patient_table.csv \ 
    --output_dir=results/urothelial_cancer_PD1_blockade_v3 \
    --data_dir=/lab/data/btc-tcr/adaptive-data-bulk/public/Contribution_of_systemic_and_somatic_factors_to_clinical_response_and_resistance_to_PD-L1_blockade_in_urothelial_cancer_An_exploratory_multi-omic_analysis
```
