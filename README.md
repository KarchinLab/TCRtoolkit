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
[!IMPORTANT]
If having an error similar to `qemu: uncaught target signal 11 (Segmentation fault) - core dumped` (observed on Mac OS Tahoe i.e. `26.5.1`), rebuild the Docker container locally, using the provided Dockerfile, and use that in the params, e.g. `--container <my image>`. It features arch selectors to pick the right `quarto` installation.