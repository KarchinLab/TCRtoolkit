FROM mambaorg/micromamba:1.5.8

# Ensure we run as root for apt
USER root

# Update the conda base environment with required packages
COPY env.yml /tmp/env.yml
WORKDIR /tmp

RUN apt-get update && apt-get install -y \
        # runtime CLIs (KEEP)
        curl \
        wget \
        git \
        unzip \
        zip \
        jq \
        \
        # build-only deps (REMOVE LATER)
        build-essential \
        gcc \
        g++ \
    && micromamba install -y -n base -f /tmp/env.yml \
    && micromamba clean -afy \
    \
    # R packages (need compilers)
    && micromamba run -n base Rscript -e "remotes::install_github('HetzDra/turboGliph')" \
    && micromamba run -n base Rscript -e "remotes::install_github('kalaga27/tcrpheno')" \
    \
    # R cleanup
    && rm -rf /tmp/Rtmp* /root/.cache/R \
    \
    # REMOVE build deps ONLY
    && apt-get purge -y \
        build-essential \
        gcc \
        g++ \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install GIANA, patch shebang, symlink for PATH command availability
RUN git init /opt/GIANA && \
    cd /opt/GIANA && \
    git remote add origin https://github.com/s175573/GIANA.git && \
    git fetch --depth 1 origin d38aaf508c204d331f329b2f48f8b247448674bd && \
    git checkout FETCH_HEAD && \
    sed -i '1s|^#!.*|#!/usr/bin/env python3|' /opt/GIANA/GIANA4.1.py && \
    chmod +x /opt/GIANA/GIANA4.1.py && \
    ln -s /opt/GIANA/GIANA4.1.py /usr/local/bin/GIANA

# Install quarto
RUN mkdir -p /opt/quarto/1.6.42 \
    && curl -o /tmp/quarto.tar.gz -L \
        "https://github.com/quarto-dev/quarto-cli/releases/download/v1.6.42/quarto-1.6.42-linux-amd64.tar.gz" \
    && tar -zxvf quarto.tar.gz \
        -C "/opt/quarto/1.6.42" \
        --strip-components=1 \
    && rm /tmp/quarto.tar.gz 

# Install VDJmatch and symlink
RUN mkdir -p /opt/vdjmatch/1.3.1 \
    && curl -L -o vdjmatch.zip \
        "https://github.com/antigenomics/vdjmatch/releases/download/1.3.1/vdjmatch-1.3.1.zip" \
    && unzip vdjmatch.zip -d /opt/vdjmatch/1.3.1 \
    && rm vdjmatch.zip \
    && ln -s /opt/vdjmatch/1.3.1/vdjmatch-1.3.1/vdjmatch-1.3.1.jar /usr/local/bin/vdjmatch.jar

# Add to PATH
ENV PATH="/opt/quarto/1.6.42/bin:${PATH}"

# Add LD_LIBRARY_PATH for pandas
ENV LD_LIBRARY_PATH=/opt/conda/lib
