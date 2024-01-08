FROM condaforge/miniforge3:23.3.1-1

# copy the environment file into /tmp
COPY env.yml /tmp/env.yml

# update the mamba base environment with required packages
WORKDIR /tmp
RUN mamba env update -n base --file env.yml

# install curl
RUN apt-get update \
    && apt-get install -y curl \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# install quarto 
RUN mkdir -p /opt/quarto/1.3.433 \
    && curl -o quarto.tar.gz -L \
        "https://github.com/quarto-dev/quarto-cli/releases/download/v1.3.433/quarto-1.3.433-linux-amd64.tar.gz" \
    && tar -zxvf quarto.tar.gz \
        -C "/opt/quarto/1.3.433" \
        --strip-components=1 \
    && rm quarto.tar.gz 

# install R
RUN mamba install -y r-base=4.1.0 \
    && mamba clean -afy

# add quarto to the path
ENV PATH="/opt/quarto/1.3.433/bin:${PATH}"
