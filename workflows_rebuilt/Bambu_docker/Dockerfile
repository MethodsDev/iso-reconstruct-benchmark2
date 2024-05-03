FROM continuumio/miniconda3

RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    curl \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libpaper-utils \
    && rm -rf /var/lib/apt/lists/*

RUN conda install -c conda-forge mamba
RUN mamba install -c conda-forge r-base


# Install bambu using conda
RUN conda install -c bioconda bioconductor-bambu
RUN conda install conda-forge::r-biocmanager


# Check if bambu is installed
RUN R -e "if (!('bambu' %in% installed.packages()[,'Package'])) stop('bambu not installed')"