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

RUN R -e "install.packages('BiocManager', repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install('FLAMES')"

ENTRYPOINT ["Rscript"]