# Use miniconda3 as base image
FROM continuumio/miniconda3

# Install mamba
RUN conda install -c conda-forge mamba -y && conda clean -a -y

# Create a new environment and install FLAIR, samtools, and bam2Bed12
RUN mamba create -n flair_env -c bioconda -c conda-forge flair=2.0.0 samtools && \
    conda clean -a -y

# Activate the environment and set it as the default
SHELL ["conda", "run", "-n", "flair_env", "/bin/bash", "-c"]

# Set the working directory
WORKDIR /usr/src/app

# Set the PATH to include the conda environment
ENV PATH /opt/conda/envs/flair_env/bin:$PATH

# Verify installations
RUN flair --help && \
    samtools --version && \
    bam2Bed12 --help


WORKDIR /usr/local/src
RUN wget https://github.com/gpertea/gffread/releases/download/v0.12.7/gffread-0.12.7.Linux_x86_64.tar.gz && \
    tar xvf gffread-0.12.7.Linux_x86_64.tar.gz && \
    mv gffread-0.12.7.Linux_x86_64/gffread /usr/local/bin


COPY flair-runner.py /usr/local/bin

