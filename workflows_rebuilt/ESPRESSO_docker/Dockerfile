# Use miniconda3 as base image
FROM continuumio/miniconda3

# Install mamba
RUN conda install -c conda-forge mamba -y && conda clean -a -y

# Create a new environment and install ESPRESSO
RUN mamba create -n espresso_env -c bioconda -c conda-forge espresso=1.5.0-0 && \
    conda clean -a -y

# Activate the environment and set it as the default
SHELL ["conda", "run", "-n", "espresso_env", "/bin/bash", "-c"]

# Set the working directory
WORKDIR /usr/src/app

# Copy the ESPRESSO repository (if needed)
# RUN git clone https://github.com/Xinglab/espresso.git

# Set the PATH to include the conda environment
ENV PATH /opt/conda/envs/espresso_env/bin:$PATH

# Verify installations
RUN perl -e 'use threads; print("Perl threading support: ok\n")' && \
    perl -e 'use Storable; print("Storable version: $Storable::VERSION\n")' && \
    hmmscan -h && \
    blastn -version && \
    samtools --version

COPY ESPRESSO-runner.py /usr/local/bin

