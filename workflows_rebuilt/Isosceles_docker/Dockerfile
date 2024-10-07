# Use a Bioconductor Docker image as a parent image
FROM bioconductor/bioconductor_docker:devel

# Set environment variable to avoid interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# Install additional dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    libseccomp-dev \
    pkg-config \
    squashfs-tools \
    cryptsetup \
    wget \
    tzdata \
    && rm -rf /var/lib/apt/lists/*

# Install necessary R packages and Isosceles
RUN R -e "install.packages(c('BiocManager', 'devtools'))" && \
    R -e "BiocManager::install(c('scran', 'scater', 'uwot', 'dittoSeq', 'DEXSeq', 'Nebulosa', 'ggbio', 'BiocStyle'))" && \
    R -e "devtools::install_github('Genentech/Isosceles', dependencies = TRUE, upgrade = TRUE, INSTALL_opts = '--install-tests')" && \
    R -e "install.packages('irlba')" && \
    R -e "devtools::install_github('powellgenomicslab/Nebulosa', upgrade = FALSE)" && \
    R -e "install.packages('optparse', repos = 'http://cran.us.r-project.org')"

# Copy the R script into the container
COPY run_isosceles.r /opt/run_isosceles.r

# Create a wrapper script
RUN echo '#!/bin/bash\nRscript /opt/run_isosceles.r "$@"' > /usr/local/bin/isosceles && \
    chmod +x /usr/local/bin/isosceles

# Set the entrypoint to bash for interactive mode
ENTRYPOINT ["/bin/bash"]