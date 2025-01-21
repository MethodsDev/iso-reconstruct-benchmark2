FROM bioconductor/bioconductor_docker:RELEASE_3_20

RUN Rscript -e 'BiocManager::install("bambu")'


# Check if bambu is installed
RUN R -e "if (!('bambu' %in% installed.packages()[,'Package'])) stop('bambu not installed')"

RUN Rscript -e 'BiocManager::install(c("argparse"))'

COPY bambu-runner.Rscript /usr/local/bin/
