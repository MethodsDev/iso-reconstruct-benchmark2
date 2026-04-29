FROM bioconductor/bioconductor_docker:RELEASE_3_22

# Install runtime dependency, pin bambu to v3.12.1, and verify in one layer
RUN Rscript -e 'BiocManager::install(c("argparse"), ask = FALSE, update = FALSE)' \
 && Rscript -e 'install.packages("remotes", repos = "https://cloud.r-project.org")' \
 && Rscript -e 'remotes::install_github("GoekeLab/bambu@v3.12.1", upgrade = "never", dependencies = TRUE)' \
 && Rscript -e 'stopifnot(as.character(packageVersion("bambu")) == "3.12.1")'

COPY bambu-runner.Rscript /usr/local/bin/
