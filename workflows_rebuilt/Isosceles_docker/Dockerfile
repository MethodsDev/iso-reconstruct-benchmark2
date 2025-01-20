FROM bioconductor/bioconductor_docker:RELEASE_3_15

RUN Rscript -e 'BiocManager::install(c("scran", "scater", "uwot", "dittoSeq", "DEXSeq", "Nebulosa", "slingshot", "pdftools", "ggbio", "ggnewscale", "BiocStyle", "fastmatch", "RcppArmadillo", "tidyverse"), update = TRUE, ask = FALSE)'

RUN Rscript -e 'install.packages(c("poolr", "irlba", "ggrastr"), repos = "https://cran.r-project.org")'

RUN Rscript -e 'devtools::install_github("powellgenomicslab/Nebulosa", upgrade = FALSE)'

RUN Rscript -e 'devtools::install_github("Genentech/Isosceles", ref = "0.2.0", dependencies = TRUE, upgrade = TRUE, INSTALL_opts = "--install-tests")'

RUN Rscript -e 'BiocManager::install(c("argparse"))'


COPY Isosceles-runner.Rscript /usr/local/bin

ENTRYPOINT ["/usr/local/bin/Isosceles-runner.Rscript"]

