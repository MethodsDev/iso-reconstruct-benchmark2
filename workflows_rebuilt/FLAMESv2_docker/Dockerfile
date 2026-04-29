FROM ghcr.io/mritchielab/flames:devel

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get -y install --no-install-recommends \
    samtools \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'install.packages("argparse", repos="https://cran.rstudio.com")'

COPY FLAMESv2-runner.Rscript /usr/local/bin/

WORKDIR /workspace

ENTRYPOINT ["FLAMESv2-runner.Rscript"]
