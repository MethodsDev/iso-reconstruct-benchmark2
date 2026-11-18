FROM condaforge/mambaforge

ENV SRC=/usr/local/src
ENV BIN=/usr/local/bin

WORKDIR $SRC

RUN mamba install -y -c bioconda samtools=1.21
RUN mamba install -y -c bioconda isoseq=4.2.0
RUN mamba install -y -c bioconda pbmm2=1.16.99
RUN mamba install -y -c bioconda pbpigeon=1.3.0

COPY  IsoSeq-runner.py /usr/local/bin/

