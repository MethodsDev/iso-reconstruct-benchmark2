FROM condaforge/mambaforge

ENV SRC=/usr/local/src
ENV BIN=/usr/local/bin

WORKDIR $SRC

RUN mamba install -c bioconda samtools
RUN mamba install -c bioconda isoseq
RUN mamba install -c bioconda pbmm2
RUN mamba install -c bioconda pbpigeon

COPY  IsoSeq-runner.py /usr/local/bin/

