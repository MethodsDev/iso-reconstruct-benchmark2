# Use miniconda3 as base image
FROM continuumio/miniconda3

# Set the working directory
WORKDIR /usr/src/app

# Install necessary packages in one command
RUN conda install -c defaults -c bioconda -c conda-forge \
    perl=5.26.2 \
    libgcc-ng \
    hmmer \
    samtools

# Download and install ucsc-faToTwoBit from source
RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit && \
    chmod +x faToTwoBit && \
    mv faToTwoBit /usr/local/bin/

# Update GPG keys and install necessary tools for building from source
RUN apt-get clean && \
    apt-get update && \
    apt-get install -y gnupg build-essential wget git && \
    apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 04EE7237B7D453EC 648ACFD622F3D138

# Download and install ucsc-bedgraphtobigwig from source
RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig && \
    chmod +x bedGraphToBigWig && \
    mv bedGraphToBigWig /usr/local/bin/

# Download and install minimap2 from source
RUN wget https://github.com/lh3/minimap2/releases/download/v2.22/minimap2-2.22_x64-linux.tar.bz2 && \
    tar -xvf minimap2-2.22_x64-linux.tar.bz2 && \
    mv minimap2-2.22_x64-linux/minimap2 /usr/local/bin/

# Download and install a specific version of BLAST
RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.11.0/ncbi-blast-2.11.0+-x64-linux.tar.gz -O blast.tar.gz && \
    tar -zxvf blast.tar.gz && \
    mv ncbi-blast-2.11.0+ /usr/local/ && \
    ln -s /usr/local/ncbi-blast-2.11.0+/bin/* /usr/local/bin/

# Clone the Espresso repository
RUN git clone https://github.com/Xinglab/espresso.git

# Set the appropriate permissions
RUN chmod -R 755 espresso

# Set PATH
ENV PATH /usr/local/bin:$PATH