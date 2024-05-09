# Start from the sqanti3 image
FROM us-central1-docker.pkg.dev/methods-dev-lab/lrtools-sqanti3/lrtools-sqanti3-plus

# Set the working directory in the container to /app
WORKDIR /app

# Install necessary packages
RUN apt-get update --fix-missing && apt-get install -y \
    python3.8 \
    python3-pip \
    wget \
    bzip2 \
    gcc \
    libncurses5-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev

# Install Python libraries
RUN python3.8 -m pip install --no-cache-dir pandas numpy matplotlib matplotlib_venn

# Copy the specific file into the container at /app
COPY LRQuant.py /app

# Download and install minimap2 from source
RUN wget https://github.com/lh3/minimap2/releases/download/v2.22/minimap2-2.22_x64-linux.tar.bz2 && \
    tar -xvf minimap2-2.22_x64-linux.tar.bz2 && \
    mv minimap2-2.22_x64-linux/minimap2 /usr/local/bin/

# Install gffcompare
RUN wget http://ccb.jhu.edu/software/stringtie/dl/gffcompare-0.11.2.Linux_x86_64.tar.gz && \
    tar -xvzf gffcompare-0.11.2.Linux_x86_64.tar.gz && \
    cp gffcompare-0.11.2.Linux_x86_64/gffcompare /usr/local/bin/

# Install minimap2
RUN wget https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17_x64-linux.tar.bz2 && \
    tar -xvf minimap2-2.17_x64-linux.tar.bz2 && \
    cp minimap2-2.17_x64-linux/minimap2 /usr/local/bin/

# Install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
    tar jxf samtools-1.9.tar.bz2 && \
    rm samtools-1.9.tar.bz2 && \
    cd samtools-1.9 && \
    ./configure --prefix $(pwd) && \
    make

# Install gffread
RUN wget http://ccb.jhu.edu/software/stringtie/dl/gffread-0.12.1.tar.gz
RUN tar -zxvf gffread-0.12.1.tar.gz
WORKDIR /app/gffread-0.12.1
RUN make

# Go back to the main working directory
WORKDIR /app

# Add installed tools to PATH
ENV PATH="/usr/local/bin:/app/gffread-0.12.1:${PATH}"

# Create a shell script to run LRQuant.py with python3.8
RUN echo '#!/bin/sh' > /usr/local/bin/LRQuant && \
    echo 'python3.8 /app/LRQuant.py "$@"' >> /usr/local/bin/LRQuant && \
    chmod +x /usr/local/bin/LRQuant