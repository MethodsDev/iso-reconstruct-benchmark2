# Use the official Rust nightly image as a base
FROM rustlang/rust:nightly

# Set the current working directory in the container
WORKDIR /usr/src/oarfish

# Set environment variables
ENV SRC=/usr/local/src
ENV BIN=/usr/local/bin
ENV PATH="/usr/src/oarfish/target/release:/opt/conda/bin:${PATH}"

# Clone the Oarfish repository
RUN git clone https://github.com/COMBINE-lab/oarfish.git .

# Build the project
RUN cargo +nightly build --release

# Download and install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda && \
    rm Miniconda3-latest-Linux-x86_64.sh

# Download and install minimap2 from source
RUN wget -qO- https://github.com/lh3/minimap2/releases/download/v2.22/minimap2-2.22_x64-linux.tar.bz2 | tar -xj && \
    mv minimap2-2.22_x64-linux/minimap2 /usr/local/bin/
    
# samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2 && \
	tar xvf samtools-1.16.1.tar.bz2 && \
	cd samtools-1.16.1 && \
	./configure && \
	make && \
	make install

# Install gffread and samtools using conda
RUN conda install -y -c bioconda gffread
RUN conda clean -a