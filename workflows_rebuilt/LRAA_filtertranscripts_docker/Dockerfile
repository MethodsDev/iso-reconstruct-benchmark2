FROM condaforge/mambaforge

ENV SRC=/usr/local/src
ENV BIN=/usr/local/bin

WORKDIR $SRC

# Install Biopython (and Python, implicitly) using mamba from the bioconda channel
RUN mamba install -c conda-forge -c bioconda biopython

# Copy the local FilterTranscripts.py file to the container
COPY FilterTranscripts.py $SRC/FilterTranscripts.py

# Make the script executable
RUN chmod +x $SRC/FilterTranscripts.py

# Create a symbolic link to make the script callable with 'FilterTranscripts'
RUN ln -s $SRC/FilterTranscripts.py $BIN/FilterTranscripts

# Use ENTRYPOINT to specify the executable script
ENTRYPOINT ["python", "/usr/local/src/FilterTranscripts.py"]