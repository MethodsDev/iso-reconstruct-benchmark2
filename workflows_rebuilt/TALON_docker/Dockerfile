# Specify the base image
FROM python:3.7

# Update package lists
RUN apt-get update

# Install dependencies
RUN apt-get -qq -y install python3
RUN apt-get -qq -y install python3-dev
RUN apt-get -qq -y install python3-distutils
RUN apt-get -qq -y install python3-pip
RUN apt-get -qq -y install python3-setuptools
RUN apt-get -qq -y install wget
RUN apt-get -qq -y install zlib1g-dev
RUN apt-get -qq -y install unzip

# Install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2 && \
    tar xvf samtools-1.16.1.tar.bz2 && \
    cd samtools-1.16.1 && \
    ./configure && \
    make && \
    make install

# Link python3 to python
RUN ln -sf /usr/bin/python3 /usr/bin/python

# Install Cython and pytest
RUN pip install cython pytest

# Install bedtools
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary && \
    chmod 755 bedtools.static.binary && \
    mv bedtools.static.binary /usr/local/bin/bedtools

# Install TALON
RUN wget https://github.com/mortazavilab/TALON/archive/refs/heads/master.zip && \
    unzip master.zip && \
    cd TALON-master && \
    pip install .

# Keep the container running and provide a shell
CMD ["bash"]