FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get -qq update
RUN apt-get -qq -y install automake
RUN apt-get -qq -y install build-essential
RUN apt-get -qq -y install curl
RUN apt-get -qq -y install libbz2-dev
RUN apt-get -qq -y install libcurl4-openssl-dev
RUN apt-get -qq -y install liblzma-dev
RUN apt-get -qq -y install libncurses-dev
RUN apt-get -qq -y install libncurses5-dev
RUN apt-get -qq -y install python3
RUN apt-get -qq -y install python3-dev
RUN apt-get -qq -y install python3-distutils
RUN apt-get -qq -y install python3-pip
RUN apt-get -qq -y install python3-setuptools
RUN apt-get -qq -y install wget
RUN apt-get -qq -y install zlib1g-dev
RUN apt-get -qq -y install git

ENV SRC=/usr/local/src
ENV BIN=/usr/local/bin

WORKDIR $SRC


RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
	tar jxf samtools-1.9.tar.bz2 && \
	rm samtools-1.9.tar.bz2 && \
	cd samtools-1.9 && \
	./configure --prefix $(pwd) && \
	make

ENV PATH=${PATH}:/usr/local/src/samtools-1.9


#RUN wget https://github.com/christopher-vollmers/Mandalorion/archive/refs/tags/v4.3.1.tar.gz && \
#    tar xvf v4.3.1.tar.gz && \
#    mv Mandalorion-4.3.1 Mandalorion && \
#    cd Mandalorion && \
#	chmod +x setup.sh && \
#	sh setup.sh && \
#	/usr/local/src/Mandalorion/Mando.py -h

#RUN git clone https://github.com/christopher-vollmers/Mandalorion.git && \

RUN git clone https://github.com/brianjohnhaas/Mandalorion.git && \
    cd Mandalorion && \
	chmod +x setup.sh && \
	sh setup.sh && \
	/usr/local/src/Mandalorion/Mando.py -h

COPY Mandalorian-runner.py /usr/local/bin/
