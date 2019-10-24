from continuumio/anaconda3:latest
RUN apt-get -y update && apt-get -y install build-essential
RUN apt-get install -y --no-install-recommends apt-utils

# install essential packages for pysam
RUN apt-get -y install zlib1g-dev && apt-get -y install libbz2-dev
RUN apt-get -y install libcurl4-gnutls-dev && apt-get -y install libssl-dev

# install methylpy, bowtie/bowtie2 and other methylpy dependencies
RUN pip install methylpy
RUN conda install -y -c bioconda bowtie && conda install -y -c bioconda bowtie2
RUN pip install cutadapt && conda install -y -c bioconda samtools
RUN mkdir -p /usr/share/man/man1
RUN apt-get -y install openjdk-11-jdk


# install dependencies for DMRfind
RUN apt-get -y install libgsl-dev
RUN ln -s /usr/lib/x86_64-linux-gnu/libgsl.so.23 lib/libgsl.so.0
