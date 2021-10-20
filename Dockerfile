FROM rocker/tidyverse:4.1.0
LABEL version="1.01"
LABEL software="HTGTS"
MAINTAINER Sergej Andrejev <sandrejev@gmail.com>

ARG http_proxy
ARG https_proxy
ENV DESTINATION=/bin/transloc_pipeline
ENV MAKEFLAGS='-j 32'
ENV PATH=$DESTINATION:$DESTINATION/HTGTS/bin:$DESTINATION/HTGTS/R:$DESTINATION/HTGTS/tools:${PATH}
ENV GENOME_DB=./genomes
ENV BOWTIE2_INDEXES=./genomes
ENV PASSWORD=s215v
ENV PERL_MM_USE_DEFAULT=1
ARG DEBIAN_FRONTEND=noninteractive
WORKDIR /mount

RUN mkdir $DESTINATION

RUN http_proxy="${http_proxy}" https_proxy="${https_proxy}" apt-get update
RUN http_proxy="${http_proxy}" https_proxy="${https_proxy}" apt-get install -y \
    python3-dev \
    libgsl-dev \
    bowtie2 \
    python3-tqdm \
    python3-pybedtools \
    python3-pandas \
    python3-pip \
    python3-numpy \
    bioperl \
    libbam-dev

RUN http_proxy="${http_proxy}" https_proxy="${https_proxy}" install2.r --error --deps TRUE argparser plyr data.table bedr ggplot2 forcats VennDiagram sortable rstatix
RUN http_proxy="${http_proxy}" https_proxy="${https_proxy}" /usr/local/lib/R/site-library/littler/examples/installBioc.r GenomicRanges BSgenome

# MACS2
RUN http_proxy="${http_proxy}" https_proxy="${https_proxy}" python3.8 -m pip install macs2

#
## Install EA-utils
RUN http_proxy="${http_proxy}" https_proxy="${https_proxy}" wget https://storage.googleapis.com/google-code-archive-downloads/v2/code.google.com/ea-utils/ea-utils.1.1.2-537.tar.gz && \
    tar xzvf ea-utils.1.1.2-537.tar.gz && \
    cd ea-utils.1.1.2-537 && \
    make && \
    PREFIX=/usr make install

#
# Install old version of SAMTOOLS 0.1.18. Either use UBUNTU libbam-dev or a version compiled from EA-utils
#
ENV SAMTOOLS=/usr/include/samtools
RUN ln /usr/lib/libbam.a -s /usr/include/samtools/libbam.a
#RUN mkdir /usr/include/samtools && \
#    cp samtools/libbam.a /usr/include/samtools/libbam.a && \
#    cp samtools/bam.h /usr/include/samtools/bam.h

# Install Perl modules.
RUN export PERL_MM_USE_DEFAULT=1
RUN http_proxy="${http_proxy}" https_proxy="${https_proxy}" cpan install Getopt::Long Text::CSV List::MoreUtils Data::GUID Interpolation IPC::System::Simple Storable Switch::Plain File::Which
RUN http_proxy="${http_proxy}" https_proxy="${https_proxy}" cpan install Bio::DB::Sam

## Install SeqPrep
RUN http_proxy="${http_proxy}" https_proxy="${https_proxy}" git clone https://github.com/jstjohn/SeqPrep.git && \
    cd SeqPrep && \
    make && \
    HOME=/usr make install && \
    ln /usr/bin/SeqPrep -s /usr/bin/seqprep
#
## Fastq-multx
#RUN cd /bin && \
#    http_proxy="${http_proxy}" https_proxy="${https_proxy}" git clone https://github.com/brwnj/fastq-multx.git fastq-multx_source && \
#    cd fastq-multx_source && \
#    make
#
COPY ./transloc_pipeline_v3 $DESTINATION/HTGTS
COPY Entrypoint  $DESTINATION/Entrypoint
RUN chmod 755 $DESTINATION/Entrypoint $DESTINATION/HTGTS/bin/*  $DESTINATION/HTGTS/tools/* $DESTINATION/HTGTS/R/*

#ENTRYPOINT ["Entrypoint"]