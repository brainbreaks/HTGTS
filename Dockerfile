FROM rocker/tidyverse:4.1.0
LABEL version="1.01"
LABEL software="HTGTS"
MAINTAINER Sergej Andrejev <sandrejev@gmail.com>

ARG http_proxy
ARG https_proxy
ENV DESTINATION=/bin/transloc_pipeline
ENV MAKEFLAGS='-j 32'
ENV PATH=$DESTINATION:$DESTINATION/HTGTS/bin:$DESTINATION/HTGTS/R:${PATH}
ENV GENOME_DB=./genomes
ENV BOWTIE2_INDEXES=~/mount
ENV PASSWORD=s215v
ENV PERL_MM_USE_DEFAULT=1
ARG DEBIAN_FRONTEND=noninteractive
WORKDIR /mount

RUN mkdir $DESTINATION

RUN http_proxy="${http_proxy}" https_proxy="${https_proxy}" apt-get update
RUN http_proxy="${http_proxy}" https_proxy="${https_proxy}" apt-get install -y bioperl libbio-samtools-perl libtext-csv-perl libfile-which-perl libipc-system-simple-perl ea-utils seqprep bowtie2 python3-tqdm python3-pybedtools python3-pandas

RUN http_proxy="${http_proxy}" https_proxy="${https_proxy}" install2.r --error --deps TRUE argparser plyr
RUN http_proxy="${http_proxy}" https_proxy="${https_proxy}" /usr/local/lib/R/site-library/littler/examples/installBioc.r GenomicRanges BSgenome

# HTGTS v2
# RUN cd $DESTINATION && git clone https://github.com/robinmeyers/transloc_pipeline.git HTGTS

RUN http_proxy="${http_proxy}" https_proxy="${https_proxy}" cpan install Interpolation Storable Switch::Plain
RUN cd /bin && \
    http_proxy="${http_proxy}" https_proxy="${https_proxy}" git clone https://github.com/brwnj/fastq-multx.git fastq-multx_source && \
    cd fastq-multx_source && \
    make

COPY ./transloc_pipeline_v3 $DESTINATION/HTGTS
COPY preprocess/download.py  $DESTINATION/download.py
RUN echo '#!/bin/bash \npython3 $DESTINATION/download.py "$@"' > $DESTINATION/download; chmod 755 $DESTINATION/download
COPY preprocess/gff_longest_transcript.py  $DESTINATION/gff_longest_transcript.py
RUN echo '#!/bin/bash \npython3 $DESTINATION/gff_longest_transcript.py "$@"' > $DESTINATION/longest-transcript; chmod 755 $DESTINATION/longest-transcript
COPY Entrypoint  $DESTINATION/Entrypoint

#ENTRYPOINT ["Entrypoint"]