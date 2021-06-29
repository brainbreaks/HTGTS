HTGTS pipeline package (Wei, Pei-Chi group; DKFZ)
====================================================
High-throughput, genome-wide, translocation sequencing (HTGTS) pipeline. More documentation about the pipeline itself can be found [here](https://robinmeyers.github.io/transloc_pipeline)


Table of Contents
=================

  * Using HTGTS docker image
    * [Singularity](#singularity) 
      * [Pull HTGTS image from DockerHUB](#singularity-pull)
      * [Run container](#singularity-run)
      * [Inspect container](#singularity-inspect)
    * [Docker](#docker) 
      * [Pull HTGTS image from DockerHUB](#docker-pull)
      * [Run container](#docker-run)
      * [Inspect container](#docker-inspect)
    * [Build HTGTS docker image](#build) 
      * [Build docker image](#build-build)
      * [Push docker image to Docker HUB](#build-build)
      * [Convert cached docker image to singularity (for local testing)](#build-convert)

<a name="singularity">Singularity</a>
====================================================

<a name="singularity-pull">Pull HTGTS image from data server</a>
----------------------------------------------------
```console
singularity pull docker://sandrejev/htgts:latest
```

<a name="singularity-run">Run container</a>
----------------------------------------------------
Before running HTGTS pipeline you will need to obtain genome (fasta), bowtie index (*.bt2), chromosome sizes (*.chrom.sizes), annotation (ChromInfo.txt, cytoBand.txt, refGene.bed). For mm9, mm10 and hg19 these can be downloaded automatically with `download` command. 
The default location is to download these files to `./genomes` folder in the current directory, but this can be changed using `GENOME_DB` environment variable. You can either set this variable in your `.bashrc` or `.bash_profile` file or prepend all commands in this 
tutorial like this `GENOME_DB=./new_location singularity [command]`
```console
singularity exec -B `pwd` htgts_latest.sif download mm10
```

Run HTGTS pipeline. Keep in mind that refGene annotation file (refGene.bed) is created automatically for you from geneRef.gtf 
```console
singularity exec -B `pwd` htgts_latest.sif TranslocPreprocess.pl tutorial_metadata.txt preprocess --read1 pooled_R1.fastq.gz --read2 pooled_R2.fastq.gz
singularity exec -B `pwd` htgts_latest.sif TranslocWrapper.pl tutorial_metadata.txt preprocess/ results/ --threads 2
```

Detect translocation peaks using MACS2
```console
singularity exec -B `pwd` htgts_latest.sif tlx2BED-MACS.pl
singularity exec -B `pwd` htgts_latest.sif macs2
```

<a name="singularity-inspect">Inspect container</a>
----------------------------------------------------
```console
singularity shell htgts_latest.sif
```

<a name="docker">Docker</a>
====================================================

<a name="docker-pull">Pull HTGTS image from DockerHUB</a>
----------------------------------------------------
```console
docker pull sandrejev:htgts
```

<a name="docker-run">Run container</a>
----------------------------------------------------
Before running HTGTS pipeline you will need to obtain genome(fasta), bowtie index (*.bt2), chromosome sizes (*.chrom.sizes) and annotation. For mm9, mm10 and hg19 these can be downloaded automatically with `download` command. The only exception being
annotation *.bed file.
```console
docker run -v ${PWD}:/mount -u $(id -g ${USER}):$(id -g ${USER}) -it --entrypoint download htgts mm10
```

Run HTGTS pipeline. Keep in mind that refGene annotation file (refGene.bed) is created automatically for you from geneRef.gtf 
```console
docker run -v ${PWD}:/mount -u $(id -g ${USER}):$(id -g ${USER}) -it --entrypoint TranslocPreprocess.pl htgts tutorial_metadata.txt preprocess --read1 pooled_R1.fastq.gz --read2 pooled_R2.fastq.gz
docker run -v ${PWD}:/mount -u $(id -g ${USER}):$(id -g ${USER}) -it --entrypoint TranslocWrapper.pl htgts tutorial_metadata.txt preprocess/ results/ --threads 2
```

Detect translocation peaks using MACS2
```console
docker run -v ${PWD}:/mount -u $(id -g ${USER}):$(id -g ${USER}) -it --entrypoint tlx2BED-MACS.pl htgts
docker run -v ${PWD}:/mount -u $(id -g ${USER}):$(id -g ${USER}) -it --entrypoint macs2 htgts
```

<a name="docker-inspect">Inspect container</a>
----------------------------------------------------
```console
docker run -v ${PWD}:/mount -u $(id -g ${USER}):$(id -g ${USER}) -it --entrypoint bash htgts
```

<a name="build">Building HTGTS docker image</a>
====================================================

<a name="build-build">Build docker image</a>
----------------------------------------------------
To build Docker image you need to execute
```console
docker build --build-arg http_proxy="http://www.inet.dkfz-heidelberg.de:80" --build-arg https_proxy="http://www.inet.dkfz-heidelberg.de:80" --rm -t sandrejev/htgts:latest .
```

<a name="build-push">Push docker image to Docker HUB</a>
----------------------------------------------------
```console
docker login
docker push sandrejev/htgts:latest
```

<a name="build-convert">Convert cached docker image to singularity (for local testing)</a>
----------------------------------------------------
```console
singularity pull docker-daemon:sandrejev/htgts:latest
```
