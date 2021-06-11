# Installation

## Download code

```
cd ~
git clone https://github.com/robinmeyers/transloc_pipeline.git
```

Add directories to $PATH in ~/.profile or ~/.bash_profile

```
echo 'export PATH=~/transloc_pipeline/bin:~/transloc_pipeline/R:$PATH' >> ~/.bash_profile
```

## Reference Genomes

```
mkdir -p ~/genomes/bowtie2_indexes
cd ~/genomes/bowtie2_indexes
wget ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/hg19.zip
unzip hg19.zip
mkdir ~/genomes/hg19
bowtie2-inspect hg19 > ~/genomes/hg19/hg19.fa 
```

Add environment variables to in ~/.profile or ~/.bash_profile

```
echo 'export BOWTIE2_INDEXES=~/genomes/bowtie2_indexes' >> ~/.bash_profile
echo 'export GENOME_DB=~/genomes' >> ~/.bash_profile
```



## Dependencies

### Bowtie2

Install [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

If using homebrew:
```
$ brew tap homebrew/science
$ brew install bowtie2
```

### Perl >= 5.16

Install [Bioperl](http://www.bioperl.org/wiki/Installing_BioPerl_on_Unix)

```
$ cpan
```
Find the most recent version

```
cpan> d /bioperl/
Distribution    CJFIELDS/BioPerl-1.6.901.tar.gz
Distribution    CJFIELDS/BioPerl-1.6.923.tar.gz
Distribution    CJFIELDS/BioPerl-1.6.924.tar.gz
```

Now install:

```
cpan> install CJFIELDS/BioPerl-1.6.924.tar.gz
```
or
```
cpan> force install CJFIELDS/BioPerl-1.6.924.tar.gz
```

Other modules:

- 


# Running the Pipeline


## Pre-processing libraries

The Alt Lab primarily uses in-line barcodes (sequenced by the MiSeq at the head of the forward read) and then deconvolutes pooled libraries using this program. This script also calls an external tool to trim Illumina adapter sequence. If using Illumina multi-plex barcoding strategy, this script will not be useful except for trimming adapters, which is still recommended.

### Starting from pooled library fastq files
Deconvolutes and trims adapters
```
$ cd ~/transloc_pipeline/data
$ TranslocPreprocess.pl metadata.txt preprocess/ --read1 pooled_R1.fq.gz --read2 pooled_R2.fq.gz
```

### Starting from deconvoluted fastq files
Just trims adapters
```
$ cd ~/transloc_pipeline/data
$ TranslocPreprocess.pl metadata.txt preprocess/ --indir ./
```

## Running the pipeline

```
$ TranslocWrapper.pl metadata.txt preprocess/ results/ --threads 2
```

# Filtering

The TranslocPipeline main output is the *.tlx file. One master tlx file will be generated per library. It can be thought of as similar to a sam file for NGS libraries. It contains output for every read in the library, and will serve as the starting point for all down-stream analyses. However, for most analyses, it will need to be filtered into only the required reads for the specific analysis. There is a default filtering that the pipeline will generate automatically. For any other filtering regime, the user will have to re-filter the master tlx file.

## Available filters

**- unaligned:** No OCS alignments.

**- baitonly:** Bait alignment is either the only alignent in the OCS or  only followed by adapter alignment.

**- uncut:** Bait alignment runs greater than some number of bases past the cutsite.

**- misprimed:** Bait alignment runs fewer than some number of bases past the primer.

**- freqcut:** Restriction enzyme site within some number of bases from the junction. (Fairly depricated.)

**- largegap:** More than some number of bases between the bait and prey alignments.

**- mapqual:** OCS had a competing prey junction.

**- breaksite:** Prey alignment maps into non-endogenous breaksite cassette. (Fairly depricated.)

**- sequential:** Junction occurs downstream on read from first bait-prey junction.


## Re-filtering a library

```
$ TranslocFilter.pl results/RF204_Alt055/RF204_Alt055.tlx results/RF204_Alt055_refiltered.tlx --filters ""
```

### Ex. 1 Keep duplicate junctions

```
$ TranslocFilter.pl results/RF204_Alt055/RF204_Alt055.tlx results/RF204_Alt055_refiltered.tlx --filters ""
```

### Ex. 2 Keep all un-translocated reads

```
$ TranslocFilter.pl results/RF204_Alt055/RF204_Alt055.tlx results/RF204_Alt055_refiltered.tlx --filters ""
```

### Using a config file

# Hotspot Detection

Two methods of detecting hotspots.

## Using MACS2 to do translocaiton peak detection

Must have MACS2 installed.

```
$ tlx2BED-MACS.pl
$ macs2
```

## Using scan statistics script to call peaks
```
$ TranslocHotspots.R
```

