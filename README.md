singularity pull docker-daemon:htgts:latest
singularity exec -B `pwd` htgts_latest.sif download mm9
singularity exec -B `pwd` htgts_latest.sif download hg19
singularity exec -B `pwd` htgts_latest.sif TranslocPreprocess.pl tutorial_metadata.txt preprocess --read1 pooled_R1.fastq.gz --read2 pooled_R2.fastq.gz
singularity exec -B `pwd` htgts_latest.sif TranslocWrapper.pl tutorial_metadata.txt preprocess/ results/ --threads 2
