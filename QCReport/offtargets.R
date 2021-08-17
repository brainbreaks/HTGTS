library(dplyr)
library(Rsamtools)
library(IRanges)
library(tidyr)

analyze.offtargets = function(sequences, genomes_path, model) {
  #
  # Find primers position in genome
  #
  # TODO: Find multiple alignments

  # bwa index -b 100000000 mm10.fa
  # bowtie-build --threads 30 mm10.fa.gz data/mm10/bowtie1/mm10
  # bwa mem -t 30 -k 5 -C -a -T 10 -B 1 -O 100 -E 100 -L 100 -D 0 -c 70000 -y 60000 mm10.fa baits.fa > tmp/primers.sam


  #
  # @ This should be created in download script
  #
  sequences_fasta = tempfile()
  sequences_bwa_sam = tempfile()
  sequences_bwa_bam = tempfile()
  sequences_bow_sam = tempfile()
  sequences_bow_bam = tempfile()
  model_fasta = file.path(genomes_path, model, paste0(model, ".fa"))
  model_path = file.path(genomes_path, model, model)
  writeLines(paste0(">query", seq_along(sequences), "\n", sequences), con=sequences_fasta)

  system(paste("bwa mem -t 30 -k 5 -C -a -T 10 -B 1 -O 100 -E 100 -L 100 -D 0 -c 70000 -y 60000", model_fasta, sequences_fasta, ">", sequences_bwa_sam))
  system(paste("samtools sort -@ 30", sequences_bwa_sam, ">", sequences_bwa_bam))

  system(paste("bowtie --threads 30 --tryhard --all -v 3 --seedlen 5 -f --sam --no-unal", model_path, sequences_fasta, ">", sequences_bow_sam))
  system(paste("samtools sort -@ 30", sequences_bow_sam, ">", sequences_bow_bam))

  # Analyze genome position and give final output
  primers_alignments_param = Rsamtools::ScanBamParam(tag=c("AS", "NM"), what=c("qname", "rname", "strand", "flag", "pos", "qwidth",  "cigar", "mapq", "qual"))
  primers_alignments = data.frame()
  for(bam_file in c(sequences_bow_bam, sequences_bwa_bam)) {
    primers_alignments.d = lapply(Rsamtools::scanBam(bam_file, param=primers_alignments_param), function(z) {
      d = data.frame(z[names(z)!="tag"])
      d$primer_alignment_mismatches=z$tag$NM
      cbind(d, do.call(rbind, lapply(d$flag, SamSeq::samFlags))) })[[1]] %>%
      dplyr::mutate(primer_sequence_name=qname, primer_alignment_end=pos+qwidth-1) %>%
      dplyr::mutate(primer_alignment_primary=!NOT_PRIMARY_ALIGNMENT) %>%
      dplyr::filter(!READ_UNMAPPED & rname!="chrM") %>%
      dplyr::select(primer_sequence_name, primer_alignment_chrom=rname, primer_alignment_strand=strand, primer_alignment_start=pos, primer_alignment_end, primer_alignment_flag=flag, primer_alignment_cigar=cigar, primer_alignment_len=qwidth, primer_alignment_mismatches, primer_alignment_primary)
    primers_alignments = rbind(primers_alignments, primers_alignments.d)
  }
  primers_alignments_ranges = GenomicRanges::makeGRangesFromDataFrame(primers_alignments %>% dplyr::mutate(chr=primer_alignment_chrom, start=primer_alignment_start, end=primer_alignment_end, primer_alignment_strand))
  primers_alignments$offtarget_sequence = get_seq(model_fasta, primers_alignments_ranges)$sequence
  primers_alignments = primers_alignments %>%
    dplyr::mutate(
      offtarget_sequenceN_total=nchar(gsub("[^N]", "", offtarget_sequence)),
      offtarget_sequence_N_max=max(nchar(stringr::str_split(offtarget_sequence, "[^N]")[[1]]))
    )

  primers_alignments = primers_alignments %>%
    dplyr::filter(primer_alignment_len==nchar(sequences) & primer_alignment_mismatches<=3 & offtarget_sequenceN_total<=1) %>%
    dplyr::distinct(primer_sequence_name, primer_alignment_chrom, primer_alignment_strand, primer_alignment_start, primer_alignment_end, .keep_all=T)

  primers_alignments %>%
    dplyr::filter(primer_alignment_mismatches<=2)

  table(primers_alignments$primer_alignment_mismatches)


  validated_offtargets_ranges = with(validated_offtargets, IRanges::IRanges(start=offtarget_start, end=offtarget_end))
  primers_targets_ranges = with(primers_alignments,  IRanges::IRanges(start=primer_alignment_start, end=primer_alignment_end))
  validated_offtargets_overlapping = validated_offtargets %>%
    dplyr::inner_join(as.data.frame(IRanges::mergeByOverlaps(primers_targets_ranges, validated_offtargets_ranges)), by=c("offtarget_start"="validated_offtargets_ranges.start", "offtarget_end"="validated_offtargets_ranges.end")) %>%
    dplyr::mutate(offtarget_validated=T) %>%
    dplyr::select(bait_name, bait_chrom, offtarget_is_primary, offtarget_chrom, offtarget_strand, offtarget_start, offtarget_end, offtarget_sequence, offtarget_start, offtarget_end, primer_alignment_start=primers_targets_ranges.start, primer_alignment_end=primers_targets_ranges.end, offtarget_validated)

  primers_targets = primers %>%
    dplyr::inner_join(primers_alignments, by="primer_sequence_name") %>%
    dplyr::left_join(validated_offtargets_overlapping, by=c("bait_name", "bait_chrom", "primer_alignment_strand"="offtarget_strand", "primer_alignment_chrom"="offtarget_chrom", "primer_alignment_start", "primer_alignment_end")) %>%
    dplyr::filter(primer_alignment_len==primer_len) %>%
    dplyr::mutate(offtarget_validated=tidyr::replace_na(offtarget_validated, F)) %>%
    # dplyr::select(-primer_sequence_name, -primer_alignment_id) %>%
    dplyr::mutate(primer_alignment_mismatches_rate=primer_alignment_mismatches/primer_len) %>%
    dplyr::filter(primer_alignment_mismatches<=10 & primer_alignment_has_pam & !grepl("_", primer_alignment_chrom)) %>%
    dplyr::filter(primer_name=="sgRNA" & primer_alignment_has_pam)

  primers_targets %>%
    dplyr::mutate(offtarget_validated=ifelse(offtarget_validated, 1, 0)) %>%
    dplyr::select(bait_name, bait_chrom, offtarget_chrom=primer_alignment_chrom, offtarget_strand=primer_alignment_strand, offtarget_start=primer_alignment_start, offtarget_end=primer_alignment_end, offtarget_mismatches=primer_alignment_mismatches, offtarget_primer_sequence=primer_sequence, offtarget_sequence=primer_alignment_sequence) %>%
    readr::write_tsv(file="data/offtargets_predicted.tsv")
}
