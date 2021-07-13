library(stringr)
library(dplyr)
library(readr)
library(bedr)

macs_cols = cols(
  macs_chrom=col_character(), macs_start=col_double(), macs_end=col_double(), macs_length=col_character(), macs_summit_abs=col_double(),
  macs_pileup=col_double(), macs_pvalue=col_double(), macs_fc=col_double(), macs_qvalue=col_double(), macs_name=col_character(), macs_comment=col_character()
)

bed_read = function(path) {
  bed = rtracklayer::import.bed(path)
  GenomicRanges::start(bed) = GenomicRanges::start(bed)-1
  bed
}

# one-based index
get_seq = function(fasta, ranges) {
  bed_df = data.frame(
    chr=as.character(GenomicRanges::seqnames(ranges)),
    start=as.numeric(GenomicRanges::start(ranges))-1,
    end=as.numeric(GenomicRanges::end(ranges)),
    strand=as.character(GenomicRanges::strand(ranges))) %>%
    dplyr::mutate(name="", score=0) %>%
    dplyr::select(chr, start, end, name, score, strand)
  bed_df.order = with(bed_df, order(chr, start, end, strand))
  res = bedr::get.fasta(bed_df[bed_df.order,], fasta=fasta, strand=T, check.zero.based=F, check.chr=F, check.valid=F, check.sort=T, check.merge=F)
  ranges$sequence = res$sequence[match(1:nrow(bed_df), bed_df.order)]

  ranges
}

repeatmasker_cols = cols(
  repeatmasker_bin=readr::col_double(),
  repeatmasker_score=readr::col_double(),
  repeatmasker_mismatches_per_kb=readr::col_double(),
  repeatmasker_deletions_per_kb=readr::col_double(),
  repeatmasker_insertions_per_kb=readr::col_double(),
  repeatmasker_chrom=readr::col_character(),
  repeatmasker_start=readr::col_double(),
  repeatmasker_end=readr::col_double(),
  repeatmasker_genoLeft=readr::col_double(),
  repeatmasker_strand=readr::col_character(),
  repeatmasker_name=readr::col_character(),
  repeatmasker_class=readr::col_character(),
  repeatmasker_family=readr::col_character(),
  repeatmasker_repStart=readr::col_character(),
  repeatmasker_repEnd=readr::col_character(),
  repeatmasker_repLeft=readr::col_character(),
  repeatmasker_id=readr::col_character()
)


macs2 = function(name, sample, control=NULL, qvalue=0.01, extsize=200, slocal=1000, llocal=10000000, output_dir="data/macs2") {
  bed_sample = paste("-t", sample)
  bed_control = ifelse(is.null(control), "", paste("-c", control))

  cmd = stringr::str_glue("macs2 callpeak {bed_sample} {bed_control} --seed 123 -f BED -g hs --keep-dup all -n {name} --outdir {output_dir} --nomodel --slocal {slocal} --extsize {extsize} -q {qvalue} --llocal {llocal} --bdg --trackline", bed_sample=bed_sample, bed_control=bed_control, name=name, output_dir=output_dir, extsize=extsize, qvalue=qvalue, llocal=sprintf("%0.0f", llocal), slocal=sprintf("%0.0f", slocal))
  print(cmd)
  system(cmd)

  readr::read_tsv(paste0(output_dir, "/", name, "_peaks.xls"), comment="#", col_names=names(macs_cols$cols), col_types=macs_cols) %>%
    dplyr::slice(-1) %>%
    dplyr::select(-macs_comment)
}

join_offtarget2bait = function(offtargets_df, baits_df, genome_path) {
  baits_ranges = GenomicRanges::makeGRangesFromDataFrame(baits_df %>% dplyr::mutate(seqnames=bait_chrom, start=bait_start, end=bait_end, strand=bait_strand))
  offtargets_ranges = GenomicRanges::makeGRangesFromDataFrame(offtargets_df %>% dplyr::mutate(seqnames=offtarget_chrom, start=offtarget_start, end=offtarget_end, strand=offtarget_strand))

  # Combine baits and offtarget ranges so that sequences can be retrieve in a single call to bedtools (performance optimization)
  ranges_seq = get_seq(fasta=genome_path, ranges=BiocGenerics::append(baits_ranges, offtargets_ranges))
  baits_df$bait_sequence = ranges_seq$sequence[1:length(baits_ranges)]
  offtargets_df$offtarget_sequence = ranges_seq$sequence[-(1:length(baits_ranges))]

  offtarget2bait_df = offtargets_df %>%
    tidyr::crossing(baits_df) %>%
    dplyr::mutate(bait2offtarget_alignment=Biostrings::pairwiseAlignment(offtarget_sequence, bait_sequence, type="global", scoreOnly=T)) %>%
    dplyr::arrange(bait2offtarget_alignment) %>%
    dplyr::distinct(offtarget_chrom, offtarget_start, offtarget_end, offtarget_strand, .keep_all=T)

  offtarget2bait_df
}

offtargets_read = function(path) {
  offtargets_ranges = bed_read(path)
  as.data.frame(offtargets_ranges) %>%
    dplyr::select(offtarget_chrom=seqnames, offtarget_start=start, offtarget_end=end, offtarget_strand=strand)
}