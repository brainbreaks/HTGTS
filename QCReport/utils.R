library(stringr)
library(dplyr)
library(readr)
library(bedr)
library(digest)

log = function(..., collapse=NULL) {
  msg = paste0(">>> ", paste0(..., collapse=collapse))
  message(msg)
  # cat(file=stderr(), paste0(msg, "\n"))
}

blank_tibble = function(cols) {
  stopifnot(class(tlx_cols)=="col_spec")

  x = tibble::tibble()
  for(n in names(cols$cols)) {
    n_class = class(cols$cols[[n]])
    n_type = match.fun(gsub("collector_", "", n_class[grepl("collector_", n_class)]))
    x[[n]] = n_type()
  }

  x
}

log_input = function(input) {
  vals = reactiveValuesToList(input)
  vals_filter = !grepl("shinyActionButtonValue", sapply(vals, function(z) {paste(class(z), collapse="|")}))
  vals = vals[vals_filter]
  vals_pretty = sapply(vals, function(z) {
    if(is.data.frame(z) & "name" %in% colnames(z)) z = paste(z$name, collapse=", ")
    if(is.null(z)) z = "NULL"
    as.character(z)
  })

  log(names(vals_pretty), "=", vals_pretty, collapse="\n")
}

macs_cols = cols(
  macs_chrom=col_character(), macs_start=col_double(), macs_end=col_double(), macs_length=col_character(), macs_summit_abs=col_double(),
  macs_pileup=col_double(), macs_pvalue=col_double(), macs_fc=col_double(), macs_qvalue=col_double(), macs_name=col_character(), macs_comment=col_character()
)

macs_blank = function() {
  blank_tibble(macs_cols) %>% dplyr::mutate(macs_group=NA_character_)
}

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


get_seq2 = function(fasta, chrom, start, end, strand) {
  bed_df = data.frame(
    chr=chrom,
    start=start-1,
    end=end,
    strand=strand) %>%
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
  repeatmasker_chrom=readr::col_factor(),
  repeatmasker_start=readr::col_double(),
  repeatmasker_end=readr::col_double(),
  repeatmasker_genoLeft=readr::col_double(),
  repeatmasker_strand=readr::col_factor(),
  repeatmasker_name=readr::col_factor(),
  repeatmasker_class=readr::col_factor(),
  repeatmasker_family=readr::col_factor(),
  repeatmasker_repStart=readr::col_character(),
  repeatmasker_repEnd=readr::col_character(),
  repeatmasker_repLeft=readr::col_character(),
  repeatmasker_id=readr::col_factor()
)

repeatmasker_read = function(path, columns=c("repeatmasker_chrom", "repeatmasker_start", "repeatmasker_end", "repeatmasker_strand", "repeatmasker_name", "repeatmasker_class", "repeatmasker_family")) {
  cache_path = paste0("tmp/", digest::digest(paste(path, paste(columns, collapse="")), algo="md5"), "_repeatmasker_df.rda")
  if(file.exists(cache_path)) {
    load(cache_path)
    return(repeatmasker_df %>% dplyr::mutate(repeatmasker_id=1:n()))
  }
  else {
    if(!dir.exists("tmp")) {
      dir.create("tmp", recursive=T)
    }

    repeatmasker_df = readr::read_tsv(path, col_names=names(repeatmasker_cols$cols), col_types=repeatmasker_cols, skip=1) %>%
      dplyr::filter(repeatmasker_chrom %in% paste0("chr", c(1:99, "X", "Y"))) %>%
      data.frame()
    repeatmasker_df = repeatmasker_df[columns]


    save(repeatmasker_df, file=cache_path)
    return(repeatmasker_df %>% dplyr::mutate(repeatmasker_id=1:n()))
  }
}

macs2 = function(name, sample, effective_size, control=NULL, maxgap=NULL, qvalue=0.01, extsize=2000, slocal=50000, llocal=10000000, output_dir="data/macs2") {
  bed_sample = paste("-t", sample)
  bed_control = ifelse(is.null(control), "", paste("-c", control))
  maxgap = ifelse(is.null(maxgap), "", paste("--max-gap", sprintf("%0.0f", maxgap)))

  cmd = stringr::str_glue("macs2 callpeak {bed_sample} {bed_control} --seed 123 {maxgap} -f BED -g {effsize} --keep-dup all -n {name} --outdir {output_dir} --nomodel --slocal {slocal} --extsize {extsize} -q {qvalue} --llocal {llocal} --bdg --trackline", bed_sample=bed_sample, bed_control=bed_control, name=name, output_dir=output_dir, extsize=extsize, qvalue=qvalue, maxgap=maxgap, llocal=sprintf("%0.0f", llocal), slocal=sprintf("%0.0f", slocal), effsize=effective_size)
  log(cmd)
  output = system(paste(cmd, " 2>&1"), intern = T)
  output = paste0(output, collapse="\n")
  log(output)

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
    dplyr::distinct(bait_sample, offtarget_chrom, offtarget_start, offtarget_end, offtarget_strand, .keep_all=T) %>%
    data.frame()

  offtarget2bait_df
}

offtargets_read = function(path) {
  offtargets_ranges = bed_read(path)
  as.data.frame(offtargets_ranges) %>%
    dplyr::select(offtarget_chrom=seqnames, offtarget_start=start, offtarget_end=end, offtarget_strand=strand)
}