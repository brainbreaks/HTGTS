source("utils.R")
library(readr)
library(dplyr)
library(GenomicRanges)

tlx_cols = cols(
  Qname=readr::col_character(), JuncID=readr::col_character(), Rname=readr::col_character(), Junction=readr::col_double(),
  Strand=readr::col_character(), Rstart=readr::col_double(), Rend=readr::col_double(),
  B_Rname=readr::col_character(), B_Rstart=readr::col_double(), B_Rend=readr::col_double(), B_Strand=readr::col_double(),
  B_Qstart=readr::col_double(), B_Qend=readr::col_double(), Qstart=readr::col_double(), Qend=readr::col_double(), Qlen=readr::col_double(),
  B_Cigar=readr::col_character(), Cigar=readr::col_character(), Seq=readr::col_character(), J_Seq=readr::col_character(), Barcode=readr::col_logical(),
  unaligned=readr::col_double(), baitonly=readr::col_double(), uncut=readr::col_double(), misprimed=readr::col_double(), freqcut=readr::col_double(),
  largegap=readr::col_double(), mapqual=readr::col_double(), breaksite=readr::col_double(), sequential=readr::col_double(), repeatseq=readr::col_double(), duplicate=readr::col_double()
)

tlx_read = function(path, sample) {
  readr::read_tsv(path, comment="#", skip=16, col_names=names(tlx_cols$cols), col_types=tlx_cols) %>%
      dplyr::mutate(tlx_sample=sample)
}

tlx_identify_baits = function(tlx_df, size=20, junsize=200) {
  # @todo: Check this arithmetic with "-" strand break
  tlx_df %>%
    dplyr::group_by(B_Rname, B_Rstart, B_Rend, B_Strand) %>%
    dplyr::mutate(n=n(), B_RstartGroup=floor(B_Rstart/junsize)*junsize) %>%
    dplyr::group_by(tlx_sample, B_RstartGroup) %>%
    dplyr::arrange(dplyr::desc(n)) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(bait_strand=ifelse(B_Strand=="1", "+", "-"), bait_start=ifelse(bait_strand=="+", B_Rend-2, B_Rstart+2), bait_end=bait_start+size-1) %>%
    dplyr::select(bait_sample=tlx_sample, bait_chrom=B_Rname, bait_strand, bait_start, bait_end)
}

tlx_mark_bait_junctions = function(tlx_df, bait_region) {
  tlx_df %>%
    dplyr::select(-dplyr::matches("tlx_is_bait_junction")) %>%
    dplyr::mutate(tlx_is_bait_junction=B_Rname==Rname & (abs(B_Rstart-Rstart)<=bait_region/2 | abs(Rend-B_Rend)<=bait_region/2))
}

tlx_mark_offtargets = function(tlx_df, offtarget2bait_df) {
  # @todo: Change 100 to something meaningful?
  tlx_df$tlx_id = 1:nrow(tlx_df)
  tlx_bait_ranges = GenomicRanges::makeGRangesFromDataFrame(tlx_df %>% dplyr::mutate(seqnames=B_Rname, start=B_Rstart-100, end=B_Rend+100), keep.extra.columns=T, ignore.strand=T)
  tlx_junc_ranges = GenomicRanges::makeGRangesFromDataFrame(tlx_df %>% dplyr::mutate(seqnames=Rname, start=Rstart-100, end=Rend+100), keep.extra.columns=T, ignore.strand=T)
  offtarget2bait_bait_ranges = GenomicRanges::makeGRangesFromDataFrame(offtarget2bait_df %>% dplyr::mutate(seqnames=bait_chrom, start=bait_start, end=bait_end), keep.extra.columns=T, ignore.strand=T)
  offtarget2bait_offt_ranges = GenomicRanges::makeGRangesFromDataFrame(offtarget2bait_df %>% dplyr::mutate(seqnames=offtarget_chrom, start=offtarget_start, end=offtarget_end), keep.extra.columns=T, ignore.strand=T)

  tlx_offtarget_ids = as.data.frame(IRanges::findOverlaps(tlx_bait_ranges, offtarget2bait_bait_ranges)) %>%
    dplyr::rename(tlx_id="queryHits", o2b_id="subjectHits") %>%
    dplyr::inner_join(as.data.frame(IRanges::findOverlaps(tlx_junc_ranges, offtarget2bait_offt_ranges)), by=c(tlx_id="queryHits", o2b_id="subjectHits")) %>%
    .$tlx_id

  tlx_df$tlx_is_offtarget = tlx_df$tlx_id %in% tlx_offtarget_ids

  tlx_df
}

tlx_mark_repeats = function(tlx_df, repeatmasker_df) {
  repeatmasker_ranges = GenomicRanges::makeGRangesFromDataFrame(repeatmasker_df %>% dplyr::mutate(seqnames=repeatmasker_chrom, start=repeatmasker_start, end=repeatmasker_end), keep.extra.columns=T)
  tlx_df = tlx_df %>%
    dplyr::mutate(tlx_id=1:n()) %>%
    dplyr::select(-dplyr::matches("tlx_repeatmasker_")) %>%
    dplyr::mutate(seqnames=Rname, start=Rstart, end=Rend)
  tlx_ranges = GenomicRanges::makeGRangesFromDataFrame(tlx_df, keep.extra.columns=T, ignore.strand=T)
  as.data.frame(IRanges::findOverlaps(tlx_ranges, repeatmasker_ranges)) %>%
    dplyr::inner_join(repeatmasker_df, by=c("subjectHits"="repeatmasker_id")) %>%
    dplyr::group_by(queryHits) %>%
    dplyr::summarise(tlx_repeatmasker_name=paste(unique(repeatmasker_name), collapse=", "), tlx_repeatmasker_class=paste(unique(repeatmasker_class), collapse=", "), tlx_repeatmasker_family=paste(unique(repeatmasker_family), collapse=", "), tlx_repeatmasker_length=max(abs(repeatmasker_end-repeatmasker_start))) %>%
    dplyr::right_join(tlx_df, by=c("queryHits"="tlx_id")) %>%
    dplyr::select(-queryHits)
}

tlx_macs2 = function(tlx_df, name=NULL, qvalue=0.01, extsize=200, slocal=1000, llocal=10000000, exclude_bait_region=F, exclude_repeats=F, exclude_offtargets=F) {
  if(exclude_offtargets && !("tlx_is_offtarget" %in% colnames(tlx_df))) {
    stop("tlx_is_offtarget is not found in tlx data frame")
  }
  if(exclude_repeats && !("tlx_repeatmasker_class" %in% colnames(tlx_df))) {
    stop("tlx_repeatmasker_class is not found in tlx data frame")
  }
  if(exclude_bait_region && !("tlx_is_bait_junction" %in% colnames(tlx_df))) {
    stop("tlx_is_bait_junction is not found in tlx data frame")
  }
  if(is.null(name) && !("tlx_sample" %in% colnames(tlx_df))) {
    stop("tlx_sample is not found in tlx data frame")
  }

  if(is.null(name)) {
    name = tlx_df$tlx_sample[1]
  }

  f_bed = tempfile()
  print(paste0("Convert to BED (", f_bed, ")"))
  tlx_df %>%
    dplyr::filter(!exclude_offtargets | !tlx_is_offtarget) %>%
    dplyr::filter(!exclude_bait_region | !tlx_is_bait_junction) %>%
    dplyr::filter(!exclude_repeats | is.na(tlx_repeatmasker_class)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(bed_start=ifelse(Strand=="1", Junction-junsize, Junction-1), bed_end=ifelse(Strand=="1", Junction, Junction+junsize-1), bed_strand=ifelse(Strand=="1", "-", "+")) %>%
    dplyr::select(Rname, bed_start, bed_end, Qname, 0, bed_strand) %>%
    readr::write_tsv(file=f_bed, na="", col_names=F)
  print("Running MACS")
  macs_df = macs2(name=basename(f_bed), sample=f_bed, extsize=extsize, qvalue=qvalue, slocal=slocal, llocal=llocal, output_dir=dirname(f_bed)) %>%
    dplyr::mutate(macs_sample=name)
}