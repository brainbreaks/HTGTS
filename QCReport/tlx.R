source("utils.R")
library(readr)
library(dplyr)
library(GenomicRanges)
library(data.table)

tlx_cols = cols(
  Qname=readr::col_character(), JuncID=readr::col_character(), Rname=readr::col_character(), Junction=readr::col_double(),
  Strand=readr::col_character(), Rstart=readr::col_double(), Rend=readr::col_double(),
  B_Rname=readr::col_character(), B_Rstart=readr::col_double(), B_Rend=readr::col_double(), B_Strand=readr::col_double(),
  B_Qstart=readr::col_double(), B_Qend=readr::col_double(), Qstart=readr::col_double(), Qend=readr::col_double(), Qlen=readr::col_double(),
  B_Cigar=readr::col_character(), Cigar=readr::col_character(), Seq=readr::col_character(), J_Seq=readr::col_character(), Barcode=readr::col_logical(),
  unaligned=readr::col_double(), baitonly=readr::col_double(), uncut=readr::col_double(), misprimed=readr::col_double(), freqcut=readr::col_double(),
  largegap=readr::col_double(), mapqual=readr::col_double(), breaksite=readr::col_double(), sequential=readr::col_double(), repeatseq=readr::col_double(), duplicate=readr::col_double()
)

tlx_blank = function() {
  as.data.frame(t(sapply(names(tlx_cols$cols), function(z) NA))) %>%
    dplyr::mutate(tlx_sample=NA, tlx_path=NA, tlx_group=NA, tlx_control=NA) %>%
    dplyr::slice(0)
}

tlx_read = function(path, sample, group="", control=F) {
  readr::read_tsv(path, comment="#", skip=16, col_names=names(tlx_cols$cols), col_types=tlx_cols) %>%
      dplyr::mutate(tlx_sample=sample, tlx_path=path, tlx_group=group, tlx_control=control)
}

tlx_read_many = function(samples_df) {
  tlx_df.all = data.frame()
  for(f in 1:nrow(samples_df)) {
    log("Reading tlx file ", f, ":",  samples_df$path[f])
    tlx_df.f = tlx_read(samples_df$path[f], sample=samples_df$sample[f], control=samples_df$control[f], group=samples_df$group[f])
    tlx_df.all = dplyr::bind_rows(tlx_df.all, tlx_df.f)
  }

  tlx_df.all
}

tlx_remove_rand_chromosomes = function(tlx_df) {
  tlx_df %>%
    dplyr::filter(Rname %in% paste0("chr", c(1:40, "X", "Y")))
}


tlx_identify_baits = function(tlx_df, breaksite_size=19) {
  if(is.null(tlx_df) || nrow(tlx_df)==0) {
    return(data.frame(bait_sample=NA, bait_chrom=NA, bait_strand=NA, bait_start=NA, bait_end=NA) %>% dplyr::slice(0))
  }

  baits_df = tlx_df %>%
    dplyr::group_by(tlx_group, tlx_sample, B_Rname, B_Strand) %>%
    dplyr::do((function(z){
      misprimed_max = max(z$misprimed-z$uncut)
      if(z$B_Strand[1]<0) bait_start = unique(z$B_Rstart + z$misprimed - misprimed_max-2)
      else bait_start = unique(z$B_Rend - z$misprimed + misprimed_max)
      data.frame(bait_start=bait_start, bait_end=bait_start+breaksite_size - 1)
    })(.)) %>%
    dplyr::mutate(bait_strand=ifelse(B_Strand<0, "-", "+")) %>%
    dplyr::ungroup() %>%
    dplyr::select(bait_group=tlx_group, bait_sample=tlx_sample, bait_chrom=B_Rname, bait_strand, bait_start, bait_end)

  baits_df
}

tlx_mark_bait_chromosome = function(tlx_df) {
  tlx_df %>%
    dplyr::select(-dplyr::matches("tlx_is_bait_chromosome")) %>%
    dplyr::mutate(tlx_is_bait_chromosome=B_Rname==Rname)
}

tlx_mark_bait_junctions = function(tlx_df, bait_region) {
  tlx_df %>%
    dplyr::select(-dplyr::matches("tlx_is_bait_junction")) %>%
    dplyr::mutate(tlx_is_bait_junction=B_Rname==Rname & (abs(B_Rstart-Rstart)<=bait_region/2 | abs(Rend-B_Rend)<=bait_region/2))
}

tlx_mark_offtargets = function(tlx_df, offtarget2bait_df) {
  # @todo: Change 100 to something meaningful?
  tlx_df$tlx_id = 1:nrow(tlx_df)
  tlx_bait_ranges = GenomicRanges::makeGRangesFromDataFrame(tlx_df %>% dplyr::mutate(seqnames=B_Rname, start=B_Rstart-50000, end=B_Rend+50000), keep.extra.columns=T, ignore.strand=T)
  tlx_junc_ranges = GenomicRanges::makeGRangesFromDataFrame(tlx_df %>% dplyr::mutate(seqnames=Rname, start=Rstart-50000, end=Rend+50000), keep.extra.columns=T, ignore.strand=T)
  offtarget2bait_bait_ranges = GenomicRanges::makeGRangesFromDataFrame(offtarget2bait_df %>% dplyr::mutate(seqnames=bait_chrom, start=bait_start, end=bait_end), ignore.strand=T)
  offtarget2bait_offt_ranges = GenomicRanges::makeGRangesFromDataFrame(offtarget2bait_df %>% dplyr::mutate(seqnames=offtarget_chrom, start=offtarget_start, end=offtarget_end), ignore.strand=T)

  tlx_offtarget_ids = as.data.frame(IRanges::findOverlaps(tlx_bait_ranges, offtarget2bait_bait_ranges)) %>%
    dplyr::rename(tlx_id="queryHits", o2b_id="subjectHits") %>%
    dplyr::inner_join(as.data.frame(IRanges::findOverlaps(tlx_junc_ranges, offtarget2bait_offt_ranges)), by=c(tlx_id="queryHits", o2b_id="subjectHits")) %>%
    dplyr::distinct(tlx_id) %>%
    .$tlx_id

  tlx_df$tlx_is_offtarget = tlx_df$tlx_id %in% tlx_offtarget_ids

  tlx_df
}

tlx_mark_repeats = function(tlx_df, repeatmasker_df) {
  # @todo: make group_by faster using data.table
  repeatmasker_ranges = GenomicRanges::makeGRangesFromDataFrame(repeatmasker_df %>% dplyr::mutate(seqnames=repeatmasker_chrom, start=repeatmasker_start, end=repeatmasker_end), keep.extra.columns=T)
  tlx_df = tlx_df %>%
    dplyr::mutate(tlx_id=1:n()) %>%
    dplyr::select(-dplyr::matches("tlx_repeatmasker_")) %>%
    dplyr::mutate(seqnames=Rname, start=Rstart, end=Rend)
  tlx_ranges = GenomicRanges::makeGRangesFromDataFrame(tlx_df, keep.extra.columns=T, ignore.strand=T)
  r1 = as.data.frame(IRanges::findOverlaps(tlx_ranges, repeatmasker_ranges)) %>%
    dplyr::inner_join(repeatmasker_df, by=c("subjectHits"="repeatmasker_id"))
  data.table::setDT(r1)[,.(tlx_repeatmasker_class=paste0(unique(repeatmasker_class),collapse=", ")), by = .(queryHits)] %>%
    dplyr::right_join(tlx_df, by=c("queryHits"="tlx_id")) %>%
    dplyr::select(-queryHits) %>%
    data.frame()
}

tlx_macs2 = function(tlx_df, qvalue=0.01, pileup=1, extsize=2000, slocal=50000, llocal=10000000, exclude_bait_region=F, exclude_repeats=F, exclude_offtargets=F, exttype=c("along", "symmetrical")) {
  if(exclude_bait_region && !("tlx_is_bait_junction" %in% colnames(tlx_df))) {
    stop("tlx_is_bait_junction is not found in tlx data frame")
  }

  if(exclude_offtargets) {
    if(!("tlx_is_offtarget" %in% colnames(tlx_df))) {
      stop("tlx_is_offtarget is not found in tlx data frame")
    }
    tlx_df = tlx_df %>% dplyr::filter(!tlx_is_offtarget)
  }
  if(exclude_repeats) {
    if(!("tlx_repeatmasker_class" %in% colnames(tlx_df))) {
      stop("tlx_repeatmasker_class is not found in tlx data frame")
    }
    tlx_df = tlx_df %>% dplyr::filter(is.na(tlx_repeatmasker_class))
  }

  tlx_df = tlx_df %>%
    dplyr::filter(!exclude_bait_region | !tlx_is_bait_junction) %>%
    dplyr::mutate(bed_strand=ifelse(Strand=="1", "-", "+"), bed_start=Junction, bed_end=Junction+1)

  # @TODO: I think macs does this internally
  # if(exttype=="along") {
  #     tlx_df = tlx_df %>% dplyr::mutate(bed_start=ifelse(Strand=="-1", Junction-junsize, Junction-1), bed_end=ifelse(Strand=="-1", Junction, Junction+junsize-1))
  # } else {
  #     tlx_df = tlx_df %>% dplyr::mutate(bed_start=Junction-ceiling(junsize/2), bed_end=Junction+ceiling(junsize/2))
  # }

  macs_df.all = data.frame()
  for(gr in unique(tlx_df$tlx_group)) {
    tlx_df.gr = tlx_df %>% dplyr::filter(tlx_group==gr)

    f_input_bed = tempfile()
    f_control_bed = tempfile()

    tlx_df.gr %>%
      dplyr::filter(!tlx_control) %>%
      dplyr::select(Rname, bed_start, bed_end, Qname, 0, bed_strand) %>%
      readr::write_tsv(file=f_input_bed, na="", col_names=F)

    if(any(tlx_df.gr$tlx_control)) {
      tlx_df.gr %>%
        dplyr::filter(tlx_control) %>%
        dplyr::select(Rname, bed_start, bed_end, Qname, 0, bed_strand) %>%
        readr::write_tsv(file=f_control_bed, na="", col_names=F)

      print("Running MACS with control")
      macs_df = macs2(name=basename(f_input_bed), sample=f_input_bed, control=f_control_bed, extsize=extsize, qvalue=qvalue, slocal=slocal, llocal=llocal, output_dir=dirname(f_input_bed)) %>%
        dplyr::mutate(macs_group=gr)
    } else {
      print("Running MACS without control")
      macs_df = macs2(name=basename(f_input_bed), sample=f_input_bed, extsize=extsize, qvalue=qvalue, slocal=slocal, llocal=llocal, output_dir=dirname(f_input_bed)) %>%
        dplyr::mutate(macs_group=gr)
    }

    macs_df.all = rbind(macs_df.all, macs_df %>% dplyr::mutate(macs_group=gr))
  }

  macs_df.all = macs_df.all %>% dplyr::filter(macs_pileup>=pileup)

  macs_df.all
}