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
  blank_tibble(tlx_cols) %>%
    dplyr::mutate(tlx_sample=NA_character_, tlx_path=NA_character_, tlx_group=NA_character_, tlx_control=NA) %>%
    dplyr::mutate(tlx_is_bait_chromosome=NA, tlx_is_bait_junction=NA, tlx_is_offtarget=NA) %>%
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


tlx_write_wig = function(tlx_df, file, extsize) {
  tld_df_mod = tlx_df %>%
    dplyr::mutate(Rstart=ifelse(Strand=="-1", Junction-extsize, Junction), Rend=ifelse(Strand=="1", Junction+extsize, Junction))
  tlx_coverage(tld_df_mod, group="none") %>%
    dplyr::select(tlxcov_chrom, tlxcov_start, tlxcov_end, tlxcov_pileup) %>%
    readr::write_tsv(file=file, col_names=F)
}


tlx_write_bed = function(tlx_df, file) {
  tlx_df %>%
    dplyr::mutate(strand=ifelse(Strand=="-1", "-", "+"), start=ifelse(Strand=="-1", Junction-1, Junction), end=ifelse(Strand=="-1", Junction, Junction+1)) %>%
    dplyr::select(Rname, start, end, Qname, mapqual, strand) %>%
    readr::write_tsv(file=file, col_names=F)
}

tlx_coverage = function(tlx_df, group=c("none", "group", "sample", "path")) {
  tlx_coverage_ = function(x) {
    x_ranges = GenomicRanges::makeGRangesFromDataFrame(x %>% dplyr::mutate(seqnames=Rname, start=Rstart, end=Rend), ignore.strand=T)
    cov_ranges = as(GenomicRanges::coverage(x_ranges), "GRanges")
    as.data.frame(cov_ranges) %>%
      dplyr::rename(tlxcov_chrom="seqnames", tlxcov_start="start", tlxcov_end="end", tlxcov_pileup="score") %>%
      dplyr::select(matches("tlxcov_"))
  }

  if(group=="none") {
    return(tlx_coverage_(tlx_df))
  }
  if(group=="group") {
    return(tlx_df %>% group_by(tlx_group) %>% do(tlx_coverage_(.)) %>% dplyr::ungroup())
  }

  if(group=="sample") {
    return(tlx_df %>% group_by(tlx_group, tlx_sample) %>% do(tlx_coverage_(.)) %>% dplyr::ungroup())
  }

  if(group=="path") {
    return(tlx_df %>% group_by(tlx_group, tlx_sample, tlx_path) %>% do(tlx_coverage_(.)) %>% dplyr::ungroup())
  }

  stop("Unknown group")
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


tlx_test_hits = function(tlx_df, hits_ranges) {
  hits_df = as.data.frame(hits_ranges) %>% dplyr::mutate(compare_chrom=seqnames, compare_start=start, compare_end=end)
  hits_ranges = GenomicRanges::makeGRangesFromDataFrame(hits_df, keep.extra.columns=T)
  tlx_ranges  = GenomicRanges::makeGRangesFromDataFrame(tlx_df %>% dplyr::mutate(seqnames=Rname, start=Rstart, end=Rend), ignore.strand=T, keep.extra.columns=T)
  tlxsum_df = tlx_df %>%
    dplyr::group_by(tlx_group, tlx_sample, .drop=F) %>%
    dplyr::summarize(total_n=sum(!tlx_is_bait_junction)) %>%
    dplyr::ungroup()

  as.data.frame(IRanges::mergeByOverlaps(hits_ranges, tlx_ranges)) %>%
    dplyr::inner_join(tlxsum_df, by=c("tlx_group", "tlx_sample")) %>%
    dplyr::group_by(compare_chrom, compare_start, compare_end, tlx_group, tlx_sample, total_n, .drop=F) %>%
    dplyr::summarize(n=n(), total_n=tidyr::replace_na(total_n[1], 0)) %>%
    dplyr::group_by(compare_chrom, compare_start, compare_end) %>%
    dplyr::do((function(z){
      zz<<-z
      do.call(rbind, apply(combn(unique(z$tlx_group), 2), 2, function(cgr) {
        zz.cgr<<-cgr

        paired=F
        if(paired) {
          z.gr = expand.grid(gr1=which(z$tlx_group==cgr[1]), gr2=which(z$tlx_group==cgr[2]))
        } else {
          z.gr = data.frame(gr1=which(z$tlx_group==cgr[1]), gr2=which(z$tlx_group==cgr[2]))
        }
        tables.gr = lapply(1:nrow(z.gr), function(g) z[as.numeric(z.gr[g, c("gr1", "gr2")]), c("n", "total_n")] )
        tables.gr = abind::abind(tables.gr, along=3)

        apply(tables.gr, 3, odds.ratio)
        i.test = mantelhaen.test(tables.gr)
        i.test$p.value


        i.contignency = as.data.frame(z[cgr, c("n", "total_n")])
        i.test = fisher.test(i.contignency)
        i.meta = z[cgr,] %>%
          dplyr::mutate(group_id=1:n()) %>%
          reshape2::melt(measure.vars=c("total_n", "n")) %>%
          dplyr::mutate(variable=paste0("compare_", variable, group_id)) %>%
          tibble::column_to_rownames("variable") %>%
          dplyr::select(value) %>%
          t() %>%
          data.frame()

        cbind(compare_group1=z$tlx_group[cgr[1]], compare_group2=z$tlx_group[cgr[2]], compare_pvalue=i.test$p.value, compare_odds=i.test$estimate, compare_fold=z$n[cgr[1]] / z$n[cgr[2]], i.meta)
      }))
    })(.)) %>%
    data.frame()
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

tlx_macs2 = function(tlx_df, effective_size, maxgap=NULL, qvalue=0.01, pileup=1, extsize=2000, slocal=50000, llocal=10000000, exclude_bait_region=F, exclude_repeats=F, exclude_offtargets=F, exttype=c("along", "symmetrical", "none")) {
  if(exclude_bait_region && !("tlx_is_bait_junction" %in% colnames(tlx_df))) {
    stop("tlx_is_bait_junction is not found in tlx data frame")
  }

  macs2_tlx_df = tlx_df

  if(exclude_offtargets) {
    if(!("tlx_is_offtarget" %in% colnames(macs2_tlx_df))) {
      stop("tlx_is_offtarget is not found in tlx data frame")
    }
    macs2_tlx_df = macs2_tlx_df %>% dplyr::filter(!tlx_is_offtarget)
  }
  if(exclude_repeats) {
    if(!("tlx_repeatmasker_class" %in% colnames(macs2_tlx_df))) {
      stop("tlx_repeatmasker_class is not found in tlx data frame")
    }
    macs2_tlx_df = macs2_tlx_df %>% dplyr::filter(is.na(tlx_repeatmasker_class))
  }

  macs2_tlx_df = macs2_tlx_df %>%
    dplyr::filter(!exclude_bait_region | !tlx_is_bait_junction) %>%
    dplyr::mutate(bed_strand=ifelse(Strand=="1", "-", "+"))

  # @TODO: I think macs does this internally
  if(exttype[1]=="along") {
    macs2_tlx_df = macs2_tlx_df %>% dplyr::mutate(bed_start=ifelse(Strand=="-1", Junction-extsize, Junction-1), bed_end=ifelse(Strand=="-1", Junction, Junction+extsize-1))
  } else {
    if(exttype[1]=="symmetrical") {
      macs2_tlx_df = macs2_tlx_df %>% dplyr::mutate(bed_start=Junction-ceiling(extsize/2), bed_end=Junction+ceiling(extsize/2))
    } else {
      macs2_tlx_df = macs2_tlx_df %>% dplyr::mutate(bed_start=Junction, bed_end=Junction+1)
    }
  }

  if(is.null(maxgap) || maxgap==0 || maxgap=="") {
    maxgap = NULL
  }

  macs_df.all = data.frame()
  for(gr in unique(macs2_tlx_df$tlx_group)) {
    tlx_df.gr = macs2_tlx_df %>% dplyr::filter(tlx_group==gr)

    f_input_bed = tempfile()
    f_control_bed = tempfile()
    # f_input_bed = "tmp/input.bed"
    # f_control_bed = "tmp/control.bed"

    tlx_df.gr %>%
      dplyr::filter(!tlx_control) %>%
      dplyr::select(Rname, bed_start, bed_end, Qname, 0, bed_strand) %>%
      readr::write_tsv(file=f_input_bed, na="", col_names=F)

    if(any(tlx_df.gr$tlx_control)) {
      tlx_df.gr %>%
        dplyr::filter(tlx_control) %>%
        dplyr::select(Rname, bed_start, bed_end, Qname, 0, bed_strand) %>%
        readr::write_tsv(file=f_control_bed, na="", col_names=F)

      log("Running MACS with control")
      macs_df = macs2(name=basename(f_input_bed), sample=f_input_bed, control=f_control_bed, maxgap=maxgap, effective_size=effective_size, extsize=extsize, qvalue=qvalue, slocal=slocal, llocal=llocal, output_dir=dirname(f_input_bed)) %>%
        dplyr::mutate(macs_group=gr)
    } else {
      log("Running MACS without control")
      macs_df = macs2(name=basename(f_input_bed), sample=f_input_bed, maxgap=maxgap, effective_size=effective_size, extsize=extsize, qvalue=qvalue, slocal=slocal, llocal=llocal, output_dir=dirname(f_input_bed)) %>%
        dplyr::mutate(macs_group=gr)
    }

    macs_df.all = rbind(macs_df.all, macs_df %>% dplyr::mutate(macs_group=gr))
  }

  macs_df.all = macs_df.all %>% dplyr::filter(macs_pileup>=pileup)

  macs_df.all
}