source("utils.R")
library(readr)
library(dplyr)
library(GenomicRanges)
library(data.table)
library(rstatix)

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

tlx_read = function(path, sample, group="", group_i=1, control=F) {
  readr::read_tsv(path, comment="#", skip=16, col_names=names(tlx_cols$cols), col_types=tlx_cols) %>%
      dplyr::mutate(tlx_sample=sample, tlx_path=path, tlx_group=group, tlx_group_i=group_i, tlx_control=control)
}

tlx_read_many = function(samples_df) {
  tlx_df.all = data.frame()
  for(f in 1:nrow(samples_df)) {
    log("Reading tlx file ", f, ":",  samples_df$path[f])
    tlx_df.f = tlx_read(samples_df$path[f], sample=samples_df$sample[f], control=samples_df$control[f], group=samples_df$group[f], group_i=samples_df$group_i[f])
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


tlx_test_hits = function(tlx_df, hits_ranges, paired_samples=T, paired_control=T, extsize=10000, exttype="along") {

  # save(tlx_df, hits_ranges, file="d.rda")
  # load("d.rda")
  # macs_ranges=hits_ranges

  hits_df = as.data.frame(hits_ranges) %>% dplyr::mutate(compare_chrom=seqnames, compare_start=start, compare_end=end)
  hits_ranges = GenomicRanges::makeGRangesFromDataFrame(hits_df, keep.extra.columns=T)

  tlxsum_df = tlx_df %>%
    dplyr::group_by(tlx_sample, .drop=F) %>%
    dplyr::summarize(compare_total=sum(!tlx_is_bait_junction)) %>%
    dplyr::ungroup()

  # Prepare overlap counts table (add compare_n=0 for missing entries)
  counts_df_incomplete = as.data.frame(IRanges::mergeByOverlaps(hits_ranges, tlx_ranges)) %>%
    dplyr::rename(compare_group="tlx_group", compare_group_i="tlx_group_i", compare_sample="tlx_sample") %>%
    dplyr::group_by(compare_chrom, compare_start, compare_end, compare_group, compare_group_i, compare_sample, tlx_control, .drop=F) %>%
    dplyr::summarize(compare_n=n())
  counts_df = dplyr::bind_rows(
    counts_df_incomplete,
    hits_df %>%
      dplyr::select(compare_chrom, compare_start, compare_end) %>%
      tidyr::crossing(tlx_df %>% dplyr::distinct(compare_group=tlx_group, compare_group_i=tlx_group_i, compare_sample=tlx_sample, tlx_control)) %>%
      dplyr::mutate(compare_n=0)) %>%
    dplyr::distinct(compare_chrom, compare_start, compare_end, compare_group, compare_group_i, compare_sample, tlx_control, .keep_all=T) %>%
    dplyr::inner_join(tlxsum_df, by=c("compare_sample"="tlx_sample")) %>%
    data.frame()

  #
  # Calculate breaks count adjusted with control (by substracting control breaks)
  #
  counts_df.input = counts_df %>% dplyr::filter(!tlx_control)
  counts_df.control = counts_df %>% dplyr::filter(tlx_control)
  if(paired_control) {
   normcounts_df = counts_df.input %>%
     dplyr::inner_join(counts_df.control %>% select(compare_chrom, compare_start, compare_end, compare_group, compare_group_i, compare_total.control=compare_total, compare_n.control=compare_n), by=c("compare_chrom", "compare_start", "compare_end", "compare_group", "compare_group_i")) %>%
     dplyr::mutate(compare_n.control_adj=compare_n.control*(compare_total/compare_total.control),  compare_n.norm=compare_n-compare_n.control_adj, compare_frac.norm=compare_n.norm/compare_total, compare_frac=compare_n/compare_total) %>%
     dplyr::arrange(compare_group, compare_group_i)
  } else {
   normcounts_df = counts_df.input %>%
     dplyr::left_join(counts_df.control %>% dplyr::group_by(compare_chrom, compare_start, compare_end, compare_group) %>% summarize(compare_total.control=sum(compare_total), compare_n.control=sum(compare_n)), by=c("compare_chrom", "compare_start", "compare_end", "compare_group")) %>%
     dplyr::mutate(compare_total.control=ifelse(!is.na(compare_total.control), compare_total.control, compare_total), compare_n.control=tidyr::replace_na(compare_n.control, 0)) %>%
     dplyr::mutate(compare_n.control_adj=compare_n.control*(compare_total/compare_total.control),  compare_n.norm=compare_n-compare_n.control_adj, compare_frac.norm=compare_n.norm/compare_total, compare_frac=compare_n/compare_total) %>%
     dplyr::arrange(compare_group, compare_group_i)
  }

  if(paired_samples)
  {
    normcounts_sum_df = normcounts_df %>%
      dplyr::select(compare_chrom, compare_start, compare_end, compare_group, compare_group_i, compare_n.norm, compare_n, compare_total, compare_n.control, compare_total.control, compare_n.control_adj)
  } else {
    normcounts_sum_df = normcounts_df %>%
      dplyr::group_by(compare_chrom, compare_start, compare_end, compare_group) %>%
      dplyr::summarize(compare_group_i=1, compare_n=sum(compare_n), compare_n.norm=sum(compare_n.norm), compare_total=sum(compare_total), compare_n.control=sum(compare_n.control), compare_total.control=sum(compare_total.control)) %>%
      dplyr::mutate(compare_n.control_adj=compare_n.control*(compare_total/compare_total.control))
  }

  z_sum.test = as.data.frame(t(apply(combn(unique(normcounts_sum_df$compare_group), 2), 2, sort))) %>%
    dplyr::rename(compare_group1="V1", compare_group2="V2") %>%
    dplyr::inner_join(normcounts_sum_df %>% dplyr::rename(compare_n.norm1="compare_n.norm", compare_n1="compare_n", compare_total1="compare_total", compare_n.control1="compare_n.control", compare_n.control_adj1="compare_n.control_adj", compare_total.control1="compare_total.control"), by=c("compare_group1"="compare_group")) %>%
    dplyr::inner_join(normcounts_sum_df %>% dplyr::rename(compare_n.norm2="compare_n.norm", compare_n2="compare_n", compare_total2="compare_total", compare_n.control2="compare_n.control", compare_n.control_adj2="compare_n.control_adj", compare_total.control2="compare_total.control"), by=c("compare_group2"="compare_group", "compare_group_i", "compare_chrom", "compare_start", "compare_end")) %>%
    # dplyr::filter(compare_group1=="group 1" & compare_group2=="group 2" & compare_chrom=="chr6" & compare_start==77128688 & compare_end==77564403) %>%
    dplyr::group_by(compare_group1, compare_group2, compare_chrom, compare_start, compare_end) %>%
    dplyr::do((function(z){
      zz<<-z

      z.groups = c(z$compare_group1[1], z$compare_group2[1])
      z.fold = mean(z$compare_n1/z$compare_total1) / mean(z$compare_n2/z$compare_total2)

      # Reapeated measures ANOVA
      z.test_data = z %>%
        reshape2::melt(measure.vars=c("compare_n1", "compare_n.control_adj1", "compare_n2", "compare_n.control_adj2")) %>%
        dplyr::select(compare_group_i, treatment=variable, breaks=value) %>%
        dplyr::mutate(group=z.groups[as.numeric(gsub(".*([0-9])$", "\\1", treatment))], treatment=gsub("([0-9])$", "", treatment)) %>%
        dplyr::mutate(treatment=c(compare_n.control="control", compare_n="treatment", compare_n.control_adj="control")[treatment]) %>%
        dplyr::mutate(group=factor(group), treatment=factor(treatment), compare_group_i=factor(compare_group_i)) %>%
        tibble::tibble()
      z.aov = rstatix::anova_test(data=z.test_data, dv=breaks, wid=compare_group_i, within=c(treatment, group))
      z.aov_pval = data.frame(z.aov) %>% dplyr::filter(Effect=="group") %>% .$p

      i.contignency = lapply(split(z, 1:nrow(z)), function(y) matrix(as.numeric(y[c("compare_n1", "compare_n2", "compare_total1", "compare_total2")]), ncol=2))
      if(length(i.contignency)>=2) {
        i.contignency = abind::abind(i.contignency, along=3)
        i.test = mantelhaen.test(i.contignency)
      } else {
        i.test = fisher.test(i.contignency[[1]])
      }

      z %>%
        dplyr::slice(1) %>%
        dplyr::mutate(compare_pvalue=i.test$p.value, compare_odds=i.test$estimate, compare_aov_pvalue=z.aov_pval, compare_fold=z.fold) %>%
        dplyr::select(compare_group1, compare_group2, compare_chrom, compare_start, compare_end, compare_pvalue, compare_odds, compare_aov_pvalue, compare_fold)
    })(.)) %>%
    data.frame()

  list(test=z_sum.test, data=normcounts_df)
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