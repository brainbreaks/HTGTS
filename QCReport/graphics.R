library(circlize)
library(dplyr)
library(RColorBrewer)

scale_breaks = function(x) {
    breaks = seq(0, max(x), 1e8)
    names(breaks) = paste0(breaks/1e6, "Mb")
    breaks
}

scale_breaks_1 = function(x) {
    breaks = seq(floor(min(x)), ceiling(max(x)), 1)
    names(breaks) = as.character(breaks)
    breaks
}

scale_breaks_2 = function(x) {
    breaks = seq(floor(min(x)), ceiling(max(x)), 1)
    names(breaks) = ifelse((breaks %% 2)>0, as.character(breaks), "")
    breaks
}

scale_breaks_sub1k = function(x) {
  breaks = seq(0, max(x), 500)
  names(breaks) = paste0(breaks/1e3, "")
  breaks
}

plot_barcodes = function(tlx_df) {
  tlx_df %>%
    dplyr::mutate(barcode_seq=substr(Seq, 1, 5)) %>%
    dplyr::select(tlx_sample, barcode_seq) %>%
    split(x=as.character(.$barcode_seq), f=.$tlx_sample) %>%
    ggseqlogo()
}

theme_brain = function() {
   ggplot2::theme_classic(base_size=18) +
   ggplot2::theme(legend.key.size=unit(20, 'pt'), plot.title=element_text(hjust=0.5), strip.background=element_blank())
}

plot_homology = function(tlx_df) {
  homology_df = tlx_df %>%
    dplyr::group_by(tlx_group, tlx_sample, B_Rname) %>%
    dplyr::mutate(misprimed_max=max(misprimed-uncut))  %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Pray=B_Qend-Qstart+1, Bait=misprimed_max-misprimed) %>%
    dplyr::mutate(homology_location=ifelse(tlx_is_bait_junction, "Inside bait region", "Outside bait region")) %>%
    dplyr::filter(is.na(tlx_repeatmasker_class)) %>%
    reshape2::melt(measure.vars=c("Pray", "Bait"), variable.name="homology_part", value.name="homology_size") %>%
    #     # dplyr::filter(!tlx_is_bait_junction & homology_size>=0) %>%
    dplyr::group_by(tlx_group, tlx_sample, homology_location, homology_part, homology_size) %>%
    dplyr::summarize(homology_abundance_abs=n()) %>%
    dplyr::group_by(tlx_group, tlx_sample, homology_part, homology_location) %>%
    dplyr::mutate(homology_abundance_rel=homology_abundance_abs/sum(homology_abundance_abs)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(homology_size>=0) %>%
    dplyr::mutate(homology_abundance_rel=ifelse(homology_part=="Bait", -1, 1)*homology_abundance_rel)

  p = list()
  for(gr in unique(homology_df$tlx_group)) {
    homology_df.gr = homology_df %>% dplyr::filter(tlx_group==gr)
    homology_ylim = homology_df.gr %>%
      dplyr::group_by(tlx_group) %>%
      dplyr::summarize(homology_abundance_rel=max(abs(homology_abundance_rel)))
    p[[length(p)+1]] = ggplot(homology_df.gr) +
      geom_rect(aes(ymin=0, ymax=homology_abundance_rel*100, xmin=0, xmax=100), fill="#FF000005", data=homology_ylim) +
      geom_rect(aes(ymin=-homology_abundance_rel*100, ymax=0, xmin=0, xmax=100), fill="#0000FF05", data=homology_ylim) +
      geom_text(aes(x=20, y=homology_abundance_rel*100/2), label="Pray", data=homology_ylim, size=8, hjust=1) +
      geom_text(aes(x=20, y=-homology_abundance_rel*100/2), label="Bait", data=homology_ylim, size=8, hjust=1) +
      geom_hline(yintercept=0) +
      ggplot2::geom_line(aes(x=homology_size, y=homology_abundance_rel*100, color=tlx_sample, group=paste0(tlx_sample, homology_location, homology_part))) +
      ggplot2::labs(x="Deletion size", y="% of junctions", color="Sample") +
      # ggplot2::scale_x_continuous(breaks=scale_breaks_2) +
      ggplot2::scale_y_continuous(labels=abs) +
      ggplot2::coord_cartesian(xlim=c(1, 20)) +
      ggplot2::facet_grid(homology_location~., scales="free") +
      ggplot2::theme_classic(base_size=18) +
      ggplot2::theme(legend.key.size=unit(20, 'pt'), plot.title=element_text(hjust=0.5), strip.background=element_blank())
  }

  gridExtra::grid.arrange(grobs=p, top=textGrob("Junctions homology stats\n(Excluding junctions overlapping with repeats)", gp=gpar(fontsize=20,font=20)))
}

plot_venn = function(x, main, pallete="Pastel2", size=500) {
  p = VennDiagram::venn.diagram(
    main=main,
    x = x, height=size, width=size,
    margin=0.1,
    cat.cex=size/50, cex=size/30, main.cex=size/30,
    lwd=2, lty='blank', fill=RColorBrewer::brewer.pal(length(x), pallete), cat.fontface="bold", filename=NULL)
  grid::grid.draw(p)
}

plot_circos = function(data, title, cytoband_path, annotations=NULL, links=NULL, circos_bw=1e6-1, cex=5) {
  cytoband = circlize::read.cytoband(cytoband_path)
  cytoband_df = data.frame(cytoband$chr.len) %>% tibble::rownames_to_column("chrom") %>% dplyr::rename(chrom_length="cytoband.chr.len")

  unknown_chroms = list("data (chrom1)"=unique(setdiff(data$chrom, cytoband_df$chrom)))
  if(!is.null(annotations)) {
    unknown_chroms[["hits (chrom)"]] = unique(setdiff(hits$chrom, cytoband_df$chrom))
  }
  if(!is.null(links)) {
    unknown_chroms[["links (chrom1)"]] = unique(setdiff(links$chrom1, cytoband_df$chrom))
    unknown_chroms[["links (chrom2)"]] = unique(setdiff(links$chrom2, cytoband_df$chrom))
  }

  if(sum(sapply(unknown_chroms, length))>0) {
    unknown_chroms = sapply(unknown_chroms[sapply(unknown_chroms, length)>0], paste, collapse=",")
    unknown_chroms_err = paste("Unknown chromosomes provided when ploting circos plot...\n", paste(paste0(names(unknown_chroms), ":", unknown_chroms), collapse="\n"))
    stop(unknown_chroms_err)
  }

  links_bw = 2e6
  if(!("color" %in% colnames(links))) {
    links$color = "#F46D4380"
  }
  links_sum = links %>%
    dplyr::mutate(oldstart1=start1, oldend1=end1, oldstart2=start2, oldend2=end2) %>%
    dplyr::mutate(start1=floor(start1/links_bw)*links_bw, end1=start1+links_bw-1) %>%
    dplyr::mutate(start2=floor(start2/links_bw)*links_bw, end2=start2+links_bw-1) %>%
    dplyr::distinct(chrom1, start1, end1, chrom2, start2, end2, color, .keep_all=T) %>%
    dplyr::inner_join(cytoband_df %>% setNames(., paste0(colnames(.), "1")), by="chrom1") %>%
    dplyr::inner_join(cytoband_df %>% setNames(., paste0(colnames(.), "2")), by="chrom2") %>%
    dplyr::mutate(start1=pmax(1, start1), end1=pmin(end1, chrom_length1)) %>%
    dplyr::mutate(start2=pmax(1, start2), end2=pmin(end2, chrom_length2))

  # @todo: calculate pileup here
  data_sum = data %>%
    dplyr::group_by(chrom) %>%
    dplyr::do((function(d){
      h = hist(d$start, plot=F, breaks=c(seq(1, max(d$start), by=circos_bw), max(d$start)))
      data.frame(start=h$breaks[-length(h$breaks)], end=h$breaks[-1]-1, count=h$count, count_log10=ifelse(h$count>0, log10(h$count), 0))
    })(.)) %>%
    dplyr::inner_join(cytoband_df, by="chrom") %>%
    dplyr::mutate(start=pmax(1, start), end=pmin(end, chrom_length-1)) %>%
    data.frame()

  circos_ylim = c(0, max(data_sum$count_log10))
  circos_ylabels = sort(expand.grid(power=10**seq(0, ceiling(circos_ylim[2]), 1), prec=c(1, 2, 5)) %>% dplyr::mutate(y=power*prec) %>% .$y)
  circos_ylabels = circos_ylabels[circos_ylabels<=10**max(data_sum$count_log10)]
  circos_yaxis = log10(circos_ylabels)
  circos_yaxis_pal = circlize::colorRamp2(circos_yaxis, colorRampPalette(rev(RColorBrewer::brewer.pal(5, "Blues")))(length(circos_yaxis)), transparency=0.8)

  # , plotType=c("axis", "labels")
  circos.par("gap.degree" = c(rep(1, nrow(cytoband_df)-1), 5))
  par(cex=cex)
  circlize::circos.initializeWithIdeogram(cytoband=cytoband_path)
  circlize::circos.genomicTrack(data_sum, bg.border=NA, ylim=circos_ylim,
      panel.fun = function(region, value, ...) {
        if(circlize::get.current.chromosome() == cytoband$chromosome[1]) {
          circlize::circos.yaxis(at=circos_yaxis[circos_ylabels>=2], labels=circos_ylabels[circos_ylabels>=2], labels.cex=cex*0.4)
        }
        if(length(circos_yaxis)>1) {
          for(ax in 2:length(circos_yaxis)) {
            circlize::circos.rect(xleft=0, xright=cytoband$chr.len[circlize::get.current.chromosome()], ybottom=circos_yaxis[ax-1], ytop=circos_yaxis[ax], col=circos_yaxis_pal(circos_yaxis[ax-1]), border="#00000000")
          }
        }
        circlize::circos.rect(xleft=region$start, xright=region$end, ybottom=0, ytop=value$count_log10, col="#333333", border="#333333")
  })
  if(!is.null(annotations)) {
    if(is.data.frame(annotations)) {
      annotations = list("annotations"=annotations)
    }
    circlize::circos.genomicTrack(annotations, bg.border=NA, ylim=c(0,1), track.height=0.01, cell.padding=c(0,0),
        panel.fun = function(region, value, ...) {
          circlize::circos.rect(xleft=region$start, xright=region$end, ybottom=0, ytop=1, col="#FF3333", border="#FF3333")
    })
  }
  circlize::circos.genomicLink(
    region1=links_sum %>% dplyr::select(chr=chrom1, start=start1, end=end1),
    region2=links_sum %>% dplyr::select(chr=chrom2, start=start2, end=end2),
    col=links_sum$color, border=NA)
  title(title)

  circlize::circos.clear()
}