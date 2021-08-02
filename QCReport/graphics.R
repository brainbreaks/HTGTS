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

plot_macs2_pileups = function(tlx_df, macs_df) {
    macs_df = macs_df %>% dplyr::mutate(seqnames=macs_chrom, start=macs_start, end=macs_end, macs_group_i=as.numeric(factor(macs_group)))
    macs_ranges = GenomicRanges::makeGRangesFromDataFrame(macs_df, ignore.strand=T, keep.extra.columns=T)

    tlx_df = tlx_df %>%
      dplyr::mutate(tlx_group=factor(tlx_group, unique(as.character(tlx_group))))

    tlxsum_df = tlx_df %>%
      dplyr::group_by(tlx_group, .drop=F) %>%
      dplyr::summarize(total_n=sum(!tlx_is_bait_junction)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(correction=min(total_n)/total_n)

    macs_reduced_df = as.data.frame(GenomicRanges::reduce(macs_ranges)) %>%
      dplyr::mutate(reduced_chrom=seqnames, reduced_start=start, reduced_end=end) %>%
      dplyr::mutate(reduced_hit=paste0(reduced_chrom, ":", reduced_start, "-", reduced_end))
    macs_reduced_ranges = GenomicRanges::makeGRangesFromDataFrame(macs_reduced_df, keep.extra.columns=T)

    facet_df = macs_df %>%
      dplyr::ungroup() %>%
      dplyr::group_by(macs_chrom) %>%
      dplyr::do((function(z){
        data.frame(facet_chrom=z$macs_chrom[1], facet_start=seq(0, 1e9, 5e5)) %>%
          dplyr::mutate(facet_end=facet_start+5e5, seqnames=facet_chrom, start=facet_start, end=facet_end) %>%
          dplyr::mutate(facet_hit=paste0(facet_chrom, ":", facet_start/1e6, "-", facet_end/1e6, "M"))
      })(.))
    facet_ranges = GenomicRanges::makeGRangesFromDataFrame(facet_df, keep.extra.columns=T)

    tlxcov_df = tlx_coverage(tlx_df, group="group")
    tlxcov_ranges = GenomicRanges::makeGRangesFromDataFrame(tlxcov_df %>% dplyr::mutate(seqnames=tlxcov_chrom, start=tlxcov_start, end=tlxcov_end), keep.extra.columns=T)
    tlxcov2macsred_df = as.data.frame(mergeByOverlaps(tlxcov_ranges, macs_reduced_ranges))
    tlxcov2macsred_ranges = GenomicRanges::makeGRangesFromDataFrame(tlxcov2macsred_df %>% dplyr::mutate(seqnames=tlxcov_chrom, start=tlxcov_start, end=tlxcov_end), keep.extra.columns=T)
    ggplot_df = as.data.frame(mergeByOverlaps(tlxcov2macsred_ranges, facet_ranges)) %>%
      dplyr::inner_join(tlxsum_df, by="tlx_group") %>%
      dplyr::mutate(tlxcov_pileup=tlxcov_pileup*correction) %>%
      dplyr::arrange(tlxcov_chrom, tlxcov_start) %>%
      dplyr::group_by(facet_hit) %>%
      dplyr::filter(cumsum(tlxcov_pileup)>0) %>%
      dplyr::ungroup()
    ggplot_reduced_df = ggplot_df %>%
      dplyr::group_by(facet_hit) %>%
      dplyr::mutate(ymin=-max(tlxcov_pileup)*0.1, ymax=-max(tlxcov_pileup)*0.05) %>%
      dplyr::ungroup() %>%
      dplyr::distinct(facet_hit, reduced_chrom, reduced_start, reduced_end, ymin, ymax)

    ggplot(ggplot_df) +
      geom_step(aes(x=tlxcov_start, y=tlxcov_pileup, color=tlx_group)) +
      geom_rect(aes(xmin=reduced_start, xmax=reduced_end, ymin=ymin, ymax=ymax), data=ggplot_reduced_df) +
      labs(x="", y="", color="Group") +
      facet_wrap(~facet_hit, scales="free") +
      ggplot2::scale_x_continuous(breaks=scale_breaks) +
      theme_brain() +
      theme(legend.position="bottom") +
      guides(fill=guide_legend(nrow=2, byrow=T))
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

plot_circos = function(input, control, title, cytoband_path, chromosomes, annotations=NULL, links=NULL, circos_bw=1e6-1, cex=5, colors=c(neutral="#999999", input="#FFCB00", control="#FF5700")) {
  cytoband = circlize::read.cytoband(cytoband_path)
  cytoband_df = data.frame(cytoband$chr.len) %>% tibble::rownames_to_column("chrom") %>% dplyr::rename(chrom_length="cytoband.chr.len")

  if(!is.null(control)) {
    has_control = nrow(control)>0
  } else {
    has_control = F
  }

  unknown_chroms = list("data (chrom1)"=unique(setdiff(unique(c(input$chrom, control$chrom)), cytoband_df$chrom)))

  if(!is.null(annotations)) {
    if(is.data.frame(annotations)) {
      annotations = list("annotations"=annotations)
    }

    for(n in names(annotations)) {
      unknown_chroms[[paste0("annotations (", n, ")")]] = unique(setdiff(annotations[[n]]$chrom, cytoband_df$chrom))
    }
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

  if(has_control) {
    scale = c(input=1, control=nrow(input)/nrow(control))
    data_sum = rbind(input %>% dplyr::mutate(circos_signal="input"), control %>% dplyr::mutate(circos_signal="control"))
  } else {
    scale = c(input=1)
    data_sum = input %>% dplyr::mutate(circos_signal="input")
  }
  data_sum = data_sum %>%
    dplyr::inner_join(cytoband_df, by="chrom") %>%
    dplyr::group_by(circos_signal, chrom_length, chrom) %>%
    dplyr::do((function(d){
      # x_ranges = GenomicRanges::makeGRangesFromDataFrame(d, ignore.strand=T)
      # cov_ranges = as(GenomicRanges::coverage(x_ranges), "GRanges")
      # as.data.frame(cov_ranges) %>%
      #   dplyr::mutate(count=score, count_log10=ifelse(count>0, log10(count), 0)) %>%
      #   dplyr::select(start, end, count=score, count_log10)

      h = hist(d$start, plot=F, breaks=c(seq(1, d$chrom_length[1], by=circos_bw), d$chrom_length[1]))
      data.frame(start=h$breaks[-length(h$breaks)], end=h$breaks[-1]-1, count=scale[d$circos_signal[1]]*h$count, count_log10=ifelse(h$count>0, log10(h$count), 0))
    })(.)) %>%
    dplyr::mutate(start=pmax(1, start), end=pmin(end, chrom_length-1)) %>%
    dplyr::mutate(circos_col=paste0("count_log10.", circos_signal)) %>%
    dplyr::filter(end<=chrom_length) %>%
    reshape2::dcast(chrom + start + end + chrom_length ~ circos_col, value.var="count_log10")

  if(has_control) {
    circos_ymax = max(c(data_sum$count_log10.input, data_sum$count_log10.control))
  } else {
    data_sum$count_log10.input
  }

  circos_ylim = c(0, circos_ymax)
  circos_ylabels = sort(expand.grid(power=10**seq(0, ceiling(circos_ylim[2]), 1), prec=c(1, 2, 5)) %>% dplyr::mutate(y=power*prec) %>% .$y)
  circos_ylabels = circos_ylabels[circos_ylabels<=10**circos_ymax]
  circos_yaxis = log10(circos_ylabels)
  circos_yaxis_pal = circlize::colorRamp2(circos_yaxis, colorRampPalette(rev(RColorBrewer::brewer.pal(5, "Blues")))(length(circos_yaxis)), transparency=0.8)

  circos.par("gap.degree"=c(rep(1, length(chromosomes)-1), 5))
  par(cex=cex, cex.main=cex)
  circlize::circos.initializeWithIdeogram(cytoband=cytoband_path, chromosome.index=chromosomes, plotType=c("axis", "labels"))
  circlize::circos.genomicTrack(data_sum, bg.border=NA, ylim=circos_ylim,
      panel.fun = function(region, value, ...) {
        rr<<-region
        vv<<-value
        if(circlize::get.current.chromosome() == cytoband$chromosome[1]) {
          circlize::circos.yaxis(at=circos_yaxis[circos_ylabels>=2], labels=circos_ylabels[circos_ylabels>=2], labels.cex=cex*0.4)
        }
        if(length(circos_yaxis)>1) {
          for(ax in 2:length(circos_yaxis)) {
            circlize::circos.rect(xleft=0, xright=cytoband$chr.len[circlize::get.current.chromosome()], ybottom=circos_yaxis[ax-1], ytop=circos_yaxis[ax], col=circos_yaxis_pal(circos_yaxis[ax-1]), border="#00000000")
          }
        }

        if(has_control) {
          count_log10.pmin = pmin(value$count_log10.input, value$count_log10.control)
          circlize::circos.rect(xleft=region$start, xright=region$end+1, ybottom=0, ytop=count_log10.pmin, col=colors["neutral"], border=colors["neutral"])

          f = value$count_log10.input>count_log10.pmin
          if(any(f)) circlize::circos.rect(xleft=region$start[f], xright=region$end[f]+1, ybottom=count_log10.pmin[f], ytop=value$count_log10.input[f], col=colors["input"], border=colors["input"])

          f = value$count_log10.control>count_log10.pmin
          if(any(f)) circlize::circos.rect(xleft=region$start[f], xright=region$end[f]+1, ybottom=count_log10.pmin[f], ytop=value$count_log10.control[f], col=colors["control"], border=colors["control"])
        } else {
          circlize::circos.rect(xleft=region$start, xright=region$end+1, ybottom=0, ytop=value$count_log10.input, col="#999999", border="#999999")
        }
  })

  if(!is.null(annotations) && length(annotations)>0) {
    circlize::circos.genomicTrack(annotations, bg.border=NA, ylim=c(0,1), track.height=0.04, cell.padding=c(0,0),
        panel.fun = function(region, value, ...) {
          circlize::circos.rect(xleft=region$start, xright=region$end, ybottom=0, ytop=1, col="#330000", border="#330000")
    })
  }

  if(!is.null(links)) {
    if(!("color" %in% colnames(links))) {
      links$color = "#F46D4380"
    }


    links_sum = links %>%
      dplyr::filter(chrom1 %in% chromosomes & chrom2 %in% chromosomes) %>%
      dplyr::inner_join(cytoband_df %>% setNames(., paste0(colnames(.), "1")), by="chrom1")  %>%
      dplyr::inner_join(cytoband_df %>% setNames(., paste0(colnames(.), "2")), by="chrom2")  %>%
      dplyr::rowwise() %>%
      dplyr::do((function(d){
        dd<<-d
        breaks = seq(1, d$chrom_length1[1], by=circos_bw)
        data.frame(
          chrom1=d$chrom1,
          start1=breaks[which(breaks>d$start1)[1]-1],
          end1=breaks[rev(which(breaks<d$end1))[1]+1],
          chrom2=d$chrom2,
          start2=breaks[which(breaks>d$start2)[1]-1],
          end2=breaks[rev(which(breaks<d$end2))[1]+1],
          color=d$color
        )
      })(.))

    if(nrow(links_sum) > 0) {
      circlize::circos.genomicLink(
        region1=links_sum %>% dplyr::select(chr=chrom1, start=start1, end=end1),
        region2=links_sum %>% dplyr::select(chr=chrom2, start=start2, end=end2),
        col=links_sum$color, border=NA)
    }
  }
  title(title)
  legend("bottomleft", title="Reads", legend=names(colors), fill=colors, xjust=1, yjust=1)

  circlize::circos.clear()
}