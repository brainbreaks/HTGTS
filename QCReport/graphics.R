library(circlize)
library(dplyr)
library(RColorBrewer)

scale_breaks = function(x) {
    breaks = seq(0, max(x), 1e8)
    names(breaks) = paste0(breaks/1e6, "Mb")
    breaks
}

plot_circos = function(data, cytoband_path, hits=NULL, links=NULL, circos_bw=5e6-1) {
  cytoband = circlize::read.cytoband(cytoband_path)
  cytoband_df = data.frame(cytoband$chr.len) %>% tibble::rownames_to_column("chrom") %>% dplyr::rename(chrom_length="cytoband.chr.len")

  unknown_chroms = list("data (chrom1)"=unique(setdiff(data$chrom, cytoband_df$chrom)))
  if(!is.null(hits)) {
    unknown_chroms[["hits (chrom)"]] = unique(setdiff(hits$chrom, cytoband_df$chrom))
  }
  if(!is.null(hits)) {
    unknown_chroms[["links (chrom1)"]] = unique(setdiff(links$chrom1, cytoband_df$chrom))
    unknown_chroms[["links (chrom2)"]] = unique(setdiff(links$chrom2, cytoband_df$chrom))
  }

  if(sum(sapply(unknown_chroms, length))>0) {
    unknown_chroms = sapply(unknown_chroms[sapply(unknown_chroms, length)>0], paste, collapse=",")
    unknown_chroms_err = paste("Unknown chromosomes provided when ploting circos plot...\n", paste(paste0(names(unknown_chroms), ":", unknown_chroms), collapse="\n"))
    writeLines(unknown_chroms_err)
    return()
  }

  data_sum = data %>%
    dplyr::mutate(start=floor(start/circos_bw)*(circos_bw+1), Rend=start+circos_bw-1) %>%
    dplyr::group_by(chrom, start, end) %>%
    dplyr::summarise(count=n()) %>%
    dplyr::ungroup() %>%
    dplyr::inner_join(cytoband_df, by="chrom") %>%
    dplyr::mutate(start=pmax(1, start), end=pmin(end, chrom_length-1)) %>%
    dplyr::select(chrom, start, end, count) %>%
    dplyr::filter(count>1) %>%
    dplyr::right_join(cytoband_df, by="chrom") %>%
    dplyr::mutate(count=tidyr::replace_na(count, 0), start=tidyr::replace_na(start, 1), end=tidyr::replace_na(end, 1)) %>%
    dplyr::mutate(count_log10=ifelse(is.na(count) | count==0, 0, log10(count))) %>%
    data.frame()

  circos_ylim = c(0, max(data_sum$count_log10))
  circos_yaxis = seq(floor(circos_ylim[1]), ceiling(circos_ylim[2]), 1)
  circos_yaxis_pal = circlize::colorRamp2(seq(floor(circos_ylim[1]), ceiling(circos_ylim[2]), length.out=5), rev(RColorBrewer::brewer.pal(5, "Blues")), transparency=0.8)

  circlize::circos.initializeWithIdeogram(cytoband=cytoband_path)
  circlize::circos.genomicTrack(data_sum, bg.border=NA, ylim=circos_ylim,
      panel.fun = function(region, value, ...) {
        if(circlize::get.current.chromosome() == cytoband$chromosome[1]) {
          circlize::circos.yaxis(at=circos_yaxis, labels.cex=0.4)
        }
        if(length(circos_yaxis)>1) {
          for(ax in 2:length(circos_yaxis)) {
            circlize::circos.rect(xleft=0, xright=cytoband$chr.len[circlize::get.current.chromosome()], ybottom=circos_yaxis[ax-1], ytop=circos_yaxis[ax], col=circos_yaxis_pal(circos_yaxis[ax-1]), border="#00000000")
          }
        }
        circlize::circos.rect(xleft=region$start, xright=region$end, ybottom=0, ytop=value$count_log10, col="#333333", border="#333333")
  })
  circlize::circos.genomicTrack(hits, bg.border=NA, ylim=c(0,1), track.height=0.01, cell.padding=c(0,0),
      panel.fun = function(region, value, ...) {
        circlize::circos.rect(xleft=region$start, xright=region$end, ybottom=0, ytop=1, col="#FF3333", border="#FF3333")
  })
  circlize::circos.genomicLink(
    region1=links %>% dplyr::select(chr=chrom1, start=start1, end=end1),
    region2=links %>% dplyr::select(chr=chrom2, start=start2, end=end2),
    col=links$color, border=NA)
  circlize::circos.clear()
}