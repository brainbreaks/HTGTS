library(circlize)
library(dplyr)
library(RColorBrewer)

scale_breaks = function(x) {
    breaks = seq(0, max(x), 1e8)
    names(breaks) = paste0(breaks/1e6, "Mb")
    breaks
}

plot_circos = function(data, hits=NULL, species="hg19", circos_bw=5e6-1) {
  cytoband = circlize::read.cytoband(system.file(package = "circlize", "extdata", "cytoBand.txt"), species=species)

  x = data %>%
    dplyr::mutate(
      B_Rstart=floor(B_Rstart/circos_bw)*(circos_bw+1), B_Rend=B_Rstart+circos_bw-1,
      Rstart=floor(Rstart/circos_bw)*(circos_bw+1), Rend=Rstart+circos_bw-1
    ) %>%
    dplyr::group_by(B_Rname, B_Rstart, B_Rend, Rname, Rstart, Rend) %>%
    dplyr::summarise(count=n(), count_sense=sum(Strand>0), count_antisense=sum(Strand<0), is_offtarget=any(is_offtarget)) %>%
    dplyr::ungroup() %>%
    data.frame()

  circos_mincount = 10
  circos_ylim = range(log10(x$count))
  circos_yaxis = seq(circos_ylim[1], circos_ylim[2], 1)
  circos_yaxis_pal = circlize::colorRamp2(seq(circos_ylim[1], circos_ylim[2], length.out=5), rev(RColorBrewer::brewer.pal(5, "Blues")), transparency=0.2)
  circos_offtarget_pal = circlize::colorRamp2(c(-1, 1), RColorBrewer::brewer.pal(11, "RdYlBu")[c(3,9)], transparency=0.5)

  circlize::circos.initializeWithIdeogram(species=species)
  circlize::circos.genomicTrack(x %>% dplyr::select(Rname, Rstart, Rend, dplyr::matches("*")), bg.border=NA, ylim=circos_ylim,
      panel.fun = function(region, value, ...) {
        for(ax in 2:length(circos_yaxis)) {
          circlize::circos.rect(xleft=0, xright=cytoband$chr.len[circlize::get.current.chromosome()], ybottom=circos_yaxis[ax-1], ytop=circos_yaxis[ax], col=circos_yaxis_pal(ax), border="#00000000")
        }
        circlize::circos.yaxis(at=circos_yaxis, labels.cex=0.4)
        circlize::circos.rect(xleft=region$Rstart, xright=region$Rend, ybottom=0, ytop=log10(value$count), col="#333333", border="#333333")
          # circos.genomicLines(region, value, pch = 16, cex = 0.3)
  })
  circlize::circos.genomicTrack(hits, bg.border=NA, ylim=c(0,1), track.height=0.01, cell.padding=c(0,0),
      panel.fun = function(region, value, ...) {
        circlize::circos.rect(xleft=region$macs_start, xright=region$macs_end, ybottom=0, ytop=1, col="#FF3333", border="#FF3333")
  })
  circlize::circos.genomicLink(
    region1=x %>% dplyr::filter(count>=circos_mincount) %>% dplyr::select(chr=B_Rname, start=B_Rstart, end=B_Rend),
    region2=x %>% dplyr::filter(count>=circos_mincount) %>% dplyr::select(chr=Rname, start=Rstart, end=Rend),
    col=circos_offtarget_pal(ifelse(x %>% dplyr::filter(count>=circos_mincount) %>% .$is_offtarget, -1, 1)), border=NA)
  circlize::circos.clear()
}