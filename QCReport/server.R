library(dplyr)
library(shiny)
library(readr)
library(GenomeInfoDb)
library(BSgenome)
library(Gviz)
library(stringr)
library(shinyjs)
source("utils.R")
source("graphics.R")

log = function(x) {
  print(paste(">>>", x))
}

server <- function(input, output, session) {
  r <- shiny::reactiveValues(
    hidden_options=T,
    tlx_df=NULL
  )

  observeEvent(input$advanced_hide, {
    shinyjs::toggle(id="advanced_panel")
  })

  observeEvent(input$tlx, {
    print("input$tlx")
    r$tlx_df = readr::read_tsv(input$tlx$datapath, comment="#", skip=16, col_names=names(tlx_cols$cols), col_types=tlx_cols) %>%
      dplyr::mutate(tlx_id=1:n())
  }, ignoreInit=T)


  observeEvent(input$circos_calculate, {
    req(!is.null(isolate(r$tlx_df)))
    print("test")

    print("input$circos_calculate")
    junsize = shiny::isolate(input$junsize)
    exclude_repeats = shiny::isolate(input$exclude_repeats)
    exclude_bait_region = shiny::isolate(input$exclude_bait_region)
    bait_region = shiny::isolate(input$bait_region)
    exclude_offtargets = shiny::isolate(input$exclude_offtargets)
    offtargets = shiny::isolate(input$offtargets)
    extsize = shiny::isolate(input$extsize)
    qvalue = shiny::isolate(input$qvalue)
    slocal = shiny::isolate(input$slocal)
    llocal = shiny::isolate(input$llocal)
    model = shiny::isolate(input$model)


    size = pmin(session$clientData$output_circos_output_width, session$clientData$output_circos_output_width)
    output$circos_output = renderImage({

      # setwd("/home/s215v/Workspace/HTGTS/QCReport")
      # size = 500
      # offtargets = list(datapath="test_data/offtargets.tsv")
      # r = list(tlx_df = readr::read_tsv("test_data/JJ03_B400_012_result.tlx", comment="#", skip=16, col_names=names(tlx_cols$cols), col_types=tlx_cols) %>% dplyr::mutate(tlx_id=1:n()))
      # session = list(userData=list())
      # junsize = 300
      # exclude_repeats = F
      # exclude_bait_region = T
      # exclude_offtargets = T
      # bait_region = 500000
      # extsize = 2000
      # qvalue = 0.001
      # slocal = 1000
      # llocal = 10000000
      # model = "hg19"

      cytoband_path = file.path(genomes_path, model, "annotation/cytoBand.txt")
      cytoband = circlize::read.cytoband(cytoband_path, species=model)
      cytoband_df = data.frame(cytoband$chr.len) %>% tibble::rownames_to_column("chrom") %>% dplyr::rename(chrom_length="cytoband.chr.len")

      tlx_df = r$tlx_df %>% dplyr::filter(Rname %in% cytoband_df$chrom)

      if(exclude_repeats) {
        print("Search for overlaps with repeats")
        repeatmasker_ranges = GenomicRanges::makeGRangesFromDataFrame(repeatmasker_df %>% dplyr::mutate(seqnames=repeatmasker_chrom, start=repeatmasker_start, end=repeatmasker_end), keep.extra.columns=T)
        tlx_ranges = GenomicRanges::makeGRangesFromDataFrame(tlx_df %>% dplyr::mutate(seqnames=Rname, start=Rstart, end=Rend), keep.extra.columns=T, ignore.strand=T)
        tlx_df = as.data.frame(IRanges::findOverlaps(tlx_ranges, repeatmasker_ranges)) %>%
          dplyr::inner_join(repeatmasker_df, by=c("subjectHits"="repeatmasker_id")) %>%
          dplyr::group_by(queryHits) %>%
          dplyr::summarise(repeatmasker_name=paste(unique(repeatmasker_name), collapse=", "), repeatmasker_class=paste(unique(repeatmasker_class), collapse=", "), repeatmasker_family=paste(unique(repeatmasker_family), collapse=", "), repeatmasker_length=max(abs(repeatmasker_end-repeatmasker_start))) %>%
          dplyr::right_join(tlx_df, by=c("queryHits"="tlx_id"))
      } else {
        tlx_df$repeatmasker_class = NA_character_
      }


      if(exclude_bait_region) {
        print("Search for bait/bait junctions")
        tlx_df = tlx_df %>% dplyr::mutate(is_bait=B_Rname==Rname & (abs(B_Rstart-Rstart)<=bait_region/2 | abs(Rend-B_Rend)<=bait_region/2))
      } else {
        tlx_df$is_bait = F
      }

      if(!is.null(exclude_offtargets)) {
        print("Search for offtargets")
        tlx_df = tlx_df %>% dplyr::filter(!(B_Rname==Rname & (abs(B_Rstart-Rstart)<=bait_region/2 | abs(Rend-B_Rend)<=bait_region/2)))
        tlx_ranges = GenomicRanges::makeGRangesFromDataFrame(tlx_df %>% dplyr::mutate(seqnames=Rname, start=Rstart, end=Rend), keep.extra.columns=T, ignore.strand=T)

        offtargets_df = readr::read_tsv(offtargets$datapath)
        offtargets_ranges = GenomicRanges::makeGRangesFromDataFrame(offtargets_df %>% dplyr::mutate(seqnames=offtarget_chrom, start=offtarget_start, end=offtarget_end), keep.extra.columns=T)
        tlx_df$is_offtarget = GenomicRanges::countOverlaps(tlx_ranges, offtargets_ranges)>0
      } else {
        tlx_df$is_offtarget = F
      }

      tlx_df = tlx_df %>%
        dplyr::filter(!is_offtarget & !is_bait & is.na(repeatmasker_class))

      f_bed = tempfile()
      print(paste0("Convert to BED (", f_bed, ")"))
      tlx_df %>%
        dplyr::ungroup() %>%
        dplyr::mutate(bed_start=ifelse(Strand=="1", Junction-junsize, Junction-1), bed_end=ifelse(Strand=="1", Junction, Junction+junsize-1), bed_strand=strand_map[Strand]) %>%
        dplyr::select(Rname, bed_start, bed_end, Qname, 0, bed_strand) %>%
        readr::write_tsv(file=f_bed, na="", col_names=F)
      print("Running MACS")
      macs_df = macs2(name=basename(f_bed), sample=f_bed, extsize=extsize, qvalue=qvalue, slocal=slocal, llocal=llocal, output_dir=dirname(f_bed))

      macs_ranges = GenomicRanges::makeGRangesFromDataFrame(macs_df %>% dplyr::mutate(seqnames=macs_chrom, start=macs_start, end=macs_end), keep.extra.columns=T)
      tlx_ranges = GenomicRanges::makeGRangesFromDataFrame(tlx_df %>% dplyr::mutate(seqnames=Rname, start=Rstart, end=Rend), keep.extra.columns=T, ignore.strand=T)
      links_df = as.data.frame(IRanges::mergeByOverlaps(macs_ranges, tlx_ranges)) %>%
        dplyr::mutate(chrom1=macs_chrom, start1=macs_start, end1=macs_end) %>%
        dplyr::group_by(chrom1, start1, end1) %>%
        dplyr::summarize(chrom2=B_Rname[1], start2=floor(mean(B_Rstart)), end2=floor(mean(B_Rend)), color=ifelse(any(is_offtarget), "#74ADD180", "#F46D4380")) %>%
        dplyr::ungroup()

      data_df = tlx_df %>%
        dplyr::select(chrom=Rname, start=Rstart, end=Rend, strand=Strand) %>%
        dplyr::ungroup()

      hits_df = macs_df %>%
        dplyr::select(chrom=macs_chrom, start=macs_start, end=macs_end)

      session$userData$outfile_svg = tempfile(fileext=".svg")
      print(paste0("Plot circos ", session$userData$outfile_svg, " (w=", size, " h=", size, ")"))
      svg(session$userData$outfile_svg, width=size/72, height=size/72, pointsize=1)
      par(cex=5)
      plot_circos(data=data_df, cytoband_path=cytoband_path, hits=hits_df, links=links_df)
      table(tlx_df$Rname=="chrM")
      table(data_df$chrom=="chrM")
      dev.off()

      print("Finished")

      list(src=normalizePath(session$userData$outfile_svg), contentType='image/svg+xml', width=size, height=size, alt="TLX circos plot")
    }, deleteFile=F)
  })

}