library(dplyr)
library(shiny)
library(dqshiny)
library(readr)
library(GenomeInfoDb)
library(BSgenome)
library(Gviz)
library(stringr)
library(shinyjs)

log = function(x) {
  print(paste(">>>", x))
}

server <- function(input, output, session) {
  r <- shiny::reactiveValues(
    hidden_options=T
  )

  observeEvent(input$advanced_hide, {
    shinyjs::toggle(id="advanced_panel")
  })


  observeEvent(input$overview_calculate, {
    input = list(
      offtargets=list(datapath="test_data/offtargets.tsv"),
      tlx=list(datapath="test_data/JJ03_B400_012_result.tlx"),
      bait_region=500000,
      exclude_bait_region=T,
      exclude_repeats=T,
      junsize=300
    )

    repeatmasker_ranges = GenomicRanges::makeGRangesFromDataFrame(repeatmasker_df %>% dplyr::mutate(seqnames=repeatmasker_chrom, start=repeatmasker_start, end=repeatmasker_end), keep.extra.columns=T)


    tlx_df = readr::read_tsv(input$tlx$datapath, comment="#", skip=16, col_names=names(tlx_cols$cols), col_types=tlx_cols) %>%
      dplyr::mutate(tlx_id=1:n())
    tlx_ranges = GenomicRanges::makeGRangesFromDataFrame(tlx_df %>% dplyr::mutate(seqnames=Rname, start=Rstart, end=Rend), keep.extra.columns=T, ignore.strand=T)


      # dplyr::mutate(sample_file=sample_file) %>%
      # dplyr::inner_join(samples_df, by="sample_file") %>%
      # dplyr::inner_join(baits_df, by="bait_id") %>%
      # dplyr::left_join(offtargets_df, by="bait_id") %>%
    tlx2repeatmasker_df = as.data.frame(IRanges::findOverlaps(tlx_ranges, repeatmasker_ranges)) %>%
      dplyr::inner_join(repeatmasker_df, by=c("subjectHits"="repeatmasker_id")) %>%
      dplyr::group_by(queryHits) %>%
      dplyr::summarise(repeatmasker_name=paste(unique(repeatmasker_name), collapse=", "), repeatmasker_class=paste(unique(repeatmasker_class), collapse=", "), repeatmasker_family=paste(unique(repeatmasker_family), collapse=", "), repeatmasker_length=max(abs(repeatmasker_end-repeatmasker_start))) %>%
      dplyr::right_join(tlx_df, by=c("queryHits"="tlx_id"))

    if(input$exclude_repeats) {
      tlx2repeatmasker_df = tlx2repeatmasker_df %>% dplyr::filter(is.na(repeatmasker_class))
    }
    if(input$exclude_bait_region) {
      tlx2repeatmasker_df = tlx2repeatmasker_df %>% dplyr::filter(!(B_Rname==Rname & (abs(B_Rstart-Rstart)<=input$bait_region/2 | abs(Rend-B_Rend)<=input$bait_region/2)))
    }
    if(input$exclude_offtargets) {
      tlx2repeatmasker_df = tlx2repeatmasker_df %>% dplyr::filter(!(B_Rname==Rname & (abs(B_Rstart-Rstart)<=input$bait_region/2 | abs(Rend-B_Rend)<=input$bait_region/2)))
      tlx2repeatmasker_ranges = GenomicRanges::makeGRangesFromDataFrame(tlx2repeatmasker_df %>% dplyr::mutate(seqnames=Rname, start=Rstart, end=Rend), keep.extra.columns=T, ignore.strand=T)

      offtargets_df = readr::read_tsv(input$offtargets$datapath)
      offtargets_ranges = GenomicRanges::makeGRangesFromDataFrame(offtargets_df %>% dplyr::mutate(seqnames=offtarget_chrom, start=offtarget_start, end=offtarget_end), keep.extra.columns=T)
      tlx2repeatmasker_df$is_offtarget = GenomicRanges::countOverlaps(tlx2repeatmasker_ranges, offtargets_ranges)>0
      tlx2repeatmasker_df = tlx2repeatmasker_df %>% dplyr::filter(is_offtarget)
    }

    #
    # Convert to BED file
    #
    tlx2repeatmasker_bed = tlx2repeatmasker_df %>%
      dplyr::ungroup() %>%
      dplyr::mutate(bed_start=ifelse(Strand=="1", Junction-input$junsize, Junction-1), bed_end=ifelse(Strand=="1", Junction, Junction+input$junsize-1), bed_strand=ifelse(Strand=="1", "-", "+")) %>%
      dplyr::select(Rname, bed_start, bed_end, Qname, 0, bed_strand)

    f_bed = file.path("tmp", basename(gsub("\\.tlx$", ".bed", input$tlx$datapath)))

    dir.create("tmp", showWarnings=F)
    readr::write_tsv(tlx2repeatmasker_bed, file=f_bed, na="", col_names=F)

    islands_df = macs2(f_filtered, f_bed, extsize=input$extsize, qvalue=input$qvalue, slocal=input$slocal, llocal=input$llocal)  %>%
      dplyr::mutate(sample_file=sample_file)
    islands_all = dplyr::bind_rows(islands_all, islands_df)

  })

}