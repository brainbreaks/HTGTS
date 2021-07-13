library(dplyr)
library(shiny)
library(readr)
library(GenomeInfoDb)
library(BSgenome)
library(Gviz)
library(stringr)
library(shinyjs)
source("utils.R")
source("tlx.R")
source("graphics.R")

log = function(x) {
  print(paste(">>>", x))
}

server <- function(input, output, session) {
  r <- shiny::reactiveValues(
    hidden_options=T,
    baits_df=NULL,
    tlx_df=NULL,
    offtargets_df=NULL
  )

  observeEvent(input$advanced_hide, {
    shinyjs::toggle(id="advanced_panel")
  })

  observeEvent(input$tlx, {
    print("input$tlx")
    print(basename(input$tlx$name))
    r$tlx_df = tlx_read(input$tlx$datapath, sample=basename(input$tlx$name))
    r$baits_df = tlx_identify_baits(r$tlx_df)
  }, ignoreInit=T)

  observeEvent(input$offtargets, {
    print("offtargets")
    print(basename(input$offtargets$name))
    r$offtargets_df = offtargets_read(input$offtargets$datapath)
  }, ignoreInit=T)


  observeEvent(input$junctions_report, {
    print("input$junctions_report")
    req(!is.null(isolate(r$tlx_df)))

    print("Start calculating")

  })

  observeEvent(input$circos_calculate, {
    print("input$circos_calculate")
    req(!is.null(isolate(r$tlx_df)))
    print("Start calculating")

    junsize = shiny::isolate(input$junsize)
    exclude_repeats = shiny::isolate(input$exclude_repeats)
    exclude_bait_region = shiny::isolate(input$exclude_bait_region)
    bait_region = shiny::isolate(input$bait_region)
    extsize = shiny::isolate(input$extsize)
    qvalue = shiny::isolate(input$qvalue)
    slocal = shiny::isolate(input$slocal)
    llocal = shiny::isolate(input$llocal)
    model = shiny::isolate(input$model)

    cat(file=stderr(), paste0(
      "junsize=", junsize,
      "\nexclude_repeats=", exclude_repeats,
      "\nexclude_bait_region=", exclude_bait_region,
      "\nbait_region=", bait_region,
      "\nextsize=", extsize,
      "\nqvalue=", qvalue,
      "\nslocal=", slocal,
      "\nllocal=", llocal,
      "\nmodel=", model
    ))


    size = pmin(session$clientData$output_circos_output_width, session$clientData$output_circos_output_width)
    output$circos_output = renderImage({

      # setwd("/home/s215v/Workspace/HTGTS/QCReport")
      # input = list(
      #   tlx=list(datapath="test_data/JJ03_B400_012_result.tlx", name="test_data/JJ03_B400_012_result.tlx"),
      #   offtargets=list(datapath="test_data/offtargets.tsv", name="test_data/offtargets.tsv"))
      # r$tlx_df = tlx_read(input$tlx$datapath, sample=basename(input$tlx$name))
      # r$baits_df = tlx_identify_baits(r$tlx_df)
      # r$offtargets_df = offtargets_read(input$offtargets$datapath)
      # session = list(userData=list())
      # size = 500; junsize = 300; exclude_repeats = F; exclude_bait_region = T; bait_region = 500000; extsize = 2000; qvalue = 0.001; slocal = 1000; llocal = 10000000; model = "hg19"

      genome_path = file.path(genomes_path, model, paste0(model, ".fa"))
      cytoband_path = file.path(genomes_path, model, "annotation/cytoBand.txt")
      cytoband = circlize::read.cytoband(cytoband_path, species=model)
      cytoband_df = data.frame(cytoband$chr.len) %>% tibble::rownames_to_column("chrom") %>% dplyr::rename(chrom_length="cytoband.chr.len")

      tlx_df = r$tlx_df %>% dplyr::filter(Rname %in% cytoband_df$chrom)
      tlx_df = tlx_mark_repeats(tlx_df, repeatmasker_df)
      tlx_df = tlx_mark_bait_junctions(tlx_df, bait_region)

      offtarget2bait_df = NULL
      if(!is.null(offtargets)) {
        offtarget2bait_df = join_offtarget2bait(offtargets_df=r$offtargets_df, baits_df=r$baits_df, genome_path=genome_path)
        tlx_df = tlx_mark_offtargets(tlx_df, offtarget2bait_df)
      }

      macs_df = tlx_macs2(tlx_df, extsize=extsize, qvalue=qvalue, slocal=slocal, llocal=llocal, exclude_bait_region=exclude_bait_region, exclude_repeats=exclude_repeats)
      if(!is.null(offtarget2bait_df)) {
        macs_df = macs_df %>%
          dplyr::left_join(offtarget2bait_df, by=c("macs_sample"="bait_sample")) %>%
          dplyr::group_by(macs_chrom, macs_start, macs_end, bait_chrom, bait_start, bait_end, macs_sample) %>%
          dplyr::summarize(macs_is_offtarget=any(offtarget_chrom==macs_chrom & (offtarget_start>=macs_start & offtarget_start<=macs_end | offtarget_end>=macs_start & offtarget_start<=macs_end))) %>%
          dplyr::ungroup() %>%
          dplyr::select(dplyr::matches("^macs_"))
      }

      links_df = macs_df %>% dplyr::inner_join(r$baits_df, by=c("macs_sample"="bait_sample"))
      if("macs_is_offtarget" %in% colnames(links_df)) {
        links_df$color = ifelse(links_df$macs_is_offtarget, "#F46D43FF", "#74ADD180")
      } else {
        links_df$color = "#74ADD180"
      }
      links_df = links_df %>% dplyr::select(chrom1=macs_chrom, start1=macs_start, end1=macs_end, chrom2=bait_chrom, start2=bait_start, end2=bait_end, color)

      session$userData$outfile_svg = tempfile(fileext=".svg")
      print(paste0("Plot circos ", session$userData$outfile_svg, " (w=", size, " h=", size, ")"))
      svg(session$userData$outfile_svg, width=size/72, height=size/72, pointsize=1)
      par(cex=5)
      plot_circos(data=tlx_df %>% dplyr::select(chrom=Rname, start=Rstart, end=Rend), cytoband_path=cytoband_path, links=links_df)
      dev.off()

      print("Finished")

      list(src=normalizePath(session$userData$outfile_svg), contentType='image/svg+xml', width=size, height=size, alt="TLX circos plot")
    }, deleteFile=F)
  })

}