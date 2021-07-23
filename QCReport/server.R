library(dplyr)
library(shiny)
library(shinyjs)
library(readr)
library(GenomeInfoDb)
library(BSgenome)
library(Gviz)
library(stringr)
library(ggplot2)
library(forcats)
library(reactlog)
source("utils.R")
source("tlx.R")
source("graphics.R")

server <- function(input, output, session) {
  # @todo: properly load default
  load(paste0("tmp/hg19_repeatmasker_df.rda"))
  r <- shiny::reactiveValues(
    hidden_options=T,
    baits_df=NULL,
    tlx_df=NULL,
    offtargets_df=NULL,
    tlx_num=1,
    repeatmasker_df=repeatmasker_df,
    observers=list()
  )

  observeEvent(input$advanced_hide, {
    shinyjs::toggle(id="advanced_panel")
  })

  # observeEvent(input$genome, {
  #   log("input$genome")
  #   # load(paste0("tmp/", input$genome, "_repeatmasker_df.rda"))
  #   # r$repeatmasker_df = repeatmasker_df
  # })


  output$tlx_uploads = shiny::insertUI("tlx_uploads", where="afterEnd", ui={
    log("output$tlx_uploads")
    # req(r$tlx_num)

    tlx_num = r$tlx_num
    observers = isolate(r$observers)

    for(i in 1:tlx_num) {
      if(!(paste0("tlx", i) %in% names(observers))) {
        observers[[paste0("tlx", i)]] = observeEvent(input[[paste0("tlx", i)]], {
          log("input$tlx", i)
          inp = input[[paste0("tlx", i)]]
          r$tlx_df = tlx_read(inp$datapath, sample=basename(inp$name))
          log("TLX file ", basename(inp$name), " processed...")
        }, ignoreInit=T, autoDestroy=T)

        observers[[paste0("tlx_add", i)]] = observeEvent(input[[paste0("tlx_add", i)]], {
          log("input$tlx_add", i)
          r$tlx_num = tlx_num + 1
        }, ignoreInit=T)
      }
    }

    r$observers

    do.call(fluidRow, lapply(1:tlx_num, function(i) {
      log("tlx", i)

      inpf = shiny::fileInput(paste0("tlx", i), label="", placeholder="No file selected")
      if(i==tlx_num) {
        inp_add = shiny::actionButton(paste0("tlx_add", i), label="+")
        inp_del = shiny::actionButton(paste0("tlx_del", i), label="-")

        if(i > 1) return(list(shiny::column(10, inpf),  shiny::column(1, inp_add), shiny::column(1, inp_del)))
        else return(list(shiny::column(10, inpf),  shiny::column(2, inp_add)))
      } else {
        return(shiny::column(10, inpf))
      }
    }))
  })



  observeEvent(input$offtargets, {
    log("input$offtargets")
    log(basename(input$offtargets$name))
    r$offtargets_df = offtargets_read(input$offtargets$datapath)
  }, ignoreInit=T)

  eventReactive(input$calculate, {
    req(!is.null(isolate(r$tlx_df)))
    log("input$calculate")

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

    # setwd("/home/s215v/Workspace/HTGTS/QCReport")
    # r = list()
    # input = list(
    #   tlx=list(datapath="test_data/JJ03_B400_012_result.tlx", name="test_data/JJ03_B400_012_result.tlx"),
    #   offtargets=list(datapath="test_data/offtargets.tsv", name="test_data/offtargets.tsv"))
    # r$tlx_df = tlx_read(input$tlx$datapath, sample=basename(input$tlx$name))
    # r$baits_df = tlx_identify_baits(r$tlx_df)
    # r$offtargets_df = offtargets_read(input$offtargets$datapath)
    # session = list(userData=list())
    # size = 500; junsize = 300; exclude_repeats = F; exclude_bait_region = T; bait_region = 500000; extsize = 2000; qvalue = 0.001; slocal = 1000; llocal = 10000000; model = "hg19"; genomes_path = "/home/s215v/Workspace/HTGTS/genomes"
    # samples_df = rbind(samples_df, as.data.frame(do.call(rbind, list(
    #   c(path="test_data/JJ01_B400_012_result.tlx", sample="Control"),
    #   c(path="test_data/JJ02_B400_012_result.tlx", sample="4Gr"),
    #   c(path="test_data/JJ03_B400_012_result.tlx", sample="4Gr + PARPi")
    # ))))

    #
    # Show microhomology profile
    #
    # tlx_df = tlx_read_many(samples_df)


    #
    # Vivien's
    #
    # setwd("/home/s215v/Workspace/HTGTS/QCReport")
    # r = list(); width = 500; height=500; junsize = 300; exclude_repeats = F; exclude_bait_region = T; bait_region = 500000; extsize = 2000; qvalue = 0.001; slocal = 1000; llocal = 10000000; model = "mm10"; genomes_path = "/home/s215v/Workspace/HTGTS/genomes"
    # session = list(userData=list(
    #   repeats_summary_svg="Vivien/reports/repeats_summary.svg",
    #   junctions_venn_svg = "Vivien/reports/junctions_venn.svg",
    #   homology_svg = "Vivien/reports/homology_profile.svg",
    #   circos_svg = "Vivien/reports/circos.svg"
    # ))
    # load(paste0("tmp/mm10_repeatmasker_df.rda"))
    # r$repeatmasker_df = repeatmasker_df
    # r$baits_df = as.data.frame(bed_read("Vivien/baits.bed")) %>%
    #   setNames(paste0("bait_", colnames(.))) %>%
    #   dplyr::rename(bait_chrom="bait_seqnames") %>%
    #   tidyr::crossing(samples_df %>% dplyr::select(bait_sample=sample))
    # samples_df = data.frame(path=list.files("Vivien", pattern="*.tlx", full.names=T)) %>%
    #   dplyr::mutate(sample=gsub("_B400_0", " / ", gsub("_result.tlx", "", basename(path))))
    # r$offtargets_df = offtargets_read("Vivien/offtargets.bed")
    # tlx_df = tlx_read_many(samples_df)
    # tlx_df = tlx_remove_rand_chromosomes(tlx_df)
    # tlx_df = tlx_mark_repeats(tlx_df, r$repeatmasker_df)
    # tlx_df = tlx_mark_bait_chromosome(tlx_df)
    # tlx_df = tlx_mark_bait_junctions(tlx_df, bait_region)
    log("SHOULD NOT BE HERE")

    genome_path = file.path(genomes_path, model, paste0(model, ".fa"))
    cytoband_path = file.path(genomes_path, model, "annotation/cytoBand.txt")

    offtarget2bait_df = NULL
    if(!is.null(r$offtargets_df)) {
      offtarget2bait_df = join_offtarget2bait(offtargets_df=r$offtargets_df, baits_df=r$baits_df, genome_path=genome_path)
      tlx_df = tlx_mark_offtargets(tlx_df, offtarget2bait_df)
    }

    macs_df = tlx_macs2(tlx_df, junsize=junsize, extsize=extsize, qvalue=qvalue, slocal=slocal, llocal=llocal, exclude_bait_region=exclude_bait_region, exclude_repeats=exclude_repeats)
    if(!is.null(offtarget2bait_df)) {
      macs_df = macs_df %>%
        dplyr::left_join(offtarget2bait_df, by=c("macs_sample"="bait_sample")) %>%
        dplyr::group_by(macs_chrom, macs_start, macs_end, bait_chrom, bait_start, bait_end, macs_sample) %>%
        dplyr::summarize(macs_is_offtarget=any(offtarget_chrom==macs_chrom & (offtarget_start>=macs_start & offtarget_start<=macs_end | offtarget_end>=macs_start & offtarget_start<=macs_end))) %>%
        dplyr::ungroup() %>%
        dplyr::select(dplyr::matches("^macs_"))
    }

    #
    # Junctions overview
    #
    output$junctions_venn = renderImage({
      width = session$clientData$output_junctions_venn_width
      height = 960

      log("output$junctions_venn")
      venn_offtarget = list()
      venn_offtarget[[paste0("Bait chr.\n(", sum(tlx_df$tlx_is_bait_chromosome), ")")]] = tlx_df %>% dplyr::filter(tlx_is_bait_chromosome) %>% .$Qname
      venn_offtarget[[paste0("Bait region\n(", sum(tlx_df$tlx_is_bait_junction), ")")]] = tlx_df %>% dplyr::filter(tlx_is_bait_junction) %>% .$Qname
      venn_offtarget[[paste0("Repeats junct.\n(", sum(!is.na(tlx_df$tlx_repeatmasker_class)), ")")]] = tlx_df %>% dplyr::filter(!is.na(tlx_repeatmasker_class)) %>% .$Qname
      if("tlx_is_offtarget" %in% colnames(tlx_df)) {
        venn_offtarget[[paste0("Offtarget junct.\n(", sum(tlx_df$tlx_is_offtarget), ")")]] = tlx_df %>% dplyr::filter(tlx_is_offtarget) %>% .$Qname
      }

      session$userData$junctions_venn_svg = tempfile(fileext=".svg")
      log("Plot venn ", session$userData$junctions_venn_svg, " (w=", width, " h=", height, ")")
      svg(session$userData$junctions_venn_svg, width=width/72, height=height/72, pointsize=1)
      plot_venn(venn_offtarget, size=width)
      dev.off()

      session$userData$junctions_venn_svg = tempfile(fileext=".svg")
      log("Plot venn ", session$userData$junctions_venn_svg, " (w=", width, " h=", height, ")")

      svg(session$userData$junctions_venn_svg, width=width/72, height=height/72, pointsize=1)
      plot_venn(venn_offtarget, size=width)
      svg("Vivien/reports/junctions_count.svg", width=500/72, height=500/72, pointsize=1)
      tlx_df %>%
        dplyr::mutate(Subset=dplyr::case_when(
          tlx_is_bait_junction~"bait peak",
          tlx_is_offtarget~"offtarget peak",
          !is.na(tlx_repeatmasker_class)~"o/w repeat",
          T~"other junctions")) %>%
        ggplot() +
          geom_bar(aes(x=forcats::fct_infreq(tlx_sample), fill=Subset)) +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
          labs(x="", fill="Subset")
      dev.off()

      list(src=normalizePath(session$userData$junctions_venn_svg), contentType='image/svg+xml', width=width, height=height, alt="Junctions venn")
    }, deleteFile=F)

    #
    # Repeats summary plot
    #
    output$repeats_summary = renderImage({
      width = session$clientData$output_repeats_summary_width
      height = 960

      cat(file=stderr(), "repeats_summary")

      session$userData$repeats_summary_svg = tempfile(fileext=".svg")
      log("Plot repeats summary ", session$userData$repeats_summary_svg, " (w=", width, " h=", height, ")")

      svg(session$userData$repeats_summary_svg, width=width/72, height=height/72, pointsize=1)
      tlx_df.ggplot = dplyr::bind_rows(
        tlx_df %>% dplyr::mutate(filter="Other chromosomes") %>% dplyr::filter(!tlx_is_bait_chromosome),
        tlx_df %>% dplyr::mutate(filter="Bait chromosomes\n(excluding bait peak)") %>% dplyr::filter(tlx_is_bait_chromosome)) %>%
        dplyr::filter(!tlx_is_bait_junction) %>%
        tidyr::separate_rows(tlx_repeatmasker_class) %>%
        dplyr::filter(is.na(tlx_repeatmasker_class) | tlx_repeatmasker_class!="") %>%
        dplyr::mutate(tlx_repeatmasker_class=ifelse(is.na(tlx_repeatmasker_class), "No repeats", tlx_repeatmasker_class))
      palette_repeats = tlx_df.ggplot %>%
        dplyr::distinct(tlx_repeatmasker_class) %>%
        dplyr::mutate(color=ifelse(tlx_repeatmasker_class=="No repeats", "#CC0000", "#666666")) %>%
        dplyr::pull(color, tlx_repeatmasker_class)
      p = ggplot2::ggplot(tlx_df.ggplot) +
          ggplot2::geom_bar(ggplot2::aes(x=forcats::fct_infreq(tlx_repeatmasker_class), fill=tlx_repeatmasker_class), position="dodge") +
          ggplot2::coord_flip() +
          ggplot2::labs(x="RepeatMasker class", y="Junctions overlapping with repeats (x1000)") +
          ggplot2::facet_wrap(~filter, scales="free") +
          # ggplot2::scale_y_continuous(breaks=scale_breaks_sub1k) +
          ggplot2::scale_fill_manual(values=palette_repeats, guide="none") +
          ggplot2::theme_bw(base_size=24) +
          ggplot2::theme(axis.text.x=ggplot2::element_text(size=16))
      log(p)
      dev.off()

      list(src=normalizePath(session$userData$repeats_summary_svg), contentType='image/svg+xml', width=width, height=height, alt="Repeats summary")
    }, deleteFile=F)

    #
    # Homology plot
    #
    output$homology = renderImage({
      width = session$clientData$output_homology_width
      height = 960

      # @todo: develop single style for plots
      session$userData$homology_svg = tempfile(fileext=".svg")
      log(paste0("Plot homology plot ", session$userData$homology_svg, " (w=", size, " h=", size, ")"))
      svg(session$userData$homology_svg, width=size/72, height=size/72, pointsize=1)
      plot_homology(tlx_df)
      dev.off()

      log("Finished")

      list(src=normalizePath(session$userData$homology_svg), contentType='image/svg+xml', width=size, height=size, alt="TLX circos plot")
    }, deleteFile=F)

    #
    # Circos plot
    #
    output$circos = renderImage({
      size = pmin(session$clientData$output_circos_width, session$clientData$output_circos_width)

      links_df = macs_df %>% dplyr::inner_join(r$baits_df, by=c("macs_sample"="bait_sample"))
      if("macs_is_offtarget" %in% colnames(links_df)) {
        links_df$color = ifelse(links_df$macs_is_offtarget, "#F46D43FF", "#74ADD180")
      } else {
        links_df$color = "#74ADD180"
      }
      links_df = links_df %>%
        dplyr::select(chrom1=macs_chrom, start1=macs_start, end1=macs_end, chrom2=bait_chrom, start2=bait_start, end2=bait_end, color)

      # session$userData$circos_svg = "circos.svg"
      session$userData$circos_svg = tempfile(fileext=".svg")
      log("Plot circos ", session$userData$circos_svg, " (w=", size, " h=", size, ")")
      svg(session$userData$circos_svg, width=size/72, height=size/72, pointsize=1)
      plot_circos(data=tlx_df %>% dplyr::select(chrom=Rname, start=Rstart, end=Rend), cytoband_path=cytoband_path, links=links_df, cex=size/300)
      dev.off()

      log("Finished")

      list(src=normalizePath(session$userData$circos_svg), contentType='image/svg+xml', width=size, height=size, alt="TLX circos plot")
    }, deleteFile=F)
  })

}