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
source("utils.R")
source("tlx.R")
source("graphics.R")

download_link = function(id, data, file) {
  if(!is.null(data) && nrow(data)>0) {
    readr::write_tsv(data, file=file, na="")
  }
}

server <- function(input, output, session) {
  # @todo: properly load default
  r <- shiny::reactiveValues(
    hidden_options=T,
    baits_df=NULL,
    tlx_df=tlx_blank(),
    offtargets_df=NULL,
    tlx_num=0,
    repeatmasker_df=NULL,
    observers=list()
  )


  observeEvent(input$advanced_hide, {
    log(input$advanced_hide)
    shinyjs::toggle(id="advanced_panel")
  })

  observeEvent(r$tlx_df, {
    log("r$tlx_df")
    r$baits_df = tlx_identify_baits(r$tlx_df, breaksite_size=input$breaksite_size)
    shinyjs::toggle(id="download_baits", condition=nrow(r$tlx_df)>0)
  }, ignoreNULL=F)

  output$download_baits <- downloadHandler(filename="baits.tsv",
    content = function(file) {
      log("output$download_baits ")
      download_link(id="download_baits", data=r$baits_df, file=file)
    }
  )

  observeEvent(input$genome, {
    log("input$genome")
    r$repeatmasker_df = repeatmasker_read(file.path(genomes_path, input$genome, "annotation/ucsc_repeatmasker.tsv"))
  })

  observeEvent(input$tlx_add, {
    log("input$tlx_add")

    group_i = r$tlx_num+1
    input_id = paste0("tlx_input", group_i)
    control_id = paste0("tlx_control", group_i)

    shiny::insertUI("#tlx_files", where="beforeEnd", immediate=T, ui=shiny::div(class="tlx_group",
      fileInput(input_id, label="Input", placeholder="No file selected", multiple=T),
      fileInput(control_id, label="Control", placeholder="No file selected", multiple=T))
    )

    session$userData[[paste0("tlx", group_i, "_observer")]] = observeEvent(c(input[[input_id]], input[[control_id]]), {
      log(paste0("tlx", group_i, "_observer"))

      group_id = paste("group", group_i)
      tlx_df = isolate(r$tlx_df %>% dplyr::filter(tlx_group!=group_id))
      for(field_id2 in c(input_id, control_id)) {
        if(is.null(input[[field_id2]])) next
        if(length(input[[field_id2]][,1])==0) next

        for(i in 1:length(input[[field_id2]][,1])) {
          tlx_df = rbind(tlx_df, tlx_read(input[[field_id2]][[i, "datapath"]], sample=basename(input[[field_id2]][[i, "name"]]), group=group_id, control=grepl("control", field_id2)))
        }
      }

      if(nrow(tlx_df)) {
        tlx_df = tlx_remove_rand_chromosomes(tlx_df)
        tlx_df = tlx_mark_bait_chromosome(tlx_df)
        tlx_df = tlx_mark_bait_junctions(tlx_df, isolate(input$bait_region))
        tlx_df$tlx_is_offtarget = F
      }
      r$tlx_df = tlx_df
    }, ignoreInit=T)

    r$tlx_num = r$tlx_num + 1
    shinyjs::toggle("tlx_del", condition=r$tlx_num>0)
  }, ignoreNULL=F)


  observeEvent(input$tlx_del, {
    log("input$tlx_del")
    group_i = r$tlx_num
    input_id = paste0("tlx", r$tlx_num)
    control_id = paste0("tlx_control", group_i)

    for(field in c(input_id, control_id)) {
      if(is.null(input[[field]])) next
      if(length(input[[field]][,1])==0) next

      for(i in 1:length(input[[field]][,1])) {
        r$tlx_df = r$tlx_df %>% dplyr::filter(tlx_path!=input[[field]][[i, "datapath"]])
      }
    }

    shiny::removeUI(selector=".tlx_group:last-of-type")
    session$userData[[paste0("tlx", group_i, "_observer")]]$destroy()

    r$tlx_num = r$tlx_num - 1
    shinyjs::toggle("tlx_del", condition=r$tlx_num>0)
  }, ignoreInit=T)


  observeEvent(input$model, {
    log("input$model")
    r$repeatmasker_df = repeatmasker_read(file.path(genomes_path, input$model, "annotation/ucsc_repeatmasker.tsv"), columns=c("repeatmasker_chrom", "repeatmasker_start", "repeatmasker_end", "repeatmasker_class"))
    log("Repeatmasker file loaded!")
  })


  observeEvent(input$offtargets, {
    log("input$offtargets")
    log(basename(input$offtargets$name))
    r$offtargets_df = offtargets_read(input$offtargets$datapath)
  }, ignoreInit=T)


  observeEvent(input$calculate, {
    req(!is.null(isolate(r$tlx_df)))
    log("input$calculate")

    exclude_repeats = shiny::isolate(input$exclude_repeats)
    exclude_bait_region = shiny::isolate(input$exclude_bait_region)
    extsize = shiny::isolate(input$extsize)
    qvalue = shiny::isolate(input$qvalue)
    pileup = shiny::isolate(input$pileup)
    slocal = shiny::isolate(input$slocal)
    llocal = shiny::isolate(input$llocal)
    model = shiny::isolate(input$model)

    log_input(input)

    tlx_df = tlx_mark_repeats(shiny::isolate(r$tlx_df), r$repeatmasker_df)



    #
    # Vivien's
    #
    # setwd("/home/s215v/Workspace/HTGTS/QCReport")
    # r = list(); width = 500; height=500; junsize = 300; exclude_repeats = F; exclude_bait_region = T; bait_region = 500000; extsize = 2000; qvalue = 0.001; slocal = 1000; llocal = 10000000; model = "mm10"; genomes_path = "/home/s215v/Workspace/HTGTS/genomes"
    # session = list(userData=list(
    #   repeats_summary_svg="Vivien/reports/repeats_summary.svg",
    #   junctions_venn_svg = "Vivien/reports/junctions_venn.svg",
    #   homology_svg = "Vivien/reports/homology_profile.svg",
    #   circos_svg = "Vivien/reports/circos.svg",
    #   junctions_count_svg="Vivien/reports/junctions_count.svg"
    # ))
    # r$repeatmasker_df = repeatmasker_read("genomes/mm10/annotation/ucsc_repeatmasker.tsv", columns=c("repeatmasker_chrom", "repeatmasker_start", "repeatmasker_end", "repeatmasker_class"))
    # samples_df = data.frame(path=list.files("Vivien", pattern="*.tlx", full.names=T)) %>%
    #   dplyr::mutate(sample=gsub("_B400_0", " / ", gsub("_result.tlx", "", basename(path))), group="group1", control=grepl("VI05", path))
    # # r$offtargets_df = offtargets_read("Vivien/offtargets.bed")
    # tlx_df = tlx_read_many(samples_df)
    # tlx_df = tlx_remove_rand_chromosomes(tlx_df)
    # tlx_df = tlx_mark_bait_chromosome(tlx_df)
    # tlx_df = tlx_mark_bait_junctions(tlx_df, bait_region)
    # tlx_df = tlx_mark_repeats(tlx_df, r$repeatmasker_df)
    # tlx_df$tlx_is_offtarget = F
    # r$offtargets_df = NULL
    # r$baits_df = tlx_identify_baits(tlx_df)
    #
    # End Vivien's
    #

    genome_path = file.path(genomes_path, model, paste0(model, ".fa"))
    cytoband_path = file.path(genomes_path, model, "annotation/cytoBand.txt")

    offtarget2bait_df = NULL
    if(!is.null(r$offtargets_df)) {
      offtarget2bait_df = join_offtarget2bait(offtargets_df=r$offtargets_df, baits_df=r$baits_df, genome_path=genome_path)
      tlx_df = tlx_mark_offtargets(tlx_df, offtarget2bait_df)
    }

    macs_df = tlx_macs2(tlx_df, extsize=extsize, qvalue=qvalue, pileup=pileup, slocal=slocal, llocal=llocal, exclude_bait_region=exclude_bait_region, exclude_repeats=exclude_repeats)
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
      log("output$junctions_venn ")
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

      list(src=normalizePath(session$userData$junctions_venn_svg), contentType='image/svg+xml', width=width, height=height, alt="Junctions venn")
    }, deleteFile=F)



    output$junctions_count = renderImage({
      log("output$junctions_count")
      width = session$clientData$output_junctions_count_width
      height = 960

      session$userData$junctions_count_svg = tempfile(fileext=".svg")
      log("Plot junctions count ", session$userData$junctions_count_svg, " (w=", width, " h=", height, ")")
      svg(session$userData$junctions_count_svg, width=width/72, height=height/72, pointsize=1)
      p = tlx_df %>%
        dplyr::mutate(Subset=dplyr::case_when(
          tlx_is_bait_junction~"bait peak",
          tlx_is_offtarget~"offtarget peak",
          !is.na(tlx_repeatmasker_class)~"o/w repeat",
          T~"other junctions")) %>%
          ggplot() +
            geom_bar(aes(x=forcats::fct_infreq(tlx_sample), fill=Subset)) +
            theme_classic(base_size=18) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.key.size = unit(20, 'pt')) +
            labs(x="", fill="Subset")
        print(p)
      dev.off()

      list(src=normalizePath(session$userData$junctions_count_svg), contentType='image/svg+xml', width=width, height=height, alt="Junctions count")
    }, deleteFile=F)


    #
    # Repeats summary plot
    #
    output$repeats_summary = renderImage({
      log("output$repeats_summary")
      width = session$clientData$output_repeats_summary_width
      height = 960

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
          ggplot2::theme_classic(base_size=18)+
          ggplot2::theme(strip.background=element_blank()) +
          ggplot2::theme(axis.text.x=ggplot2::element_text(size=16))
      print(p)
      dev.off()

      list(src=normalizePath(session$userData$repeats_summary_svg), contentType='image/svg+xml', width=width, height=height, alt="Repeats summary")
    }, deleteFile=F)

    #
    # Homology plot
    #
    output$homology = renderImage({
      log("output$homology ")
      width = session$clientData$output_homology_width
      height = 960

      # @todo: develop single style for plots
      session$userData$homology_svg = tempfile(fileext=".svg")
      log(paste0("Plot homology plot ", session$userData$homology_svg, " (w=", width, " h=", height, ")"))
      svg(session$userData$homology_svg, width=width/72, height=height/72, pointsize=1)
      plot_homology(tlx_df)
      dev.off()

      log("Finished")

      list(src=normalizePath(session$userData$homology_svg), contentType='image/svg+xml', width=width, height=height, alt="TLX circos plot")
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
  }, ignoreInit=T)

}