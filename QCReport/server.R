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

download_link = function(data, file) {
  if(!is.null(data) && nrow(data)>0) {
    readr::write_tsv(data, file=file, na="")
  }
}

event = function(expr) {
  withCallingHandlers(expr, message=function(m) {
    shinyjs::html(id="logs", html=gsub("\n", "<br />", m$message), add=T)
  })
}

server <- function(input, output, session) {
  # @todo: properly load default
  r <- shiny::reactiveValues(
    hidden_options=T,
    baits_df=NULL,
    tlx_df=tlx_blank(),
    macs_df=macs_blank(),
    offtargets_df=NULL,
    tlx_num=0,
    repeatmasker_df=NULL,
    observers=list()
  )


  observeEvent(input$advanced_hide, event({
    log(input$advanced_hide)
    shinyjs::toggle(id="advanced_panel")
  }))

  observeEvent(r$tlx_df, event({
    log("r$tlx_df")
    r$baits_df = tlx_identify_baits(r$tlx_df, breaksite_size=input$breaksite_size)
    shinyjs::toggle(id="download_baits", condition=nrow(r$tlx_df)>0)
  }), ignoreNULL=F)

  observeEvent(input$logs_hide, event({
    shinyjs::toggle(id="logs")
  }), ignoreNULL=F)

  output$download_baits <- downloadHandler(filename="baits.tsv",
    content = function(file) {
      event({
        log("output$download_baits ")
        download_link(data=r$baits_df, file=file)
      })
    }
  )


  output$download_macs2 = downloadHandler(filename="macs2.tsv",
    content = function(file) {
      event({
        log("output$download_macs2")
        download_link(data=r$macs_df, file=file)
      })
    }
  )

  observeEvent(input$genome, event({
    log("input$genome")
    r$repeatmasker_df = repeatmasker_read(file.path(genomes_path, input$genome, "annotation/ucsc_repeatmasker.tsv"))
  }))

  observeEvent(input$tlx_add, event({
    log("input$tlx_add")

    group_i = r$tlx_num+1
    input_id = paste0("tlx_input", group_i)
    control_id = paste0("tlx_control", group_i)
    group_id = paste("group", group_i)

    shiny::insertUI("#tlx_files", where="beforeEnd", immediate=T, ui=shiny::div(class="tlx-group",
      fileInput(input_id, label="Input", placeholder="No file selected", multiple=T),
      shiny::downloadLink(paste0("download_pileup_", input_id), "wig"),
      shiny::downloadLink(paste0("download_bed_", input_id), "bed"),
      fileInput(control_id, label="Control", placeholder="No file selected", multiple=T),
      shiny::downloadLink(paste0("download_pileup_", control_id), "wig"),
      shiny::downloadLink(paste0("download_bed_", control_id), "bed"),
    ))

    session$userData[[paste0("tlx", group_i, "_observer")]] = observeEvent(c(input[[input_id]], input[[control_id]]), event({
      log(paste0("tlx", group_i, "_observer"))

      tlx_df = isolate(r$tlx_df %>% dplyr::filter(tlx_group!=group_id))
      for(field_id2 in c(input_id, control_id)) {
        if(is.null(input[[field_id2]])) next
        if(length(input[[field_id2]][,1])==0) next

        for(i in 1:length(input[[field_id2]][,1])) {
          tlx_df.i = tlx_read(input[[field_id2]][[i, "datapath"]], sample=basename(input[[field_id2]][[i, "name"]]), group=group_id, control=grepl("control", field_id2)) %>%
            dplyr::mutate(tlx_sample=gsub("(_result)?\\.tlx$", "", tlx_sample))
          tlx_df = dplyr::bind_rows(tlx_df, tlx_df.i)
        }
      }

      if(nrow(tlx_df)) {
        tlx_df = tlx_remove_rand_chromosomes(tlx_df)
        tlx_df = tlx_mark_bait_chromosome(tlx_df)
        tlx_df = tlx_mark_bait_junctions(tlx_df, isolate(input$bait_region))
        tlx_df$tlx_is_offtarget = F
      }

      r$tlx_df = tlx_df

      for(field_id in c(input_id, control_id)) {
        shinyjs::show(paste0("download_pileup_", field_id))
        shinyjs::show(paste0("download_bed_", field_id))
      }
    }), ignoreInit=T)

    r$tlx_num = r$tlx_num + 1
    shinyjs::toggle("tlx_del", condition=r$tlx_num>0)

    for(field_id in c(input_id, control_id)) {
      shinyjs::hide(paste0("download_pileup_", field_id))
      shinyjs::hide(paste0("download_bed_", field_id))
    }

    # @ todo: this should be possible to put into loop (be carefull about context!!!)
    output[[paste0("download_pileup_", input_id)]] = downloadHandler(filename=paste0(input_id, ".wig"),
      content=function(file) {
        event({
          log("input$download_pileup_", input_id)
          log(group_id)
          tlx_df.f = r$tlx_df %>% dplyr::filter(!tlx_control & tlx_group==group_id)
          tlx_write_wig(tlx_df.f, file=file, extsize=input$extsize)
        })
    })

    output[[paste0("download_bed_", input_id)]] = downloadHandler(filename=paste0(input_id, ".bed"),
      content=function(file) {
        event({
          log("input$download_bed_", input_id)
          tlx_df.f = r$tlx_df %>% dplyr::filter(!tlx_control & tlx_group==group_id)
          tlx_write_bed(tlx_df.f, file=file)
        })
      })

    output[[paste0("download_pileup_", control_id)]] = downloadHandler(filename=paste0(control_id, ".wig"),
      content=function(file) {
        event({
          log("input$download_pileup_", control_id)
          log(group_id)
          tlx_df.f = r$tlx_df %>% dplyr::filter(tlx_control & tlx_group==group_id)
          tlx_write_wig(tlx_df.f, file=file, extsize=input$extsize)
        })
    })

    output[[paste0("download_bed_", control_id)]] = downloadHandler(filename=paste0(control_id, ".bed"),
      content=function(file) {
        event({
          log("input$download_bed_", control_id)
          log(group_id)
          tlx_df.f = r$tlx_df %>% dplyr::filter(tlx_control & tlx_group==group_id)
          tlx_write_bed(tlx_df.f, file=file)
        })
      })

  }), ignoreNULL=F)


  observeEvent(input$tlx_del, event({
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

    shiny::removeUI(selector=".tlx-group:last-of-type")
    session$userData[[paste0("tlx", group_i, "_observer")]]$destroy()

    r$tlx_num = r$tlx_num - 1
    shinyjs::toggle("tlx_del", condition=r$tlx_num>0)
  }), ignoreInit=T)


  observeEvent(input$model, event({
    log("input$model")
    r$repeatmasker_df = repeatmasker_read(file.path(genomes_path, input$model, "annotation/ucsc_repeatmasker.tsv"), columns=c("repeatmasker_chrom", "repeatmasker_start", "repeatmasker_end", "repeatmasker_class"))

    cytoband_path = file.path(genomes_path, input$model, "annotation/cytoBand.txt")
    cytoband = circlize::read.cytoband(cytoband_path)
    shiny::updateSelectInput(inputId="circos_chromosomes", choices=cytoband$chromosome, selected=c())

    log("Repeatmasker file loaded!")
  }))


  output$input_validation <- renderText({
    log(input$extsize, " <= ", input$slocal)
    shiny::validate(need(!any(r$tlx_df$tlx_control) || input$extsize <= input$slocal, 'slocal should be less or equal to extsize when calling peaks with background'))
  })


  observeEvent(input$offtargets, event({
    log("input$offtargets")
    log(basename(input$offtargets$name))
    r$offtargets_df = offtargets_read(input$offtargets$datapath)
  }), ignoreInit=T)

  observeEvent(input$calculate, event({
    log("input$calculate")

    req(!is.null(r$tlx_df) && nrow(r$tlx_df))

    exclude_repeats = shiny::isolate(input$exclude_repeats)
    exclude_bait_region = shiny::isolate(input$exclude_bait_region)
    extsize = shiny::isolate(input$extsize)
    qvalue = shiny::isolate(input$qvalue)
    pileup = shiny::isolate(input$pileup)
    slocal = shiny::isolate(input$slocal)
    llocal = shiny::isolate(input$llocal)
    model = shiny::isolate(input$model)
    circos_bw = shiny::isolate(input$circos_bw)

    genome_path = file.path(genomes_path, model, paste0(model, ".fa"))
    cytoband_path = file.path(genomes_path, model, "annotation/cytoBand.txt")

    circos_chromosomes = shiny::isolate(input$circos_chromosomes)
    if(length(circos_chromosomes)==0) {
      cytoband = circlize::read.cytoband(cytoband_path)
      circos_chromosomes = cytoband$chromosome
    }

    log_input(input)

    tlx_df = tlx_mark_repeats(shiny::isolate(r$tlx_df), r$repeatmasker_df)
    baits_df = shiny::isolate(r$baits_df)

    #
    # Vivien's
    #
    # setwd("/home/s215v/Workspace/HTGTS/QCReport")
    # r = list(); width = 500; height=500; circos_bw=50000; circos_chromosomes="chr6"; pileup=5; exclude_repeats = F; exclude_bait_region = T; bait_region = 500000; extsize = 2000; qvalue = 0.001; slocal = 2000; llocal = 10000000; model = "mm10"; genomes_path = "/home/s215v/Workspace/HTGTS/genomes"
    # session = list(userData=list(
    #   repeats_summary_svg="Vivien/reports/repeats_summary.svg",
    #   junctions_venn_svg = "Vivien/reports/junctions_venn.svg",
    #   homology_svg = "Vivien/reports/homology_profile.svg",
    #   circos_svg = "Vivien/reports/circos.svg",
    #   junctions_count_svg="Vivien/reports/junctions_count.svg",
    #   compare_pileup_svg="Vivien/reports/compare_pileup.svg"
    # ))
    # r$repeatmasker_df = repeatmasker_read("genomes/mm10/annotation/ucsc_repeatmasker.tsv", columns=c("repeatmasker_chrom", "repeatmasker_start", "repeatmasker_end", "repeatmasker_class"))
    #
    # cytoband_path = file.path(genomes_path, model, "annotation/cytoBand.txt")
    # cytoband = circlize::read.cytoband(cytoband_path)
    # #circos_chromosomes = cytoband$chromosome
    #
    # samples_df = rbind(readr::read_tsv("Vivien/B400_011_metadata_2.txt"), readr::read_tsv("Vivien/B400_012_metadata_complete.txt")) %>%
    #   dplyr::filter(grepl("55|56|57|58", Library)) %>%
    #   dplyr::mutate(path=str_glue("Vivien/{lib}_{seq}_result.tlx", lib=Library, seq=Sequencing),
    #                 group=gsub("^[^0-9]+[0-9/]+ ?", "", gsub(".*(One.*)", "\\1", gsub(" - (DMSO|APH)", "", Description)), perl=T),
    #                 sample=paste0(Library, ifelse(grepl("promoter", Description), "-prom", "+prom")),
    #                 control=!grepl("APH", Description, perl=T)) %>%
    #   dplyr::select(path, sample, group, control, Description) %>%
    #   data.frame()
    # # r$offtargets_df = offtargets_read("Vivien/offtargets.bed")
    # tlx_df = tlx_read_many(samples_df)
    # tlx_df = tlx_remove_rand_chromosomes(tlx_df)
    # tlx_df = tlx_mark_bait_chromosome(tlx_df)
    # tlx_df = tlx_mark_bait_junctions(tlx_df, bait_region)
    # tlx_df = tlx_mark_repeats(tlx_df, r$repeatmasker_df)
    # tlx_df$tlx_is_offtarget = F
    # baits_df = tlx_identify_baits(tlx_df)
    #
    # End Vivien's
    #

    # for(p in unique(tlx_df$tlx_path)) {
    #   tlx_df.p = tlx_df %>%
    #     dplyr::filter(tlx_path==p) %>%
    #     dplyr::mutate(strand=ifelse(Strand=="-1", "-", "+"), start=Junction-200, end=Junction+200, score=0)
    #
    #   tlx_df.p %>%
    #     dplyr::select(Rname, start, end, Qname, score, strand) %>%
    #     readr::write_tsv(file=paste0("Vivien/bed/", tlx_df.p$tlx_sample[1], ".bed"), na="", col_names=F)
    # }


    groups_n = length(unique(tlx_df$tlx_group))
    if(groups_n>1) {
      groups_grid = matrix(1:(ceiling(groups_n/2)*2), ncol=2, byrow=T)
    } else {
      groups_grid = matrix(1)
    }

    offtarget2bait_df = NULL
    if(!is.null(r$offtargets_df)) {
      offtarget2bait_df = join_offtarget2bait(offtargets_df=r$offtargets_df, baits_df=baits_df, genome_path=genome_path)
      tlx_df = tlx_mark_offtargets(tlx_df, offtarget2bait_df)
    }

    r$macs_df = tlx_macs2(tlx_df, extsize=extsize, qvalue=qvalue, pileup=pileup, slocal=slocal, llocal=llocal, exclude_bait_region=exclude_bait_region, exclude_repeats=exclude_repeats)
    if(!is.null(offtarget2bait_df)) {
      r$macs_df = r$macs_df %>%
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


      session$userData$junctions_venn_svg = tempfile(fileext=".svg")
      log("Plot venn ", session$userData$junctions_venn_svg, " (w=", width, " h=", height, ")")
      svg(session$userData$junctions_venn_svg, width=width/72, height=height/72, pointsize=1)
      pushViewport(plotViewport(layout=grid.layout(nrow=nrow(groups_grid), ncol=ncol(groups_grid))))
      for(gr in unique(tlx_df$tlx_group)) {
        gr_i = which(gr==unique(tlx_df$tlx_group))
        gr_pos = which(groups_grid==gr_i, arr.ind=T)

        tlx_df.gr = tlx_df %>% dplyr::filter(tlx_group==gr)
        venn_offtarget = list()
        venn_offtarget[[paste0("Bait chr.\n(", sum(tlx_df.gr$tlx_is_bait_chromosome), ")")]] = tlx_df.gr %>% dplyr::filter(tlx_is_bait_chromosome) %>% .$Qname
        venn_offtarget[[paste0("Bait region\n(", sum(tlx_df.gr$tlx_is_bait_junction), ")")]] = tlx_df.gr %>% dplyr::filter(tlx_is_bait_junction) %>% .$Qname
        venn_offtarget[[paste0("Repeats junct.\n(", sum(!is.na(tlx_df.gr$tlx_repeatmasker_class)), ")")]] = tlx_df.gr %>% dplyr::filter(!is.na(tlx_repeatmasker_class)) %>% .$Qname
        if("tlx_is_offtarget" %in% colnames(tlx_df.gr)) {
          venn_offtarget[[paste0("Offtarget junct.\n(", sum(tlx_df.gr$tlx_is_offtarget), ")")]] = tlx_df.gr %>% dplyr::filter(tlx_is_offtarget) %>% .$Qname
        }

        pushViewport(plotViewport(layout.pos.row=gr_pos[1], layout.pos.col=gr_pos[2]))
        main = paste0(gr, " (Total junctions ", nrow(tlx_df.gr), ")")
        plot_venn(venn_offtarget, main=gr, size=width/ncol(groups_grid))
        popViewport()
      }
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
      tlx_df.sum = tlx_df %>%
        dplyr::mutate(Subset=dplyr::case_when(
          tlx_is_bait_junction~"bait peak",
          tlx_is_offtarget~"offtarget peak",
          !is.na(tlx_repeatmasker_class)~"o/w repeat",
          T~"other junctions")) %>%
        dplyr::group_by(tlx_group, tlx_sample) %>%
        dplyr::mutate(total_n=n()) %>%
        dplyr::group_by(tlx_group, tlx_sample, Subset, total_n) %>%
        dplyr::summarize(subset_n=n())

        gridExtra::grid.arrange(
          ggplot(tlx_df.sum) +
            geom_bar(aes(x=tlx_sample, y=subset_n, fill=Subset), stat="identity", show.legend=F) +
            theme_brain() +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.key.size = unit(20, 'pt')) +
            labs(x="", y="Junctions (absolute)", fill="Subset") +
            facet_wrap(~tlx_group, scales="free_x"),
          ggplot(tlx_df.sum) +
            geom_bar(aes(x=tlx_sample, y=subset_n/total_n, fill=Subset), stat="identity") +
            scale_y_continuous(labels = scales::percent) +
            theme_brain() +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.key.size = unit(20, 'pt')) +
            theme(legend.position="bottom") +
            guides(fill=guide_legend(nrow=2, byrow=T)) +
            labs(x="", y="Junctions (relative)", fill="Subset") +
            facet_wrap(~tlx_group, scales="free_x"),
          ncol=1, heights=c(0.45, 0.55)
        )
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
          ggplot2::labs(x="", y="Junctions overlapping with repeats (x1000)") +
          ggplot2::facet_grid(tlx_group~filter, scales="free") +
          # ggplot2::scale_y_continuous(breaks=scale_breaks_sub1k) +
          ggplot2::scale_fill_manual(values=palette_repeats, guide="none") +
          theme_brain()
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

      list(src=normalizePath(session$userData$homology_svg), contentType='image/svg+xml', width=width, height=height, alt="homology plot")
    }, deleteFile=F)

    #
    # Circos plot
    #
    output$circos = renderImage({
      width = session$clientData$output_circos_width
      height = (width + 100) * groups_n

      session$userData$circos_svg = tempfile(fileext=".svg")
      log("Plot circos ", session$userData$circos_svg, " (w=", width, " h=", height, ")")
      svg(session$userData$circos_svg, width=width/72, height=height/72, pointsize=1)

      layout(matrix(1:groups_n, groups_n, 1))
      for(gr in unique(tlx_df$tlx_group)) {
        tlx_df.gr = tlx_df %>% dplyr::filter(tlx_group==gr)
        macs_df.gr = r$macs_df %>% dplyr::filter(macs_group==gr)
        baits_df.gr = baits_df %>% dplyr::distinct(bait_group, .keep_all=T)
        hits_df.gr = macs_df.gr %>% dplyr::select(chrom=macs_chrom, start=macs_start, end=macs_end)

        links_df.gr = macs_df.gr %>% dplyr::inner_join(baits_df.gr, by=c("macs_group"="bait_group"))
        if("macs_is_offtarget" %in% colnames(links_df.gr)) {
          links_df.gr$color = ifelse(links_df.gr$macs_is_offtarget, "#F46D43FF", "#74ADD180")
        } else {
          links_df.gr$color = "#74ADD180"
        }
        links_df.gr = links_df.gr %>%
          dplyr::select(chrom1=macs_chrom, start1=macs_start, end1=macs_end, chrom2=bait_chrom, start2=bait_start, end2=bait_end, color)

        plot_circos(
          input=tlx_df.gr %>% dplyr::filter(!tlx_control) %>% dplyr::select(chrom=Rname, start=Rstart, end=Rend),
          control=tlx_df.gr %>% dplyr::filter(tlx_control) %>% dplyr::select(chrom=Rname, start=Rstart, end=Rend),
          title=gr,
          annotations=list("hits"=hits_df.gr),
          cytoband_path=cytoband_path,
          chromosomes=circos_chromosomes,
          links=links_df.gr,
          circos_bw=circos_bw,
          cex=width/300)
      }
      dev.off()

      log("Finished")

      list(src=normalizePath(session$userData$circos_svg), contentType='image/svg+xml', width=width, height=height, alt="TLX circos plot")
    }, deleteFile=F)


    #
    # Compare MACS hits across groups
    #
    output$compare_pileup = renderImage({
      log("output$compare_pileup")

      width = session$clientData$output_compare_pileup_width
      height = 960

      session$userData$compare_pileup_svg = tempfile(fileext=".svg")
      log(paste0("Plot compare plot ", session$userData$compare_pileup_svg, " (w=", width, " h=", height, ")"))
      svg(session$userData$compare_pileup_svg, width=width/72, height=height/72, pointsize=1)
      print(plot_macs2_pileups(tlx_df, r$macs_df))
      dev.off()
      list(src=normalizePath(session$userData$compare_pileup_svg), contentType='image/svg+xml', width=width, height=height, alt="Compare pileups for MACS2 hits")
    }, deleteFile=F)
  }))
}