library(dplyr)
library(dplyr)
library(shiny)
library(dqshiny)
library(shinycssloaders)


ui <- shiny::fluidPage(
  shinyjs::useShinyjs(),
  shiny::titlePanel("HTGTS postprocess"),
  shiny::HTML("<style> .tlx-group .shiny-file-input-progress { margin-bottom: 1px !important; } </style>"),
  shiny::HTML("<style> .tlx-group .shiny-file-input-progress { height: 15px !important; } </style>"),
  shiny::HTML("<style> .tlx-group .shiny-input-container { margin-bottom: 1px !important; } </style>"),



  shiny::sidebarLayout(
    shiny::sidebarPanel(
      #, fileInput("tlx1", label="", placeholder="No file selected")
      shiny::span(id="tlx_files"),
      shiny::actionButton("tlx_add", label="+"),
      shiny::actionButton("tlx_del", label="-"),

      shiny::br(),
      shiny::downloadLink("download_baits", "Download baits"),

      shiny::fileInput("offtargets", label="Offtargets file", placeholder = "No file selected"),
      shiny::selectInput("model", label="Model", selected="mm10", choices=c("mm10", "hg19")),
      shiny::numericInput("qvalue", label="MACS2 qvalue", value=0.001),
      shiny::numericInput("pileup", label="MACS2 smallest pileup to call a peak", value=5),
      shiny::selectInput("circos_chromosomes", label="Chromosomes", choices=c(), multiple=T),


      shiny::actionLink("advanced_hide", "advanced options"),
      shinyjs::hidden(
        shiny::wellPanel(id="advanced_panel",
          shiny::checkboxInput("exclude_repeats", label="Exclude repeats", value=F),
          shiny::fluidRow(
            shiny::column(2, shiny::checkboxInput("exclude_bait_region", label="Exclude bait region", value=T)),
            shiny::column(10, shiny::numericInput("bait_region", label="", value=500000))
          ),
          shiny::numericInput("breaksite_size", label="breaksite size", value=19),
          shiny::numericInput("extsize", label="extsize", value=300),
          shiny::numericInput("slocal", label="slocal", value=50000),
          shiny::numericInput("llocal", label="llocal", value=10000000)
        )
      ),
      shiny::HTML("<br />"),
      shiny::actionButton("calculate", label="Calculate")
    ),

    shiny::mainPanel(
      type="tabs",
      shiny::tabsetPanel(
        shiny::tabPanel(
          "Junctions",
          shiny::plotOutput("junctions_venn", height = "auto"),
          shiny::plotOutput("junctions_count", height = "auto")
        ),
        shiny::tabPanel(
          "Repeats",
          shiny::plotOutput("repeats_summary")
        ),
        shiny::tabPanel(
          "Homology profile",
          shiny::plotOutput("homology")
        ),
        shiny::tabPanel(
          "Overview",
        shiny::downloadLink("download_macs2", "Download MACS2"),
          shiny::plotOutput("circos") %>% shinycssloaders::withSpinner(color="#0dc5c1")
        ),
        shiny::tabPanel(
          "Compare hits",
          shiny::plotOutput("compare_pileup") %>% shinycssloaders::withSpinner(color="#0dc5c1")
        )
      )
    )
  )
)