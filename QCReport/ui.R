library(dplyr)
library(dplyr)
library(shiny)
library(dqshiny)
library(sortable)
library(shinycssloaders)


ui <- shiny::fluidPage(
  shinyjs::useShinyjs(),
  shiny::titlePanel("HTGTS postprocess"),
  shiny::HTML("<style>
  .tlx-group .shiny-file-input-progress { margin-bottom: 1px !important; }
  .tlx-group .shiny-file-input-progress { height: 15px !important; }
  .tlx-group .shiny-input-container { margin-bottom: 1px !important; }
  #paired_inputs { display: table-row; }
  #paired_inputs .shiny-input-container { display: table-cell; padding: 5px }
  #input_validation { color: #FF0000; font-weight: bold; background:#fffebd }
  </style>"),


  verbatimTextOutput("input_validation"),


  shiny::sidebarLayout(
    shiny::sidebarPanel(
      shiny::div(id="paired_inputs",
        shiny::checkboxInput("paired_controls", label="Paired controls", value=T),
        shiny::checkboxInput("paired_samples", label="Paired samples", value=T)
      ),
      #, fileInput("tlx1", label="", placeholder="No file selected")
      shiny::span(id="tlx_files"),
      shiny::actionButton("tlx_add", label="+"),
      shiny::actionButton("tlx_del", label="-"),

      shiny::br(),
      shiny::downloadLink("download_baits", "Download baits"),

      shiny::fileInput("offtargets", label="Offtargets file", placeholder = "No file selected"),
      shiny::selectInput("model", label="Model", selected="mm10", choices=c("mm9", "mm10", "hg19")),
      shiny::numericInput("qvalue", label="MACS2 qvalue", value=0.001),
      shiny::numericInput("pileup", label="MACS2 smallest pileup to call a peak", value=5),
      shiny::numericInput("extsize", label="extsize", value=10000),
      shiny::numericInput("maxgap", label="maxgap", value=0),

      shiny::wellPanel(
        shiny::selectInput("circos_chromosomes", label="Circos chromosomes", choices=c(), multiple=T),
        shiny::numericInput("circos_bw", label="Circos bin width", value=1e6)
      ),

      shiny::actionLink("advanced_hide", "advanced options"),
      shinyjs::hidden(
        shiny::wellPanel(id="advanced_panel",
          shiny::checkboxInput("exclude_repeats", label="Exclude repeats", value=F),
          shiny::selectInput("exttype", label="exttype", selected="along", choices=c("along", "symmetrical", "none")),
          shiny::fluidRow(
            shiny::column(2, shiny::checkboxInput("exclude_bait_region", label="Exclude bait region", value=T)),
            shiny::column(10, shiny::numericInput("bait_region", label="", value=1500000))
          ),
          shiny::numericInput("breaksite_size", label="breaksite size", value=19),
          shiny::numericInput("slocal", label="slocal", value=10000),
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
          shiny::plotOutput("test", height="auto"),
          shiny::plotOutput("junctions_venn", height="auto"),
          shiny::plotOutput("junctions_count", height="auto")
        ),
        shiny::tabPanel(
          "Repeats",
          shiny::plotOutput("repeats_summary", height="auto")
        ),
        shiny::tabPanel(
          "Homology profile",
          shiny::plotOutput("homology", height="auto")
        ),
        shiny::tabPanel(
          "Overview",
        shiny::downloadLink("download_macs2", "Download MACS2"),
          shiny::plotOutput("circos", height="auto") %>% shinycssloaders::withSpinner(color="#0dc5c1")
        ),
        shiny::tabPanel(
          "Compare hits",
          shiny::plotOutput("compare_breaks", height="auto") %>% shinycssloaders::withSpinner(color="#0dc5c1"),
          shiny::plotOutput("compare_pileup", height="auto") %>% shinycssloaders::withSpinner(color="#0dc5c1")
        )
      )
    )
  ),

  shiny::actionLink("logs_hide", "show logs"),
  textOutput("logs")
)