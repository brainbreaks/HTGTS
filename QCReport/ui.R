library(dplyr)
library(dplyr)
library(shiny)
library(dqshiny)
library(shinycssloaders)

ui <- shiny::fluidPage(
  shinyjs::useShinyjs(),
  shiny::titlePanel("title panel"),

  shiny::sidebarLayout(
    shiny::sidebarPanel(
      #, fileInput("tlx1", label="", placeholder="No file selected")
      shiny::span(id="tlx_files"),
      shiny::actionButton("tlx_add", label="+"),
      shiny::actionButton("tlx_del", label="-"),

      shiny::fileInput("offtargets", label="Offtargets file", placeholder = "No file selected"),
      shiny::selectInput("model", label="Model", selected="mm10", choices=c("mm10", "hg19")),
      shiny::numericInput("qvalue", label="MACS2 qvalue", value=0.001),

      shiny::actionLink("advanced_hide", "advanced options"),
      shinyjs::hidden(
        shiny::wellPanel(id="advanced_panel",
          shiny::checkboxInput("exclude_repeats", label="Exclude repeats", value=T),
          shiny::fluidRow(
            shiny::column(2, shiny::checkboxInput("exclude_bait_region", label="Exclude bait region", value=T)),
            shiny::column(10, shiny::numericInput("bait_region", label="", value=500000))
          ),
          shiny::numericInput("breaksite_size", label="breaksite_size size", value=19),
          shiny::numericInput("junsize", label="Junction size", value=300),
          shiny::numericInput("extsize", label="extsize", value=2000),
          shiny::numericInput("slocal", label="slocal", value=1000),
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
          shiny::plotOutput("circos") %>% shinycssloaders::withSpinner(color="#0dc5c1")
        )
      )
    )
  )
)