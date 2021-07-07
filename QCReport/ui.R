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
      shiny::fileInput("tlx", label="TLX file", placeholder = "No file selected"),
      shiny::fileInput("offtargets", label="Offtargets file", placeholder = "No file selected"),
      shiny::selectInput("model", label="Model", choices=c("hg19", "mm10")),
      shiny::numericInput("qvalue", label="MACS2 qvalue", value=0.001),

      shiny::actionLink("advanced_hide", "advanced options"),
      shiny::wellPanel(id="advanced_panel",
        shiny::checkboxInput("exclude_repeats", label="Exclude repeats", value=T),
        shiny::fluidRow(
          shiny::column(2, shiny::checkboxInput("exclude_bait_region", label="Exclude bait region", value=T)),
          shiny::column(10, shiny::numericInput("bait_region", label="", value=500000))
        ),
        shiny::numericInput("junsize", label="Junction size", value=300),
        shiny::numericInput("extsize", label="extsize", value=2000),
        shiny::numericInput("slocal", label="slocal", value=1000),
        shiny::numericInput("llocal", label="llocal", value=10000000)
      )
    ),

    shiny::mainPanel(
      type="tabs",
      shiny::tabsetPanel(
        shiny::tabPanel(
          "Overview",
          shiny::plotOutput("circos_output") %>% shinycssloaders::withSpinner(color="#0dc5c1"),
            shiny::actionButton("circos_calculate", label="Calculate")
        )
      )
    )
  )
)