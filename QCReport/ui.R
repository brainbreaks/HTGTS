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

      shiny::actionLink("advanced_hide", "advanced options"),
      shiny::wellPanel(id="advanced_panel",
        shiny::checkboxInput("exclude_repeats", label="Exclude repeats", value=F),
        shiny::fluidRow(
          shiny::column(2, shiny::checkboxInput("exclude_bait_region", label="Exclude bait region", value=F)),
          shiny::column(10, shiny::textInput("bait_region", label="", value=500000))
        ),
        shiny::textInput("junsize", label="Junction size", value=300),
        shiny::textInput("extsize", label="extsize", value=2000),
        shiny::textInput("slocal", label="slocal", value=1000),
        shiny::textInput("llocal", label="llocal", value=10000000)
      )
    ),

    shiny::mainPanel(
      type="tabs",
      shiny::tabsetPanel(
        shiny::tabPanel(
          "Overview",
          shiny::plotOutput("overview_output"),
            shiny::actionButton("overview_calculate", label="Calculate")
        )
      )
    )
  )
)