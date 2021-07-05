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

  })

}