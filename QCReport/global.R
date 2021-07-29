library(readr)
library(stringr)
library(dplyr)
library(shinyjs)
library(GenomeInfoDb)
library(shiny)
library(BSgenome)
source("utils.R")
source("graphics.R")


options(shiny.maxRequestSize=5*1024^3)
options(shiny.sanitize.errors = TRUE)

# Remove all VennDiagram report logs
for(p in list.files(pattern="^VennDiagram.*log$", full.names=T)) {
  file.remove(p)
}

genomes_path = "./genomes"