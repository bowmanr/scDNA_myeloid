#Nayla's WD
#setwd("/Users/naylaboorady/Desktop/MSK/scDNA_myeloid/shiny/scDNA_myeloid_shiny")
#Bobby's WD
#setwd("/Users/bowmanr/Projects/scDNA/scDNA_myeloid/shiny/scDNA_myeloid_shiny/")

library(shinythemes)
shinyUI( 	
  navbarPage(title = strong("AML Mutational Profiling"), windowTitle = "AML Mutational Profiling", 
             fluid = TRUE, id = "nav",inverse=FALSE,theme = shinytheme("cerulean"),
             source("tabs/sampleClonality.r", local = TRUE)$value,
             source("tabs/clonograph.r", local = TRUE)$value,
             source("tabs/networkGraph.r", local = TRUE)$value

  )
)