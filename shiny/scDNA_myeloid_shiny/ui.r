#Nayla's WD
#setwd("/Users/naylaboorady/Desktop/MSK/scDNA_myeloid/shiny/scDNA_myeloid_shiny")
#Bobby's WD
#setwd("/Users/bowmanr/Projects/scDNA/scDNA_myeloid/shiny/scDNA_myeloid_shiny/")

library(shinythemes)
# RELEVANT DATA ====
test<-read.csv("data/for_NB.csv")
sample_list<-readRDS(file="data/final_sample_summary.rds")
clone_mutations<-readRDS(file="data/clone_mutations.rds")
sample_mutations <-readRDS(file="data/sample_mutations_with_pheno.rds")

#Sample Clonality Customization
test$Final_group<- factor(test$Final_group,levels=c("CH","MPN","Signaling","DTAI","DTAI-RAS","DTAI-FLT3","DTAI-FLT3-RAS"))
test$Dx<- factor(test$Dx,levels=c("AML","CH","MPN","Other","sAML","tAML"))
test$Group<- factor(test$Group, levels = c("DTAI", "DTAI-FLT3", "DTAI-FLT3-RAS", "Signaling"))
list_final <- as.list(levels(test$Final_group))
list_dx <- as.list(levels(test$Dx))
list_group <- as.list(levels(test$Group))
available_groups <- c(list_final, list_dx, list_group)


shinyUI( 	
  
  navbarPage(title = strong("AML Mutational Profiling"), windowTitle = "AML Mutational Profiling", 
             fluid = TRUE, id = "nav",inverse=FALSE,theme = shinytheme("sandstone"),
             source("tabs/sampleClonality.r", local = TRUE)$value,
             source("tabs/clonograph.r", local = TRUE)$value,
             source("tabs/networkGraph.r", local = TRUE)$value,
             tabPanel(HTML("</a></li><li><a href=\"https://www.biorxiv.org/content/10.1101/2020.02.07.938860v1\" target=\"_blank\">Paper")),
             tabPanel(HTML("</a></li><li><a href=\"https://bowmanr.github.io/scDNA_myeloid/\" target=\"_blank\">Tutorial"))
  )
)

