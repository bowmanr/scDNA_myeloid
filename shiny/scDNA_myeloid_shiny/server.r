library(shiny)
library(magrittr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(igraph)
library(pals)
library(UpSetR)

#Nayla's WD
#setwd("/Users/naylaboorady/Desktop/MSK/scDNA_myeloid/shiny/scDNA_myeloid_shiny")
#Bobby's WD
#setwd("/Users/bowmanr/Projects/scDNA/scDNA_myeloid/shiny/scDNA_myeloid_shiny/")

source("functions.r")

# RELEVANT DATA ====
test<-read.csv("data/for_NB.csv")
test$Final_group<- factor(test$Final_group,levels=c("CH","MPN","Signaling","DTAI","DTAI-RAS","DTAI-FLT3","DTAI-FLT3-RAS"))
final_sample_summary<-readRDS(file="data/final_sample_summary.rds")
clone_mutations<-readRDS(file="data/clone_mutations.rds")
sample_mutations <-readRDS(file="data/sample_mutations_with_pheno.rds")

sample_list <-final_sample_summary


#SERVER ====
shinyServer(function(input,output,session) {
  
  plotCount_sample <- reactive({as.numeric(length(input$clonoInputSample))})
  plotHeight_sample <- reactive({ifelse(plotCount_sample()==1,400,400 * round(plotCount_sample()/2))})      
  plotWidth_sample  <- reactive({ifelse(plotCount_sample()==1,500,1000)})      

  plotCount_gene <- reactive({as.numeric(length(gene_sample_selection(input$clonoInputGene) ))})
  plotHeight_gene <- reactive({ifelse(plotCount_gene()==1,400,400 * round(plotCount_gene()/2))})      
  plotWidth_gene  <- reactive({ifelse(plotCount_gene()==1,500,1000)})      
  
  

  output$sampleClonPlotP <- renderPlot(switch(input$sc, oneC = gg_number_of_mutations(input$sampleClonGroups), 
                                             oneE = gg_number_of_clones(input$sampleClonGroups), 
                                             twoA = gg_shannon(input$sampleClonGroups), 
                                             twoB = gg_Number_of_mutations_in_Dclone(input$sampleClonGroups),
                                             threeA = gg_dominant_clone_size_function(input$sampleClonGroups)
                             ))
  
  output$sampleClonPlotC <- renderPlot(switch(input$sc, oneC = gg_number_of_mutations(input$sampleClonGroups), 
                                              oneE = gg_number_of_clones(input$sampleClonGroups), 
                                              twoA = gg_shannon(input$sampleClonGroups), 
                                              twoB = gg_Number_of_mutations_in_Dclone(input$sampleClonGroups),
                                              threeA = gg_dominant_clone_size_function(input$sampleClonGroups)
  ))
  
  output$clonalBarplot <- renderPlot(switch(input$selection_Feature,
                                            Sample= if(plotCount_sample()==1){
                                                          gg_clonograph(input$clonoInputSample)
                                                      }
                                                     else if(plotCount_sample()>1&plotCount_sample()<7){
                                                              gg_clonograph_multiplot(input$clonoInputSample)
                                                      } ,
                                            Gene =if(plotCount_gene()==1){
                                                        gg_clonograph(input$clonoInputGene)
                                                      }
                                                      else if(plotCount_gene()>1&plotCount_gene()<25){
                                                        gg_clonograph_multiplot(input$clonoInputGene)
                                                      })
  ) 
  
  output$plot.ui <- renderUI({switch(input$selection_Feature,
                                     Sample=plotOutput("clonalBarplot", height = plotHeight_sample(),
                                                       width = plotHeight_sample()),
                                     Gene =plotOutput("clonalBarplot", height = plotHeight_gene(),
                                                      width = plotHeight_gene())
                              )})
  
  
  
  output$networkPlot <- renderPlot(network_graph(input$networkInput,disease = "AML"))
  
})




