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

options(stringsAsFactors = FALSE)

#Nayla's WD
setwd("/Users/naylaboorady/Desktop/MSK/scDNA_myeloid/shiny/scDNA_myeloid_shiny")
#Bobby's WD
#setwd("/Users/bowmanr/Projects/scDNA/scDNA_myeloid/shiny/scDNA_myeloid_shiny/")




# SAMPLE CLONALITY DATA ====
test<-read.csv("data/for_NB.csv")
test$Final_group<- factor(test$Final_group,levels=c("CH","MPN","Signaling","DTAI","DTAI-RAS","DTAI-FLT3","DTAI-FLT3-RAS"))

# Number of mutations
gg_number_of_mutations<- function(selected_group) {ggplot(test%>%group_by(.data[[selected_group]])%>%
                                 summarise(mean=mean(Number_of_mutations),
                                           sd = sd(Number_of_mutations),
                                           sem = sd(Number_of_mutations)/
                                             sqrt(length(Number_of_mutations))),
                               aes(x=get(selected_group),y=mean,fill=get(selected_group)))+
  geom_bar(stat="identity",color="black")+
  geom_errorbar(aes(ymin = mean-sem, ymax = mean+sem),width=0.5,lwd=0.5)+
  theme_classic(base_size = 16)+
  ylab("Number of mutations")+xlab("")+ggtitle("")+
  scale_y_continuous(limits = c(0,9), expand = c(0, 0)) +
  theme(axis.text.x = element_text(angle=30,hjust=1)) +
  scale_fill_brewer(type="seq",palette = "Reds",aesthetics = "fill",guide=FALSE)
}
  
# Number of clones
gg_number_of_clones<- function(selected_group) {ggplot(test,aes(y=Number_of_clones,x=get(selected_group),fill=get(selected_group)))+
  geom_boxplot(outlier.shape = NA)+  
  geom_jitter(width = 0.1,size=0.5)+
  theme_classic(base_size = 16)+
  ylab("Number of clones")+
  xlab("")+
  theme(axis.text.x = element_text(angle=30,hjust=1)) +
  scale_fill_brewer(type="seq",palette = "Reds",aesthetics = "fill",guide=FALSE)
}
  
# Shannon diversity index
gg_shannon<- function(selected_group) {ggplot(test,aes(y=Shannon,x=get(selected_group),fill=get(selected_group)))+
  geom_boxplot(outlier.shape = NA)+  
  geom_jitter(width = 0.1,size=0.5)+
  theme_classic(base_size = 16)+
  ylab("Shannon diveristy index")+
  xlab("")+
  theme(axis.text.x = element_text(angle=30,hjust=1)) +
  scale_fill_brewer(type="seq",palette = "Reds",aesthetics = "fill",guide=FALSE)
}
  
# Number of mutations in each cohort
gg_Number_of_mutations_in_Dclone<- function(selected_group) {ggplot(test%>%group_by(.data[[selected_group]])%>%
                                           summarise(mean=mean(Number_of_mutations_in_dominant_clone),
                                                     sd = sd(Number_of_mutations_in_dominant_clone),
                                                     sem = sd(Number_of_mutations_in_dominant_clone)/
                                                       sqrt(length(Number_of_mutations_in_dominant_clone))),
                                         aes(x=get(selected_group),y=mean,fill=get(selected_group)))+
  geom_bar(stat="identity",color="black")+
  geom_errorbar(aes(ymin = mean-sem, ymax = mean+sem),width=0.5,lwd=0.5)+
  theme_classic(base_size = 16)+
  ylab("Number of mutations \n in dominant clone")+xlab("")+ggtitle("")+
  scale_y_continuous(limits = c(0,4.5), expand = c(0, 0)) +
  theme(axis.text.x = element_text(angle=30,hjust=1)) +
  scale_fill_brewer(type="seq",palette = "Reds",
                    aesthetics = "fill",guide=FALSE)
}

# Dominant clone size
gg_dominant_clone_size_function <- function(selected_group){
  gg_dominant_clone_size<-ggplot(test,
                                 aes(y=Dominant_clone_size,x=get(selected_group),fill=get(selected_group)))+
    geom_boxplot(outlier.shape = NA)+  
    geom_jitter(width = 0.1,size=0.5)+
    theme_classic(base_size = 16)+
    ylab("Fraction of sample \n in dominant clone")+
    xlab("")+
    theme(axis.text.x = element_text(angle=30,hjust=1)) +
    scale_fill_brewer(type="seq",palette = "Reds",aesthetics = "fill",guide=FALSE)
  return(gg_dominant_clone_size)
}

#gg_number_of_mutations$mapping


# CLONOGRAPH DATA ====
final_sample_summary<-readRDS(file="data/final_sample_summary.rds")
clone_mutations<-readRDS(file="data/clone_mutations.rds")
sample_mutations <-readRDS(file="data/sample_mutations_with_pheno.rds")

sample_list <-final_sample_summary


# Generate clonal abundance barplot
gg_clonograph <- function(sample) {
  # Extract out the sample of interest    
  clonal_abundance <-sample_list[[sample]]$Clones 
  clonal_architecture <-sample_list[[sample]]$Architecture 
  
  # Ensure the order of the clone abundance and clone architecture are the same.
  clonal_architecture$Clone <- factor(clonal_architecture$Clone, levels=rev(clonal_abundance$Clone))
  clonal_abundance$Clone <- factor(clonal_abundance$Clone, levels=levels(clonal_architecture$Clone))
  
  
  gg_clonal_barplot<-ggplot(data=clonal_abundance, aes(x=Clone, y=Count,fill=Count)) + 
    geom_col()+ 
    ggtitle(sample)+
    theme_classic(base_size=12)+
    scale_y_continuous(expand=c(0.01,0))+
    #ylim() + 
    ylab("Cell Count")+
    geom_errorbar(aes(ymin = LCI, ymax = UCI), width = 0.2)+
    scale_fill_distiller(name = "Value", palette = "Reds", direction = 1) +
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(),
          axis.line.x =element_blank(),
          legend.position = "none",
          plot.title = element_text(hjust=0.5,size=16),
          plot.margin=unit(c(0,0,0,0),"cm"))
  
  gg_heatmap<-ggplot(data=clonal_architecture,
                     aes(x=Clone, y=Mutant, fill=Genotype))+
    geom_tile() +
    scale_fill_manual(values=c("WT"=brewer.pal(7,"Reds")[1],
                               "Heterozygous"=brewer.pal(7,"Reds")[3],
                               "Homozygous"=brewer.pal(7,"Reds")[6],
                               "Unknown"="grey50"),name="Genotype")+
    theme_classic(base_size=12) +
    ylab("Mutation")+
    scale_y_discrete(limits = rev(levels(clonal_architecture$Mutant)))+
    theme(legend.position = "right", legend.direction = "vertical",
          axis.text.x = element_blank(), 
          axis.line=element_blank(),
          axis.title.x=element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin=unit(c(0,0,0,0),"cm"))
  
  return(plot_grid(gg_clonal_barplot,gg_heatmap,ncol=1,align="v",axis="lr",rel_heights=c(1,0.75)))
}


gg_clonograph_multiplot<-function(sample){
  graphs <- lapply(sample,gg_clonograph)
  grobs <- lapply(graphs,as_grob)
  return(plot_grid(plotlist = grobs,ncol=2,align="v",axis="lr"))
}




# NETWORK GRAPH ====
network_graph<-function(genes_of_interest,disease,multi_mutant_only){
  # Subset data.frame above to only the dominant clones
  dominant_clone_mutations <- clone_mutations%>%filter(Clonality=="Dominant")
  
  # For each sample, determine if the dominant clone
  sample_mutations$Match<-ifelse(sapply(sample_mutations$Sample,function(sample) {
    all(sample_mutations%>%
          filter(Sample==sample)%>%
          select(all_of(genes_of_interest))==
          dominant_clone_mutations%>%
          filter(Sample==sample)%>%
          select(all_of(genes_of_interest)))
  }) ,"Match","Absent")
  
  
  # identify sample with at least 2 DTAI mutationss
  multi_mutant<-sample_mutations%>%
    filter(grepl(`disease`,Dx))%>%
    mutate(mutations = rowSums(select(., all_of(genes_of_interest))))%>%
    filter(mutations>=2)%>%
    distinct(Sample)%>%pull(Sample)
  
  if(length(multi_mutant)<2){
    return(print("No patients Identified"))
  }
  # Identify dominant clones 
  dominant_clones_of_interest<-clone_mutations%>%filter(Sample%in%multi_mutant)%>%
    filter(Clonality=="Dominant")%>%
    select(Clone,Clone_size,Sample,all_of(genes_of_interest))%>%
    pivot_longer(cols=all_of(genes_of_interest),
                 names_to="Gene",values_to="Mutated")%>%
    filter(Mutated==1)
  
  # Now we want to know which variants are in the dominant clone, and the size of that clone. 
  # I'm sure there is a nice way to do this in dplyr, grouping on sample, but I couldn't figure it out
  # so we will use lapply.
  genes_in_each_dominant_clone<- do.call(rbind,setNames(lapply(multi_mutant,function(x){
    # Extract the genes
    dominant_variants<- dominant_clones_of_interest%>%filter(Sample==x)%>%pull(Gene)
    
    # Extract the clone size
    dominant_clone_size<- dominant_clones_of_interest%>%filter(Sample==x)%>%pull(Clone_size)
    
    # if there are more than two DTAI variants in the dominant clone make a combinatorial edgelist
    if(length(dominant_variants)>=2){
      return(setNames(data.frame(t(combn(dominant_variants,2)),dominant_clone_size[1],"Dominant"),c("to","from","size","Clonality")))} 
    # if there is only 1 mutant in the dominant clone, list it for now so we can count the mutation, 
    # but we will eventually filter it out
    else if(length(dominant_variants)==1){
      return(setNames(data.frame(t(c(dominant_variants,dominant_variants)),dominant_clone_size,"Subclone"),c("to","from","size","Clonality")))} 
    # if no DTAI mutants in the dominant clone, ignore.
    else if(length(dominant_variants)==0){
      NULL
    }
  }),multi_mutant))%>%distinct()
  
  # Now we will go for a similar process with subclones.
  sub_clones_of_interest<-clone_mutations%>%filter(Sample%in%multi_mutant)%>%
    filter(Clonality!="Dominant")%>%
    select(Clone,Clone_size,Sample,all_of(genes_of_interest))%>%
    pivot_longer(cols=all_of(genes_of_interest),
                 names_to="Gene",values_to="Mutated")%>%
    filter(Mutated==1)%>%
    # This is how we specifically select multi mutant subclone
    group_by(Clone,Sample)%>%
    add_tally()%>%filter(n>1)%>%
    ungroup()
  
  # Same process as above, but note that we decided to only plot the largest multi mutant clone.
  # Try getting rid of this and seeing how it looks.
  genes_in_each_subclone <- do.call(rbind,setNames(lapply(multi_mutant,function(x){
    subclone_variants <- sub_clones_of_interest%>%filter(Sample==x)%>%
      filter(Clone_size==max(Clone_size))%>%
      pull(Gene)
    subclone_size <- sub_clones_of_interest%>%filter(Sample==x)%>%
      filter(Clone_size==max(Clone_size))%>%
      pull(Clone_size)
    
    if(length(subclone_variants)>=2){
      return(setNames(data.frame(t(combn(rev(subclone_variants),2)),subclone_size[1],"Subclone"),c("to","from","size","Clonality")))} 
    else if(length(subclone_variants)==1){
      return(setNames(data.frame(t(c(subclone_variants,subclone_variants)),subclone_size[1],"Subclone"),c("to","from","size","Clonality")))} 
    else if(length(subclone_variants)==0){
      NULL
    }
  }),multi_mutant))%>%distinct()
  
  # Now bind these two dataframe together
  final_set<- rbind(genes_in_each_dominant_clone,genes_in_each_subclone)
  
  # And remove the edges that are self referencing. We preserve the input variable so we can represent
  # the vertex size in relation to total mutation burden in this subset of patients.
  final_set_filtered <-final_set%>%filter(to!=from)
  
  graph<-graph_from_data_frame(final_set_filtered,directed=F)%>%
    set_edge_attr("weight", value = as.numeric(final_set_filtered%>%pull(size))*3) %>%
    set_edge_attr("color", value = ifelse(final_set_filtered%>% 
                                            pull(Clonality)=="Dominant",
                                          alpha(brewer.pal(5,"Reds")[5],0.5),
                                          alpha("grey20",0.5)))
  
  mutant_counts<-table(c(as.character(final_set$to),as.character(final_set$from)))[names(V(graph))]
  scaled_mutant_counts <-mutant_counts/sum(mutant_counts)*50
  
  radian.rescale <- function(x, start=0, direction=1) {
    c.rotate <- function(x) (x + start) %% (2 * pi) * direction
    c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
  }
  
  lab.locs <- radian.rescale(x=1:length(genes_of_interest), direction=-1, start=length(genes_of_interest))
  #lab.locs[1]<- 0
  
  plot.igraph(graph,
              edge.width = E(graph)$weight,
              vertex.color=kelly(n=1+length(genes_of_interest))[-1],#brewer.pal(5,"Reds")[5],
              vertex.frame.color="black",
              vertex.size=scaled_mutant_counts[names(V(graph))], 
              vertex.label.family="Arial",
              vertex.label.font=2,
              vertex.label.color="black",
              #   vertex.label.degree=lab.locs,
              vertex.label.dist=rep(3,length(genes_of_interest)),
              layout=layout_with_graphopt)
  legend(x=1.5,y=1,legend=names(V(graph)),
         col=kelly(n=1+length(names(V(graph))))[-1],
         bty = "n", pch=20 , 
         pt.cex = 2, 
         cex =0.75, 
         text.font=2,
         text.col=kelly(n=1+length(names(V(graph))))[-1], 
         horiz = FALSE,
         ncol=2#, 
         # inset = c(0.1, 0.1)
  )
}


#SERVER ====
shinyServer(function(input,output,session) {
  
  
  plotCount <- reactive({
    as.numeric(length(input$clonoInput))
  })
  plotHeight <- reactive({ifelse(plotCount()==1,400,400 * round(plotCount()/2))})      
  plotWidth  <- reactive({ifelse(plotCount()==1,500,1000)})      
  
  
  output$sampleClonPlot <- renderPlot(switch(input$sc, oneC = gg_number_of_mutations(input$selected_group), 
                                             oneE = gg_number_of_clones(input$selected_group), 
                                             twoA = gg_shannon(input$selected_group), 
                                             twoB = gg_Number_of_mutations_in_Dclone(input$selected_group),
                                             threeA = gg_dominant_clone_size_function(input$selected_group)
  ))
  
  output$clonalBarplot <- renderPlot(  if(plotCount()==1){
                                            gg_clonograph(input$clonoInput)
                                        }
                                       else if(plotCount()>1&plotCount()<7){
                                                gg_clonograph_multiplot(input$clonoInput)
                                        } 
  ) 
  
  output$plot.ui <- renderUI({
                              plotOutput("clonalBarplot", height = plotHeight(),
                              width = plotWidth())
  })
  
  
  
  output$networkPlot <- renderPlot(network_graph(input$networkInput,disease = "AML")
  )
  
  #downloads
  # output$downloadSampleClonPlot <- downloadHandler(
  #                 filename =  function() {
  #                   paste("iris", input$var3, sep=".")
  #                 },
  #                 # content is a function with argument file. content writes the plot to the device
  #                 content = function(file) {
  #                   if(input$var3 == "png")
  #                     png(file) # open the png device
  #                   else
  #                     pdf(file) # open the pdf device
  #                     switch(input$sc, oneC = gg_number_of_mutations, 
  #                          oneE = gg_number_of_clones, 
  #                          twoA = gg_shannon, 
  #                          twoB = gg_Number_of_mutations_in_Dclone,
  #                          threeA = gg_dominant_clone_size_function(input$selected_group))
  #                     dev.off()  # turn the device off
  #     
  #                 } 
  # )
  # 
  # output$downloadClonograph <- downloadHandler(
  #                 filename =  function() {
  #                   paste("iris", input$var3, sep=".")
  #                 },
  #                 # content is a function with argument file. content writes the plot to the device
  #                 content = function(file) {
  #                   if(input$var3 == "png")
  #                     png(file) # open the png device
  #                   else
  #                     pdf(file) # open the pdf device
  #                       if(plotCount()==1){
  #                         gg_clonograph(input$clonoInput)
  #                       }
  #                       else if(plotCount()>1&plotCount()<7){
  #                         gg_clonograph_multiplot(input$clonoInput)
  #                       } 
  #                   dev.off()  # turn the device off
  #                   
  #                 } 
  # )
  

  # output$downloadNetworkGraph <- downloadHandler(
  #                 filename =  function() {
  #                   paste("iris", input$var3, sep=".")
  #                 },
  #                 # content is a function with argument file. content writes the plot to the device
  #                 content = function(file) {
  #                   if(input$var3 == "png")
  #                     png(file) # open the png device
  #                   else
  #                     pdf(file) # open the pdf device
  #                     network_graph(input$networkInput,disease = "AML")
  #                   dev.off()  # turn the device off
  # 
  #                 }
  # )
  
  
})




