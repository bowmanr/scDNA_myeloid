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
setwd("/Users/naylaboorady/Desktop/MSK/scDNA_myeloid/scDNA_myeloid_shiny/")
#Bobby's WD
#setwd("/Users/bowmanr/Projects/scDNA/scDNA_myeloid/shiny/scDNA_myeloid/")



# SAMPLE CLONALITY DATA ====
test<-read.csv("data/for_NB.csv")
test$Final_group<- factor(test$Final_group,levels=c("CH","MPN","Signaling","DTAI","DTAI-RAS","DTAI-FLT3","DTAI-FLT3-RAS"))

# Number of mutations
gg_number_of_mutations<-ggplot(test%>%group_by(Final_group)%>%
                                 summarise(mean=mean(Number_of_mutations),
                                           sd = sd(Number_of_mutations),
                                           sem = sd(Number_of_mutations)/
                                             sqrt(length(Number_of_mutations))),
                               aes(x=Final_group,y=mean,fill=Final_group))+
  geom_bar(stat="identity",color="black")+
  geom_errorbar(aes(ymin = mean-sem, ymax = mean+sem),width=0.5,lwd=0.5)+
  theme_classic(base_size = 16)+
  ylab("Number of mutations")+xlab("")+ggtitle("")+
  scale_y_continuous(limits = c(0,9), expand = c(0, 0)) +
  theme(axis.text.x = element_text(angle=30,hjust=1)) +
  scale_fill_brewer(type="seq",palette = "Reds",aesthetics = "fill",guide=FALSE)

# Number of clones
gg_number_of_clones<-ggplot(test,aes(y=Number_of_clones,x=Final_group,fill=Final_group))+
  geom_boxplot(outlier.shape = NA)+  
  geom_jitter(width = 0.1,size=0.5)+
  theme_classic(base_size = 16)+
  ylab("Number of clones")+
  xlab("")+
  theme(axis.text.x = element_text(angle=30,hjust=1)) +
  scale_fill_brewer(type="seq",palette = "Reds",aesthetics = "fill",guide=FALSE)

# Shannon diversity index
gg_shannon<-ggplot(test,aes(y=Shannon,x=Final_group,fill=Final_group))+
  geom_boxplot(outlier.shape = NA)+  
  geom_jitter(width = 0.1,size=0.5)+
  theme_classic(base_size = 16)+
  ylab("Shannon diveristy index")+
  xlab("")+
  theme(axis.text.x = element_text(angle=30,hjust=1)) +
  scale_fill_brewer(type="seq",palette = "Reds",aesthetics = "fill",guide=FALSE)

# Number of mutations in each cohort
gg_Number_of_mutations_in_Dclone<-ggplot(test%>%group_by(Final_group)%>%
                                           summarise(mean=mean(Number_of_mutations_in_dominant_clone),
                                                     sd = sd(Number_of_mutations_in_dominant_clone),
                                                     sem = sd(Number_of_mutations_in_dominant_clone)/
                                                       sqrt(length(Number_of_mutations_in_dominant_clone))),
                                         aes(x=Final_group,y=mean,fill=Final_group))+
  geom_bar(stat="identity",color="black")+
  geom_errorbar(aes(ymin = mean-sem, ymax = mean+sem),width=0.5,lwd=0.5)+
  theme_classic(base_size = 16)+
  ylab("Number of mutations \n in dominant clone")+xlab("")+ggtitle("")+
  scale_y_continuous(limits = c(0,4.5), expand = c(0, 0)) +
  theme(axis.text.x = element_text(angle=30,hjust=1)) +
  scale_fill_brewer(type="seq",palette = "Reds",
                    aesthetics = "fill",guide=FALSE)

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

gg_number_of_mutations$mapping


# CLONOGRAPH DATA ====
final_sample_summary<-readRDS(file="data/final_sample_summary.rds")

#sample <-input$clonoInput
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



# NETWORK GRAPH DATA: 2E ====
#final_sample_summary<-readRDS(file="./data/final_sample_summary.rds")
pheno<-readRDS(file="./data/pheno.rds")

sample_mutations<-do.call(rbind,lapply(names(final_sample_summary),function(sample){
  data.frame("Sample"=sample,
             "DNMT3A"=ifelse(any(grepl("DNMT3A",colnames(final_sample_summary[[sample]]$NGT))),1,0),
             "TET2"=ifelse(any(grepl("TET2",colnames(final_sample_summary[[sample]]$NGT))),1,0),
             "IDH2"=ifelse(any(grepl("IDH2",colnames(final_sample_summary[[sample]]$NGT))),1,0),
             "IDH1"=ifelse(any(grepl("IDH1",colnames(final_sample_summary[[sample]]$NGT))),1,0),
             "ASXL1"=ifelse(any(grepl("ASXL1",colnames(final_sample_summary[[sample]]$NGT))),1,0),
             "FLT3"=ifelse(any(grepl("FLT3",colnames(final_sample_summary[[sample]]$NGT))),1,0),
             "JAK2"=ifelse(any(grepl("JAK2",colnames(final_sample_summary[[sample]]$NGT))),1,0),
             "NRAS"=ifelse(any(grepl("NRAS",colnames(final_sample_summary[[sample]]$NGT))),1,0),
             "KRAS"=ifelse(any(grepl("KRAS",colnames(final_sample_summary[[sample]]$NGT))),1,0),
             "PTPN11"=ifelse(any(grepl("PTPN11",colnames(final_sample_summary[[sample]]$NGT))),1,0)
  )
}))

clone_mutations<-do.call(rbind,lapply(names(final_sample_summary),function(sample){
  
  # select the clones
  clones<-final_sample_summary[[sample]]$Clones%>%select(Clone)
  # Extract the mutations
  mutations <- colnames(final_sample_summary[[sample]]$NGT%>%select(!Clone))
  
  out<-final_sample_summary[[sample]]$Clones%>%
                  mutate(Clone_size=Count/sum(Count))%>%
                  select(Clone,Clone_size)%>%
                  separate(col=Clone,
                          into=mutations,sep="_",
                          remove=FALSE)%>%
                 pivot_longer(cols=mutations,
                               names_to="Variant",
                              values_to="Genotype")%>%
                  add_column(Sample=`sample`)%>%
                  group_by(Clone)%>%
                  mutate(WT=ifelse(all(Genotype==0),1,0))%>%
                  filter(WT==0)%>%
                  filter(Genotype!=0)%>%
                  ungroup()%>%
                  mutate(Clonality=ifelse(Clone_size==max(Clone_size),
                                          "Dominant","Subclone"))%>%
                  group_by(Clone)%>%
                  mutate(Gene=do.call(rbind,strsplit(Variant,"[\\._]"))[,1])%>%
                  mutate(DNMT3A=ifelse(any(Gene%in%"DNMT3A"),1,0),
                        TET2=ifelse(any(Gene%in%"TET2"),1,0),
                        ASXL1=ifelse(any(Gene%in%"ASXL1"),1,0),
                        IDH1=ifelse(any(Gene%in%"IDH1"),1,0),
                        IDH2=ifelse(any(Gene%in%"IDH2"),1,0),
                        FLT3=ifelse(any(Gene%in%"FLT3"),1,0),
                        NRAS=ifelse(any(Gene%in%"NRAS"),1,0),
                        KRAS=ifelse(any(Gene%in%"KRAS"),1,0),
                        PTPN11=ifelse(any(Gene%in%"PTPN11"),1,0),
                        JAK2=ifelse(any(Gene%in%"JAK2"),1,0))%>%
                 ungroup()%>%
                  select(!c(Variant,Genotype,Gene))%>%
                  distinct()
}))


# identify sample with at least 2 DTAI mutationss
multi_DTAI<-test_set%>%filter(grepl("AML",Dx))%>%
  filter((ASXL1+DNMT3A+TET2+IDH1+IDH2)>=2)%>%
  distinct(Sample)%>%pull(Sample)

# Identify dominant clones 
DTAI_dominant_clones<-clone_mutations%>%filter(Sample%in%multi_DTAI)%>%
  filter(Clonality=="Dominant")%>%
  select(Clone,Clone_size,Sample,DNMT3A,TET2,ASXL1,IDH1,IDH2)%>%
  pivot_longer(cols=c(DNMT3A,TET2,ASXL1,IDH1,IDH2),
               names_to="Gene",values_to="Mutated")%>%
  filter(Mutated==1)

# Now we want to know which variants are in the dominant clone, and the size of that clone. 
# I'm sure there is a nice way to do this in dplyr, grouping on sample, but I couldn't figure it out
# so we will use lapply.
genes_in_each_dominant_clone<- do.call(rbind,setNames(lapply(multi_DTAI,function(x){
  
  # Extract the genes
  dominant_variants<- DTAI_dominant_clones%>%filter(Sample==x)%>%pull(Gene)
  
  # Extract the clone size
  dominant_clone_size<- DTAI_dominant_clones%>%filter(Sample==x)%>%pull(Clone_size)
  
  # if there are more than two DTAI variants in the dominant clone make a combinatorial edgelist
  if(length(dominant_variants)>=2){
    return(setNames(data.frame(t(combn(dominant_variants,2)),dominant_clone_size,"Dominant"),c("to","from","size","Clonality")))} 
  # if there is only 1 mutant in the dominant clone, list it for now so we can count the mutation, 
  # but we will eventually filter it out
  else if(length(dominant_variants)==1){
    return(setNames(data.frame(t(c(dominant_variants,dominant_variants)),dominant_clone_size,"Subclone"),c("to","from","size","Clonality")))} 
  # if no DTAI mutants in the dominant clone, ignore.
  else if(length(dominant_variants)==0){
    NULL
  }
}),multi_DTAI))%>%distinct()

# Now we will go for a similar process with subclones.
DTAI_sub_clones<-clone_mutations%>%filter(Sample%in%multi_DTAI)%>%
  filter(Clonality!="Dominant")%>%
  select(Clone,Clone_size,Sample,DNMT3A,TET2,ASXL1,IDH1,IDH2)%>%
  pivot_longer(cols=c(DNMT3A,TET2,ASXL1,IDH1,IDH2),
               names_to="Gene",values_to="Mutated")%>%
  filter(Mutated==1)%>%
  # This is how we specifically select multi mutant subclone
  group_by(Clone,Sample)%>%
  add_tally()%>%filter(n>1)%>%
  ungroup()

# Same process as above, but note that we decided to only plot the largest multi mutant clone.
# Try getting rid of this and seeing how it looks.
genes_in_each_subclone <- do.call(rbind,setNames(lapply(multi_DTAI,function(x){
  subclone_variants <- DTAI_sub_clones%>%filter(Sample==x)%>%
    filter(Clone_size==max(Clone_size))%>%
    pull(Gene)
  subclone_size <- DTAI_sub_clones%>%filter(Sample==x)%>%
    filter(Clone_size==max(Clone_size))%>%
    pull(Clone_size)
  
  if(length(subclone_variants)>=2){
    return(setNames(data.frame(t(combn(rev(subclone_variants),2)),subclone_size,"Subclone"),c("to","from","size","Clonality")))} 
  else if(length(subclone_variants)==1){
    return(setNames(data.frame(t(c(subclone_variants,subclone_variants)),subclone_size,"Subclone"),c("to","from","size","Clonality")))} 
  else if(length(subclone_variants)==0){
    NULL
  }
}),multi_DTAI))%>%distinct()

# Now bind these two dataframe together
final_set<- rbind(genes_in_each_dominant_clone,genes_in_each_subclone)

# And remove the edges that are self referencing. We preserve the input variable so we can represent
# the vertex size in relation to total mutation burden in this subset of patients.
final_set_filtered <-final_set%>%filter(to!=from)

graph<-graph_from_data_frame(final_set_filtered,directed=F)%>%
  set_edge_attr("weight", value = as.numeric(final_set_filtered%>%pull(size))*3) %>%
  set_edge_attr("color", value = ifelse(final_set_filtered%>% 
                                          pull(Clonality)=="Dominant",
                                        brewer.pal(5,"Reds")[5],"grey20"))

mutant_counts<-table(c(as.character(final_set$to),as.character(final_set$from)))[names(V(graph))]
scaled_mutant_counts <-mutant_counts/sum(mutant_counts)*50

radian.rescale <- function(x, start=0, direction=1) {
  c.rotate <- function(x) (x + start) %% (2 * pi) * direction
  c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
}

lab.locs <- radian.rescale(x=1:5, direction=-1, start=5)
lab.locs[3]<- -2.5

reordered_graph<-igraph::permute(graph,c(4,3,2,1,5))





# NETWORK GRAPH DATA: 2F ====
multi_signaling<-test_set%>%filter(grepl("AML",Dx))%>%
  filter((FLT3+JAK2+NRAS+KRAS+PTPN11)>=2)%>%
  distinct(Sample)%>%pull(Sample)

signaling_dominant_clones<-clone_mutations%>%filter(Sample%in%multi_signaling)%>%
  filter(Clonality=="Dominant")%>%
  select(Clone_size,Sample,FLT3,JAK2,NRAS,KRAS,PTPN11)%>%
  pivot_longer(cols=c(FLT3,JAK2,NRAS,KRAS,PTPN11),
               names_to="Gene",values_to="Mutated")%>%
  filter(Mutated==1)

genes_in_each_dominant_clone<- do.call(rbind,setNames(lapply(multi_signaling,function(x){
  dominant_variants<- signaling_dominant_clones%>%filter(Sample==x)%>%pull(Gene)
  dominant_clone_size<- signaling_dominant_clones%>%filter(Sample==x)%>%pull(Clone_size)
  
  if(length(dominant_variants)>=2){
    return(setNames(data.frame(t(combn((dominant_variants),2)),dominant_clone_size,"Dominant"),c("to","from","size","Clonality")))} 
  else if(length(dominant_variants)==1){
    return(setNames(data.frame(t(c(dominant_variants,dominant_variants)),dominant_clone_size,"Subclone"),c("to","from","size","Clonality")))} 
  else if(length(dominant_variants)==0){
    NULL
  }
}),multi_signaling))%>%distinct()

signaling_sub_clones<-clone_mutations%>%filter(Sample%in%multi_signaling)%>%
  filter(Clonality!="Dominant")%>%
  select(Clone,Clone_size,Sample,FLT3,JAK2,NRAS,KRAS,PTPN11)%>%
  pivot_longer(cols=c(FLT3,JAK2,NRAS,KRAS,PTPN11),
               names_to="Gene",values_to="Mutated")%>%
  filter(Mutated==1)%>%
  group_by(Clone,Sample)%>%
  add_tally()%>%filter(n>1)%>%
  ungroup()

genes_in_each_subclone<- do.call(rbind,setNames(lapply(multi_signaling,function(x){
  subclone_variants<- signaling_sub_clones%>%filter(Sample==x)%>%
    filter(Clone_size==max(Clone_size))%>%pull(Gene)
  subclone_size<- signaling_sub_clones%>%filter(Sample==x)%>%
    filter(Clone_size==max(Clone_size))%>%pull(Clone_size)
  
  if(length(subclone_variants)>=2){
    return(setNames(data.frame(t(combn((subclone_variants),2)),subclone_size,"Subclone"),c("to","from","size","Clonality")))} 
  else if(length(subclone_variants)==1){
    return(setNames(data.frame(t(c(subclone_variants,subclone_variants)),subclone_size,"Subclone"),c("to","from","size","Clonality")))} 
  else if(length(subclone_variants)==0){
    NULL
  }
}),multi_signaling))%>%distinct()


final_set<- rbind(genes_in_each_dominant_clone,genes_in_each_subclone)

final_set_filtered <-final_set%>%filter(to!=from)

graph<-graph_from_data_frame(final_set_filtered,directed=F)%>%
  set_edge_attr("weight", value = as.numeric(final_set_filtered%>%pull(size))*3) %>%
  set_edge_attr("color", value = ifelse(final_set_filtered%>% 
                                          pull(Clonality)=="Dominant",
                                        brewer.pal(5,"Reds")[5],"grey20"))

mutant_counts<-table(c(as.character(final_set$to),as.character(final_set$from)))[names(V(graph))]
scaled_mutant_counts <-mutant_counts/sum(mutant_counts)*50
radian.rescale <- function(x, start=0, direction=1) {
  c.rotate <- function(x) (x + start) %% (2 * pi) * direction
  c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
}

lab.locs <- radian.rescale(x=1:5, direction=-1, start=5)
lab.locs[3]<- -2.5

reordered_graph<-igraph::permute(graph,c(5,2,4,1,3))







#SERVER ====
shinyServer(function(input,output,session) {
  
  
  plotCount <- reactive({
    as.numeric(length(input$clonoInput))
  })
  plotHeight <- reactive({ifelse(plotCount()==1,400,400 * round(plotCount()/2))})      
  plotWidth  <- reactive({ifelse(plotCount()==1,500,1000)})      
  
  
  output$sampleClonPlot <- renderPlot(switch(input$sc, oneC = gg_number_of_mutations, 
                                             oneE = gg_number_of_clones, 
                                             twoA = gg_shannon, 
                                             twoB = gg_Number_of_mutations_in_Dclone,
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
  
  
  
  output$network2e <- renderPlot(plot.igraph(reordered_graph,
                                             edge.width = E(reordered_graph)$weight,
                                             vertex.color=brewer.pal(5,"Reds")[5],
                                             vertex.frame.color=brewer.pal(5,"Reds")[5],
                                             vertex.size=scaled_mutant_counts[names(V(reordered_graph))], 
                                             vertex.label.family="Helvetica",
                                             vertex.label.color="black",
                                             vertex.label.degree=lab.locs,
                                             vertex.label.dist=c(3,4,3,7,3),
                                             layout=layout_in_circle)
  )
  
  
  output$network2f  <- renderPlot(plot.igraph(reordered_graph,
                                              edge.width = E(reordered_graph)$weight,
                                              vertex.color=brewer.pal(5,"Reds")[5],
                                              vertex.frame.color=brewer.pal(5,"Reds")[5],
                                              vertex.size=scaled_mutant_counts[names(V(reordered_graph))], 
                                              vertex.label.family="Helvetica",
                                              vertex.label.color="black",
                                              vertex.label.degree=lab.locs,
                                              vertex.label.dist=c(3,4,3,7,3),
                                              layout=layout_in_circle)
  )
  
  
})




