source("/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Scripts/mBio_processing_scripts.R")
library(UpSetR)

final_sample_summary<-readRDS(file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/final_sample_summary.rds")
clonal_sample_set_after_boostrap<-readRDS(file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/clonal_sample_set_after_boostrap")
DTAI_genes <- c("ASXL1","DNMT3A","TET2","IDH1","IDH2")

genes_in_each_sample<-setNames(mclapply(mc.cores = 4,final_sample_summary[clonal_sample_set_after_boostrap],function(x){
  aggregate_by<-do.call(rbind,strsplit(colnames(x$NGT),split="[_\\.]")) [,1][-length(colnames(x$NGT))]}), clonal_sample_set_after_boostrap)

list_of_mutants<-setNames(lapply(DTAI_genes,function(x){
  names(genes_in_each_sample)[do.call(c,setNames(lapply(genes_in_each_sample,function(y){x%in%y }),names(genes_in_each_sample)))]}),DTAI_genes)

data<-fromList(list_of_mutants)
rownames(data)<-unique(unlist(list_of_mutants))
data$Sample<-unique(unlist(list_of_mutants))
data$Group <- ifelse(grepl("CH",rownames(data)),"CH","AML")

DTAI_in_dominant_clone<-lapply(final_sample_summary[clonal_sample_set_after_boostrap],function(sample){
  variants<-colnames(sample$NGT)[1:(length(colnames(sample$NGT))-1)]
  genes <-do.call(rbind,strsplit(variants,split="[_\\.]"))[,1]
  clones<-sample$Clones[,"Clone"]
  WT_clone_to_exclude <-!apply(do.call(rbind,strsplit(clones,split="_")),1,function(x){all(x=="0")})
  dominant_clone <-sample$Clones%>%filter(!apply(do.call(rbind,strsplit(Clone,split="_")),1,function(x){all(x=="0")})) %>%filter(Count==max(Count))%>%pull(Clone)
  dominant_clone_genes<-genes[do.call(c,strsplit(dominant_clone,split="_"))!=0]
  unique(intersect(dominant_clone_genes,DTAI_genes))
})

list_of_mutants_in_dominant_clone<-setNames(lapply(DTAI_genes,function(x){
  names(DTAI_in_dominant_clone)[do.call(c,setNames(lapply(DTAI_in_dominant_clone,function(y){x%in%y}),names(DTAI_in_dominant_clone)))]}),DTAI_genes)

data_dom<-fromList(list_of_mutants_in_dominant_clone)
rownames(data_dom)<-unique(unlist(list_of_mutants_in_dominant_clone))
data_dom$Sample<-unique(unlist(list_of_mutants_in_dominant_clone))
data_dom$Match<-ifelse(sapply(rownames(data_dom),function(sample){all(data_dom[sample,1:5]==data[sample,1:5]) }),"Match","Absent")

test_set<-full_join(data,data_dom[,6:7],by="Sample")
test_set$Match[is.na(test_set$Match)]<-"Absent"

test_set[test_set$Sample%in%c("MA1715B","MA6300A","MA9521A","MA9521B"),"Group"]<- "MPN"
test_set[test_set$Sample%in%c("E4840","E4838new"),"Group"]<- "CMML"
test_set[test_set$Sample%in%c("MA4244A","MA2725","MA9521B","MA0092A","R2715","MA1715B","MA6363B"),"Group"]<- "MF"


test_set%>%filter(Group=="AML")%>%filter((ASXL1+DNMT3A+TET2+IDH1+IDH2)==1)%>%distinct(Sample)%>%summarize(Count=n())
multi_DTAI<-test_set%>%filter(Group=="AML")%>%filter((ASXL1+DNMT3A+TET2+IDH1+IDH2)>=2)%>%distinct(Sample)%>%pull(Sample)


genes_in_each_dominant_clone<-lapply(final_sample_summary[multi_DTAI],function(sample){
  variants<-colnames(sample$NGT)[1:(length(colnames(sample$NGT))-1)]
  genes <-do.call(rbind,strsplit(variants,split="[_\\.]"))[,1]
  clones<-sample$Clones[,"Clone"]
  WT_clone_to_exclude <-!apply(do.call(rbind,strsplit(clones,split="_")),1,function(x){all(x=="0")})
  dominant_clone <-sample$Clones%>%filter(!apply(do.call(rbind,strsplit(Clone,split="_")),1,function(x){all(x=="0")})) %>%filter(Count==max(Count))%>%pull(Clone)
  
  dominant_clone_size <-sample$Clones%>%filter(Clone==dominant_clone)%>%pull(Count)  /sum(sample$Clones$Count)
  dominant_clone_genes<-unique(intersect(genes[do.call(c,strsplit(dominant_clone,split="_"))!=0],DTAI_genes))
return(list(dominant_clone_genes,dominant_clone_size))
  })


DTAI_dom_clone_edge_list<-lapply(genes_in_each_dominant_clone,function(y){
  x<-y[[1]]
  if(length(x)>=2){data.frame(t(combn(x,2)),y[[2]])} 
  else if(length(x)==1){data.frame( t(c(x,x) ),y[[2]] )} 
  else if(length(x)==0){NULL}
})
DTAI_dom_clone_edge_mat<-data.frame(setNames(do.call(rbind,DTAI_dom_clone_edge_list),c("to","from","size")),
                                   "Clonality"="Dominant")




genes_in_each_subclone<-lapply(final_sample_summary[multi_DTAI],function(sample){
  variants<-colnames(sample$NGT)[1:(length(colnames(sample$NGT))-1)]
  genes <-do.call(rbind,strsplit(variants,split="[_\\.]"))[,1]
  clones<-sample$Clones[,"Clone"]
  WT_clone_to_exclude <-!apply(do.call(rbind,strsplit(clones,split="_")),1,function(x){all(x=="0")})
  subclones <-sample$Clones%>%filter(!apply(do.call(rbind,strsplit(Clone,split="_")),1,function(x){all(x=="0")})) %>%filter(Count!=max(Count))%>%pull(Clone)
  subclone_genes<-lapply(subclones,function(clone){
          subclone_genes <-unique(intersect(genes[do.call(c,strsplit(clone,split="_"))!=0],DTAI_genes))
  })
  names(subclone_genes)<- subclones
  multi_mutant_subclones <-names(subclone_genes)[do.call(c,lapply(subclone_genes,function(x){length(x)>1})) ]
  if(length(multi_mutant_subclones)>0){
   largest_multi_mutant_subclone <-sample$Clones%>%filter(Clone%in%multi_mutant_subclones)%>%filter(Count==max(Count))%>%pull(Clone)
   subclone_size <-sample$Clones%>%filter(Clone==largest_multi_mutant_subclone)%>%pull(Count)  /sum(sample$Clones$Count)
   largest_multi_mutant_subclone_genes<-unique(intersect(genes[do.call(c,strsplit(largest_multi_mutant_subclone,split="_"))!=0],DTAI_genes))
  return(list(largest_multi_mutant_subclone_genes,subclone_size))
  } else{
    NULL
  }
})

DTAI_subclone_edge_list<-lapply(genes_in_each_subclone,function(y){
  x<-y[[1]]
  if(length(x)>=2){data.frame(t(combn(x,2)),y[[2]])} 
  else if(length(x)==1){data.frame( t(c(x,x) ),y[[2]] )} 
  else if(length(x)==0){NULL}
})
DTAI_subclone_edge_mat<-data.frame(setNames(do.call(rbind,DTAI_subclone_edge_list),c("to","from","size")),
                                     "Clonality"="Subclone")

final_set <- rbind(DTAI_dom_clone_edge_mat,DTAI_subclone_edge_mat)

final_set<-final_set%>%filter(to!=from)

graph<-graph_from_data_frame(final_set,directed=F)%>%
                    set_edge_attr("weight", value = final_set$size*3) %>%
                    set_edge_attr("color", value =ifelse(final_set$Clonality=="Dominant",  brewer.pal(5,"Reds")[5],"grey20"))
        
mutant_counts<-c(table(final_set[,1])+table(final_set[,2]))[names(V(graph))]
scaled_mutant_counts <-mutant_counts/sum(mutant_counts)*50

radian.rescale <- function(x, start=0, direction=1) {
  c.rotate <- function(x) (x + start) %% (2 * pi) * direction
  c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
}

lab.locs <- radian.rescale(x=1:5, direction=-1, start=5)
lab.locs[3]<- -2.5



tiff(width=5,height=5,res=300, units = "in",file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/Figures/F2E-network_DTAI.tiff") # or other device

plot.igraph(graph,edge.width = E(graph)$weight, vertex.color=brewer.pal(5,"Reds")[5],
            vertex.frame.color=brewer.pal(5,"Reds")[5],
            vertex.size=scaled_mutant_counts, vertex.label.family="Helvetica",
            vertex.label.color="black",vertex.label.degree=lab.locs,vertex.label.dist=c(3,4,3,7,3),
            layout=layout_in_circle)

dev.off()