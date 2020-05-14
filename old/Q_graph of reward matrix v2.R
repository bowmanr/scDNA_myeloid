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
data[data$Sample%in%c("MA1715B","MA6300A","MA9521A","MA9521B"),"Group"]<- "MPN"
data[data$Sample%in%c("E4840","E4838new"),"Group"]<- "CMML"
data[data$Sample%in%c("MA4244A","MA2725","MA9521B","MA0092A","R2715","MA1715B","MA6363B"),"Group"]<- "MF"


DTAI_AML_samples <- data%>%filter(Group=="AML")%>%pull(Sample)
DTAI_AML_samples_summary <-final_sample_summary[DTAI_AML_samples]

graph_results<-list()


graph_results  <- lapply(DTAI_AML_samples, function(i){
  print(i)
  Known_mat<-do.call(cbind,strsplit(as.character(DTAI_AML_samples_summary[[i]]$Clones$Clone),split="_"))
  rownames(Known_mat)<-colnames(DTAI_AML_samples_summary[[i]]$NGT)[1:(length(colnames(DTAI_AML_samples_summary[[i]]$NGT ))-1)]
  mode(Known_mat)<-"numeric"
  weights <-DTAI_AML_samples_summary[[i]]$Clones$Count/sum(DTAI_AML_samples_summary[[i]]$Clones$Count)*100
  names(weights)<-apply(Known_mat,2,paste,sep="_",collapse="_")
  print(paste(i, names(DTAI_AML_samples_summary)[i],dim(Known_mat)[[1]]))
  graph_results[[i]]<-create_reward_matrix(Known_mat,weights)
})

names(graph_results) <-DTAI_AML_samples
saveRDS(graph_results,file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/MDP_results.rds")

graph_results<-readRDS(file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/MDP_results.rds")


final_results<-list()
for(i in 1:length(graph_results)){
  print(names(graph_results)[i])
  final_results[[i]]<-query_initiating_mutations(graph_results[[i]])
}
names(final_results) <- names(graph_results)
saveRDS(final_results,file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/MDP_trajectory_for_each_gene.rds")



