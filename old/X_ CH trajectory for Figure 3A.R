source("/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Scripts/mBio_processing_scripts.R")
library(UpSetR)

final_sample_summary<-readRDS(file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/final_sample_summary.rds")
clonal_sample_set_after_boostrap<-readRDS(file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/clonal_sample_set_after_boostrap")
DTAI_genes <- c("ASXL1","DNMT3A","TET2","IDH1","IDH2")


CH_samples <-grep("^CH",names(final_sample_summary),value=TRUE)
CH_samples_summary <-final_sample_summary[CH_samples]

graph_results<-list()


graph_results  <- lapply(CH_samples, function(i){
  print(i)
  Known_mat<-do.call(cbind,strsplit(as.character(CH_samples_summary[[i]]$Clones$Clone),split="_"))
  rownames(Known_mat)<-colnames(CH_samples_summary[[i]]$NGT)[1:(length(colnames(CH_samples_summary[[i]]$NGT ))-1)]
  mode(Known_mat)<-"numeric"
  weights <-CH_samples_summary[[i]]$Clones$Count/sum(CH_samples_summary[[i]]$Clones$Count)*100
  names(weights)<-apply(Known_mat,2,paste,sep="_",collapse="_")
  print(paste(i, names(CH_samples_summary)[i],dim(Known_mat)[[1]]))
  graph_results[[i]]<-create_reward_matrix(Known_mat,weights)
})

names(graph_results) <-CH_samples
saveRDS(graph_results,file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/MDP_results_on_CH.rds")

graph_results<-readRDS(file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/MDP_results_on_CH.rds")


final_results<-list()
for(i in 1:length(graph_results)){
  print(names(graph_results)[i])
  final_results[[i]]<-query_initiating_mutations(graph_results[[i]])
}
names(final_results) <- names(graph_results)
saveRDS(final_results,file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/MDP_trajectory_for_each_gene_in_CH.rds")

names(final_results)


plot_optimal_graph_for_trajectory("CH1414")

"CH1414","CH7557","CH9901new"

final_sample_summary[CH_samples][["CH7557"]]$Clones

colnames(final_sample_summary[CH_samples][["CH9901new"]]$NGT)


CH_samples[9]
plot_optimal_graph_for_trajectory<-function(sample){
  all_results <-final_results[[sample]]
  
  all_results_filtered<-lapply(all_results,function(storage_results){
    storage_results[lapply(storage_results,function(x){sum(x$reward,na.rm = TRUE)})==0]<-NULL
    storage_results<-lapply(storage_results,function(x){x[1:which.max(x$cumulative_reward),]})
    final<-do.call(rbind,storage_results)[,c(1,4,3)]
    nodes_to_remove <- setdiff(setdiff(final$nextState,final$states),final_sample_summary[[sample]]$Clones$Clone)
    final <- final%>%filter(!nextState%in%nodes_to_remove)
    
    nodes_to_remove <- setdiff(setdiff(final$nextState,final$states),final_sample_summary[[sample]]$Clones$Clone)
    final <- final%>%filter(!nextState%in%nodes_to_remove)
    final[!duplicated(final),]
  })
  
  all_results_filtered<-setNames(lapply(names(all_results_filtered),function(x){
    data.frame(all_results_filtered[[x]],
               "mutation"=x,
               "edge"=paste(all_results_filtered[[x]]$states,all_results_filtered[[x]]$nextState,sep="->"))
  }),names(all_results))
  
  optimal<-names(which.max(do.call(c,lapply(all_results_filtered,function(x){
    sum(x$reward)
  }))))
  print(optimal)
  final<-do.call(rbind,all_results_filtered)#%>%filter(mutation==optimal)
  final%>%group_by(mutation)%>%distinct(nextState,reward)%>%summarize(total=sum(reward))
  
  graph<-graph_from_data_frame(final,directed=T)
  weights<-final_sample_summary[[sample]]$Clones$Count/sum(final_sample_summary[[sample]]$Clones$Count)
  names(weights) <-final_sample_summary[[sample]]$Clones$Clone
  weight_subset<-weights[names(weights)%in%names(V(graph))]
  nodes_to_add_names<-setdiff(names(V(graph)),names(weights))
  nodes_to_add <- rep(0.1,length(nodes_to_add_names))
  names(nodes_to_add)<-nodes_to_add_names
  weight_final <- c(weight_subset,nodes_to_add)[names(V(graph))]
  
  mutations<-names(all_results_filtered)
  
  
  clone_colors <-ifelse(names(V(graph))%in%final_sample_summary[[sample]]$Clones$Clone,brewer.pal(5,"Reds")[5],"grey80")
  
  #par(mar=c(0,0,0,0))
  #pdf(width=8,height=8,file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/Figures/F2E-network_DTAI_test.pdf") # or other device
  #    quartz()
  (plot(graph,layout=layout_as_tree,
        vertex.color=ifelse(names(V(graph))%in%final_sample_summary[[sample]]$Clones$Clone,brewer.pal(5,"Reds")[5],"grey80"),
        vertex.frame.color=ifelse(names(V(graph))%in%final_sample_summary[[sample]]$Clones$Clone ,brewer.pal(5,"Reds")[5],"grey80"),
        vertex.size=(weight_final*50),#main=c(optimal),annotate.plot=TRUE,
        vertex.label=NA,
        edge.color="grey30",#ifelse(edge_attr(graph)$mutation%in%mutations[1],brewer.pal(5,"Oranges")[4],
        #   ifelse(edge_attr(graph)$mutation%in%mutations[2],brewer.pal(5,"Blues")[4],brewer.pal(5,"Greens")[4])), 
        edge.arrow.size=0.5,arrow.width=0.5))
  
}
