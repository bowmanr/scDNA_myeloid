source("/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Scripts/mBio_processing_scripts.R")

final_sample_summary<-readRDS(file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/final_sample_summary.rds")
clonal_sample_set_after_boostrap<-readRDS(file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/clonal_sample_set_after_boostrap")

graph_results<-readRDS(file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/MDP_results.rds")
final_results<-readRDS(file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/MDP_trajectory_for_each_gene.rds")


output<-lapply(names(final_results),function(sample){
  print(sample)
all_results <-final_results[[sample]]

all_results_filtered<-lapply(all_results,function(storage_results){
  storage_results[lapply(storage_results,function(x){sum(x$reward,na.rm = TRUE)})==0]<-NULL
  if(length(storage_results)==0){
    return(NULL)
  }
  storage_results<-lapply(storage_results,function(x){x[1:which.max(x$cumulative_reward),]})
  final<-do.call(rbind,storage_results)[,c(1,4,3)]
  nodes_to_remove <- setdiff(setdiff(final$nextState,final$states),final_sample_summary[[sample]]$Clones$Clone)
  final <- final%>%filter(!nextState%in%nodes_to_remove)
  
  nodes_to_remove <- setdiff(setdiff(final$nextState,final$states),final_sample_summary[[sample]]$Clones$Clone)
  final <- final%>%filter(!nextState%in%nodes_to_remove)
  final[!duplicated(final),]
})

all_results_filtered<-plyr::compact(all_results_filtered)#[is.null(all_results_filtered)]<-NULL

all_results_filtered<-setNames(lapply(names(all_results_filtered),function(x){
  data.frame(all_results_filtered[[x]],
             "mutation"=x,
             "edge"=paste(all_results_filtered[[x]]$states,all_results_filtered[[x]]$nextState,sep="->"))
}),names(all_results_filtered))


      final<-do.call(rbind,all_results_filtered)
      if(length(final)==0){
        return("exclude")
      } else {
        final%>%group_by(mutation)%>%distinct(nextState,reward)%>%summarize(total=sum(reward))
      }
})

favorite_genes <- c("ASXL1","DNMT3A","TET2","IDH1","IDH2","NPM1","KRAS","NRAS","PTPN11","FLT3","JAK2")

out_mat<-data.frame(do.call(rbind,output))%>%filter(mutation!="exclude")
out_mat$Gene <- do.call(rbind,strsplit(out_mat$mutation,split="[_\\.]"))[,1]
out_mat<- out_mat%>%filter(Gene%in%favorite_genes)
out_mat$total <- as.numeric(out_mat$total)

gene_order <- out_mat%>%group_by(Gene)%>%summarize(median=median(total))%>%arrange(median)%>%pull(Gene)
out_mat$Gene <- factor(out_mat$Gene,levels=gene_order)
ggA<-ggplot(out_mat,aes(x=Gene,y=total,fill=Gene))+
                      geom_boxplot(outlier.shape=NA)+
                      scale_fill_manual(values=tol.rainbow(n=length(levels(out_mat$Gene))))+
                      geom_jitter(width=0.2,size=0.5)+  coord_flip()+theme_bw(base_size=8)+
                      guides(fill=FALSE)+
                      ylab("Fraction of sample explained \n by initiating mutation")
                      

ggsave(ggA, width=2.5,height=5,
       file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/Figures/F3B-initiating_mutation_contribution.pdf")



