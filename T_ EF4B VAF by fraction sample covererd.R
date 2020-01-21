source("/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Scripts/mBio_processing_scripts.R")

final_sample_summary<-readRDS(file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/final_sample_summary.rds")
clonal_sample_set_after_boostrap<-readRDS(file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/clonal_sample_set_after_boostrap")

graph_results<-readRDS(file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/MDP_results.rds")
final_results<-readRDS(file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/MDP_trajectory_for_each_gene.rds")


output<-setNames(lapply(names(final_results),function(sample){
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
}),names(final_results))

favorite_genes <- c("ASXL1","DNMT3A","TET2","IDH1","IDH2","NPM1","KRAS","NRAS","PTPN11","FLT3","JAK2")

output<-setNames(lapply(names(output),function(x){
  data.frame(data.frame(output[[x]]),"Sample"=x)
}),names(output))

out_vaf <-setNames(lapply(names(output),function(sample){
  if(output[[sample]][,1]=="exclude"){
    return(NULL) 
  } else {
    VAF <- data.frame("VAF"=colSums(as.matrix(final_sample_summary[[sample]]$NGT[,1:(ncol(final_sample_summary[[sample]]$NGT)-1)]))/(nrow(final_sample_summary[[sample]]$NGT)*2)*100,
                      "mutation"=colnames(final_sample_summary[[sample]]$NGT[,1:(ncol(final_sample_summary[[sample]]$NGT)-1)]))
  x <-inner_join(output[[sample]],VAF)
  return(x)
  }
}),names(output))
  
out_vaf<-plyr::compact(out_vaf)
out_vaf_mat <- do.call(rbind,out_vaf)
out_vaf_mat$Gene <- do.call(rbind,strsplit(out_vaf_mat$mutation,split="[_\\.]"))[,1]
favorite_genes <- c("ASXL1","DNMT3A","TET2","IDH1","IDH2","NPM1","KRAS","NRAS","PTPN11","FLT3","JAK2")

cor(out_vaf_mat$VAF,out_vaf_mat$total)

ggA<-ggplot(out_vaf_mat%>%filter(Gene%in%favorite_genes),aes(y=VAF/100,x=total/100,color=Gene))+
          geom_point()+theme_minimal(base_size=10)+
          ylab("Variant allele frequency")+
          xlab("Fraction of sample explained \n by initiating mutation")+
          scale_color_manual(values=tol.rainbow(n=length(favorite_genes)))
ggsave(ggA, width=4,height=4,
       file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/Figures/EF4B-VAF_by_fraction_sample_explained.pdf")




DTAI_sample <-data.frame(out_vaf_mat%>%group_by(Sample)%>%filter(Gene%in%c("DNMT3A","IDH1","IDH2","TET2","ASXL1")))
sort(table(do.call(rbind,(lapply(split(DTAI_sample,f=DTAI_sample$Sample),function(x){paste(unique(x$Gene),sep="_",collapse="_")})))[,1]))

gene_combos<-data.frame("gene_combo"=do.call(rbind,(lapply(split(DTAI_sample,f=DTAI_sample$Sample),function(x){paste(unique(x$Gene),sep="_",collapse="_")}))))
gene_combos$Sample <-rownames(gene_combos)

set<-inner_join(DTAI_sample,gene_combos)


ggplot(set,aes(x=Gene,y=total,group=Sample))+geom_point()+geom_path()+facet_grid(.~gene_combo,scale="free_x")

ggA<-ggplot(test_DI2,aes(x=Gene,y=VAF,group=Sample))+geom_point()+geom_path()
ggB<-ggplot(test_DI1,aes(x=Gene,y=VAF,group=Sample))+geom_point()+geom_path()
ggC<-ggplot(test_DT,aes(x=Gene,y=VAF,group=Sample))+geom_point()+geom_path()
ggD<-ggplot(test_DI2,aes(x=Gene,y=total,group=Sample))+geom_point()+geom_path()
ggE<-ggplot(test_DI1,aes(x=Gene,y=total,group=Sample))+geom_point()+geom_path()
ggF<-ggplot(test_DT,aes(x=Gene,y=total,group=Sample))+geom_point()+geom_path()

plot_grid(ggD,ggE,ggF,ncol=3)
