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
}),names(final_results))

output<-plyr::compact(output)

out_mat<-do.call(rbind,lapply(names(output),function(sample_name){
  sample <- output[[sample_name]]
  initial_state <- sample$states[1]
  sample$Sample <- sample_name
  sample%>%filter(states==initial_state)#%>%pull(nextState)
}))
out_mat$Gene <- do.call(rbind,strsplit(out_mat$mutation,split="[_\\.]"))[,1]
favorite_genes <- c("ASXL1","DNMT3A","TET2","IDH1","IDH2","NPM1","KRAS","NRAS","PTPN11","FLT3","JAK2")

out_mat<- out_mat%>%filter(Gene%in%favorite_genes)
out_mat$present <-ifelse(out_mat$reward>0,1,0)

gene_order <- out_mat%>%group_by(Gene)%>%summarize(Count=sum(present))%>%arrange(Count)%>%pull(Gene)
out_mat$Gene <- factor(out_mat$Gene,levels=rev(gene_order))

table(do.call(c,lapply(split(out_mat,f=out_mat$Sample),function(x){
  sum(x$present)
})))


ggA<-ggplot(out_mat,aes(x=Gene,fill=factor(present)))+
  geom_bar(stat="count")+
  ylab("Number of samples")+
  xlab("")+
 # ylab("Fraction of samples with \n homozygous mutant in >10% of cells")+
  theme_classic(base_size=9)+
  scale_y_continuous(expand=c(0,0))+
  #theme(axis.text.x = element_text(angle=30,hjust=1))+
  scale_fill_manual(values=c("grey80","dodgerblue4"),guide=FALSE)


#ggsave(ggA, width=4,height=2,
 #      file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/Figures/F3C-initiating_mutation_monoallelic_presence.pdf")



out_mat$Gene <- factor(out_mat$Gene,levels=rev(gene_order))


DNMT3A<-data.frame("R882"=tally(out_mat%>%filter(Gene=="DNMT3A")%>%filter(grepl("DNMT3A.p.R882",mutation))%>%group_by(present)%>%select(Sample))%>%pull(n),
                  "Other"=tally(out_mat%>%filter(Gene=="DNMT3A")%>%filter(!grepl("DNMT3A.p.R882",mutation))%>%group_by(present)%>%select(Sample))%>%pull(n))
DNMT3A$Mono_allele <- c("Absent","Present")

ggB<-ggplot(melt(DNMT3A),aes(x=factor(variable,levels=c("Other","R882")),y=value,fill=Mono_allele))+
                    geom_bar(stat="identity")+theme_classic(base_size=9)+
                   # theme(axis.text.x = element_text(angle=30,hjust=1))+
                    xlab("")+ylab("Sample count")+
                    scale_y_continuous(expand=c(0,0))+
                    scale_fill_manual(values=c("Absent"="grey80","Present"=brewer.pal(5,"Reds")[5]),guide=FALSE)

#ggsave(ggB, width=2,height=2,
#       file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/Figures/F3D-DNMT3A_hotspot_vs_other.pdf")


plot_grid(ggA,ggB, align="hv",axis="tb",rel_widths = c(1,0.3))

ggsave(plot_grid(ggA,ggB, align="hv",axis="tb",rel_widths = c(1,0.3)), width=5.75,height=2,
       file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/Figures/F3CD-initiating_mutation_monoallelic_presence_DNMT3A.pdf")

