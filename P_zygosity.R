source("/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Scripts/mBio_processing_scripts.R")
library(UpSetR)

final_sample_summary<-readRDS(file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/final_sample_summary.rds")
clonal_sample_set_after_boostrap<-readRDS(file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/clonal_sample_set_after_boostrap")
genes <- c("ASXL1","DNMT3A","TET2","IDH1","IDH2","NPM1","FLT3","KRAS","NRAS","PTPN11","JAK2")

genes_in_each_sample<-setNames(mclapply(mc.cores = 4,final_sample_summary[clonal_sample_set_after_boostrap],function(x){
  aggregate_by<-do.call(rbind,strsplit(colnames(x$NGT),split="[_\\.]")) [,1][-length(colnames(x$NGT))]}), clonal_sample_set_after_boostrap)

list_of_mutants<-setNames(lapply(genes,function(x){
  names(genes_in_each_sample)[do.call(c,setNames(lapply(genes_in_each_sample,function(y){x%in%y }),names(genes_in_each_sample)))]}),genes)

data<-fromList(list_of_mutants)
rownames(data)<-unique(unlist(list_of_mutants))
data$Sample<-unique(unlist(list_of_mutants))
data$Group <- ifelse(grepl("CH",rownames(data)),"CH","AML")
data[data$Sample%in%c("MA1715B","MA6300A","MA9521A","MA9521B"),"Group"]<- "MPN"
data[data$Sample%in%c("E4840","E4838new"),"Group"]<- "CMML"
data[data$Sample%in%c("MA4244A","MA2725","MA9521B","MA0092A","R2715","MA1715B","MA6363B"),"Group"]<- "MF"

AML_samples <- data%>%filter(Group=="AML")%>%pull(Sample)
AML_samples_summary <-final_sample_summary[AML_samples]



merged_by_gene<-do.call(rbind,setNames(lapply(AML_samples,function(y){
  sample<-AML_samples_summary[[y]]
  genes <- do.call(rbind,strsplit(colnames(sample$NGT),split="[_\\.]"))[,1]
  genes<- genes[1:(length(genes)-1)]
  if(length(unique(genes))>1){
  weights <- sample$Clones$Count/sum(sample$Clones$Count)
  names(weights) <- sample$Clones$Clone
  mut_mat <- as.matrix(do.call(rbind,strsplit(sample$Clones$Clone,split="_")))
  mode(mut_mat) <-"numeric"
  
  test<-t(aggregate(t(mut_mat),by=list(genes),FUN="sum"))
  colnames(test) <- test[1,]
  final<-data.frame(test[-1,])
  rownames(final)<-sample$Clones$Clone
  final[final>2] <-2
  
  melt(data.frame(do.call(rbind,setNames(apply(final,2,function(x){
    data.frame("WT"=sum(0+weights[x==0]),
                "Het"=sum(0+weights[x==1]),
                "Hom"=sum(0+weights[x==2]))
  }) ,colnames(final) ) ),"Gene"=colnames(final) ))
  }
}),AML_samples))

genes_of_interest <- c("ASXL1","DNMT3A","TET2","IDH1","IDH2","NPM1","FLT3","KRAS","NRAS","PTPN11","JAK2")



merged_by_gene$Gene<-factor(merged_by_gene$Gene,levels=merged_by_gene%>%filter(variable=="Hom")%>%group_by(Gene)%>%summarize(median=median(value))%>%arrange(median)%>%pull(Gene))

homozygote_table<-setNames(data.frame(inner_join(tally(merged_by_gene%>%filter(Gene%in%genes_of_interest)%>%filter(variable=="Hom")%>%group_by(Gene)),
                        tally(merged_by_gene%>%filter(Gene%in%genes_of_interest)%>%filter(variable=="Hom"&value>0.1)%>%group_by(Gene)),by="Gene")),c("Gene","Total","Homozygote"))

homozygote_table$Not_homozygose <- homozygote_table$Total-homozygote_table$Homozygote

as.character(homozygote_table$Gene)[order(homozygote_table$Homozygote/homozygote_table$Total)]

final <- melt(homozygote_table[,-2])
final$variable <-factor(final$variable,levels=c("Not_homozygose","Homozygote"))
final$Gene <- factor(final$Gene,levels=rev(as.character(homozygote_table$Gene)[order(homozygote_table$Homozygote/homozygote_table$Total)]))

zygosity<-ggplot(final,aes(x=Gene,y=value,fill=variable))+
                      geom_bar(stat="identity",position="fill")+
                      xlab("Gene")+
                      ylab("Fraction of samples with \n homozygous mutant in >10% of cells")+
                      theme_classic(base_size=10)+
                      scale_y_continuous(expand=c(0,0))+
                      scale_fill_manual(values=c("grey80","dodgerblue4"),guide=FALSE)
                      
ggsave(zygosity, width=6,height=3,
       file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/Figures/EF4A-zygosity.pdf")

