source("/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Scripts/mBio_processing_scripts.R")


final_sample_summary<-readRDS(file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/final_sample_summary.rds")
clonal_sample_set_after_boostrap<-readRDS(file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/clonal_sample_set_after_boostrap")
sample_key<-read.delim(file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/From_LAM/OS_and_Cheatsheet/Mission_Bio_MRN_RunID_Samples.txt",sep="\t")
clinical_VAF<-read.delim(file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/From_LAM/Mission Bio_Clinical Variants_AZlist.txt",sep="\t")

clinical_VAF$Variant <- paste(clinical_VAF$Gene,clinical_VAF$AA,sep=".")

gene_level_info <- data.frame(do.call(rbind,lapply(names(final_sample_summary),function(y){
  x <- final_sample_summary[[y]]$NGT
  data.frame("Computed_VAF"=colSums(x[,!c(colnames(x)=="Clone")])/(length(rownames(x))*2),
             "Sample.Name"=y,
             "Gene"= do.call(rbind,strsplit(colnames(x)[1:(length(colnames(x))-1)],split="[_\\.]"))[,1],
             "Variant"=colnames(x)[1:(length(colnames(x))-1)])
})))

gene_level_info$Variant<-gsub("\\.$","*",gene_level_info$Variant)
gene_level_info$Variant<-gsub("fs\\.","fs*",gene_level_info$Variant)
gene_level_info$Variant<-gsub("\\.[1232456789]$","",gene_level_info$Variant)

complete_set<-inner_join(full_join(gene_level_info,sample_key),clinical_VAF,by=c("Molecular_accession_num"))

complete_set%>%Grou

set<-complete_set%>%filter(Variant.x==Variant.y | (Gene.x=="FLT3"&Gene.y=="FLT3"))                        
non_FLT3 <- set%>%filter(Gene.x!="FLT3")
FLT3 <- set%>%filter(Gene.x=="FLT3")

FLT3%>%group_by(Variant.x,Variant.y)%>%summarise(select=abs(as.numeric(Computed_VAF)-as.numeric(VAF)))
%>%filter(select==max(select))

FLT3$VAF_diff<-abs(FLT3$Computed_VAF-as.numeric(FLT3$VAF))

final<-rbind(data.frame(non_FLT3%>%group_by(Variant.x,Variant.y)),data.frame(FLT3%>%group_by(Variant.x,Variant.y)%>%filter(VAF_diff==min(VAF_diff))%>%select(-VAF_diff)))

cor.test(as.numeric(final$Computed_VAF),as.numeric(final$VAF),method="pearson")
plot(as.numeric(final$Computed_VAF),as.numeric(final$VAF))

final$Computed_VAF <- as.numeric(final$Computed_VAF)
final$VAF <- as.numeric(final$VAF)

ggCorrelation<-ggplot(final,aes(x=Computed_VAF,y=VAF,color=Gene.x))+
        geom_point(size=0.5)+xlab("Computed single cell VAF")+ylab("Bulk VAF")+theme_bw(base_size=10)#+
 # geom_smooth()

ggsave(addSmallLegend(ggCorrelation),width=5,height=4,
       file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/Figures/EF1E-correlation_plot.pdf")

