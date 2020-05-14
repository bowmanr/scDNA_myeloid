source("/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Scripts/mBio_processing_scripts.R")

pheno_mut_melted<-readRDS("/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/pheno_mut_melted.rds")
final_sample_set<-readRDS("/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/final_sample_set.rds")


final_mut_melt <-pheno_mut_melted%>%filter(Sample.ID%in%final_sample_set)
final_mut_melt$Gene<- factor(final_mut_melt$Gene,levels=names(sort(table(final_mut_melt$Gene), decreasing=TRUE)))

## Melt to get a tally of how many mutations per patient
melted_mut_mat <- data.frame(count(final_mut_melt, Gene, Sample.ID))
melted_mut_mat$Gene<- factor(melted_mut_mat$Gene,levels=names(sort(table(melted_mut_mat$Gene),decreasing=TRUE)))

### create matrix for oncoprint
mut_mat <- table(melted_mut_mat$Sample.ID,melted_mut_mat$Gene)
mut_mat <- mut_mat[,colSums(mut_mat)!=0]

### Prepare matrix for co occurence map
cooccur_mat <- cooccur(mat=t(mut_mat),type="spp_site",only_effects = FALSE,eff_matrix=TRUE,thresh=FALSE,eff_standard=FALSE,spp_names=TRUE)$results
cooccur_mat$score<-  ifelse(cooccur_mat[,"p_lt"]<=0.05,-1,ifelse(cooccur_mat[,"p_gt"]<=0.05,1,0))
cooccur_mat_subset <- rbind(cooccur_mat[,c("sp1_name","sp2_name","score")],c(setdiff(cooccur_mat$sp2_name,cooccur_mat$sp1_name),setdiff(cooccur_mat$sp1_name,cooccur_mat$sp2_name),0))
cooccur_data_mat <- dcast(cooccur_mat_subset, sp1_name  ~ sp2_name, value.var="score")
rownames(cooccur_data_mat) <- cooccur_data_mat[,1]
cooccur_data_mat2<-cooccur_data_mat[,-1]
cooccur_data_mat2[is.na(cooccur_data_mat2)]<-0
cooccur_data_mat2[lower.tri(cooccur_data_mat2)]<-NA
cooccur_data_mat2$Gene <- rownames(cooccur_data_mat2)
cooccur_data_mat_melt <- na.omit(melt(cooccur_data_mat2, 'Gene', variable_name='Gene2'))
cooccur_data_mat_melt$variable <- factor(cooccur_data_mat_melt$variable, levels=rev(levels(cooccur_data_mat_melt$variable)))


# Triangle heatmap to compare cohorts
grob_corrplot<-ggplot(cooccur_data_mat_melt, aes(Gene, variable))+
  geom_tile(aes(fill = factor(value,levels=c(-1,1,0))), color='grey90') +
  scale_fill_manual(values=c("-1"="firebrick3","1"="steelblue2","0"="white"),"Correlation",
                    labels=c("Mutually Exclusive","Mutually Inclusive",""))+
  theme_classic(base_size=8)+xlab("")+ylab("")+
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1),
        axis.line = element_blank(),
        legend.position = c(0.8,1), 
        legend.justification = c(1, 1),
        legend.direction = "vertical")+
  theme(legend.key.size = unit(0.5,"line"))


ggsave(grob_corrplot,width=3.5,height=3.5,
       file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/Figures/EF1-cooccurrence_mat_bulk.pdf")
