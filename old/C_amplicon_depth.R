source("/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Scripts/mBio_processing_scripts.R")

final_NGTs<-readRDS(file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/final_NGTs.rds")
pheno_mut_melted<-readRDS(file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/pheno_mut_melted.rds")
final_sample_set<-readRDS(file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/final_sample_set.rds")
sample_SNPS<-readRDS(file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/sample_SNPS.rds")
cells_of_interest<-readRDS(file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/cells_of_interest.rds")


sample_DPs<-setNames(lapply(as.list(names(sample_SNPS)),function(x){
  y<-extract_DP_files(x,sample_SNPS[[x]]$SNP_variant)
  rownames(y) <- paste("Cell",1:nrow(y),sep="_")
  return(y)
}),names(sample_SNPS))

#saveRDS(sample_DPs,file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/sample_DPs.rds")

sample<-1

processing_DP<-lapply(final_sample_set,function(sample){
  variants_of_interest <-rownames(sample_SNPS[[sample]])[sample_SNPS[[sample]]$protein%in%colnames(final_NGTs[[sample]])]
 if(length(variants_of_interest)>0){
   DP_subet<-  melt(as.matrix(sample_DPs[[sample]][cells_of_interest[[sample]],variants_of_interest]))
   output<-plyr::ddply(DP_subet, c( "Var2"), summarise,
                       mean = mean(value), sd = sd(value),
                       median = median(value), 
                       sem = sd(value)/sqrt(length(value)))
    } else{
   return(NULL)
 }
})

gg_coverage_per_amplicon <- data.frame("Coverage"=unlist(processing_DP))%>%ggplot(aes(x=log2(Coverage+1)))+
                                                                      geom_histogram()+
                                                                      theme_classic(base_size = 10)+
                                                                      ylab("Count")+
                                                                      ggtitle("Informative Reads per Amplicon")#+
                                                                     # scale_x_log10()
