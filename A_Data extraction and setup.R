source("/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Scripts/mBio_processing_scripts.R")

all_samples <- list.files("/Volumes/LevineLab/Levine Lab/MissionBio_Tapestri/Insights_Output/ASH_Cohort")
blacklist <- read.xls("/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/From_LAM/Blacklist2.xlsx")
whitelist <- read.xls("/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/From_LAM/Whitelist.xlsx")
pheno <- read.xls("/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/From_LAM/pheno_MBio.xlsx")
sample_set<- intersect(as.character(all_samples),as.character(pheno$Sample.ID)) 

setwd("/Volumes/LevineLab/Levine Lab/MissionBio_Tapestri/Insights_Output/ASH_Cohort")
sample_SNPS<-lapply(as.list(sample_set),function(x){
  y<-extract_SNP_info(x)
  z <- y[!(rownames(y)%in%blacklist[,"Blacklisted.coordinate"]|y[,"protein"]%in%blacklist[,"Protein.Change"]),]
  z[,"SNP_variant"] <- rownames(z)
  z[,"Sample"] <- x
  return(z)
})

names(sample_SNPS) <- sample_set
saveRDS(sample_SNPS,file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/sample_SNPS.rds")

sample_NGTs<-lapply(as.list(names(sample_SNPS)),function(x){
  extract_NGT_files(x,sample_SNPS[[x]]$SNP_variant)
})
names(sample_NGTs) <- names(sample_SNPS)

named_sample_NGTs<-lapply(sample_NGTs,function(x){
  rownames(x) <- paste("Cell",1:nrow(x),sep="_")
  return(x)
})
saveRDS(named_sample_NGTs,file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/sample_NGTs.rds")


pheno <- read.xls("/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/From_LAM/pheno_MBio.xlsx")

# Now we filter out any synonymous mutations as well
sample_SNPS_filter<-lapply(sample_SNPS,function(x){ x%>%filter(!(grepl("[=>]",protein)))%>%filter(!grepl("ASXL1:p.[DNGPIL]81",protein))%>%filter(!grepl("FLT3_INS_chr13:28602226:",protein)) }) 
colnames_to_grab <- c("Variant","protein","Gene","gene","allele_freq_loom","allele_freq","SNP_variant","Sample")
sample_SNPS_filter_cols<-lapply(sample_SNPS_filter,function(x){x[,colnames(x)%in%colnames_to_grab]}) 

sample_NGTs_filter<-lapply(names(sample_SNPS),function(x){setNames(data.frame(sample_NGTs[[x]][,sample_SNPS_filter[[x]][,"SNP_variant"]]),sample_SNPS_filter[[x]][,"SNP_variant"])}); names(sample_NGTs_filter) <- names(sample_SNPS)
final_variants_of_interest<-lapply(sample_NGTs_filter,function(x){colnames(x)[colSums(x)>0]})
names(final_variants_of_interest) <-names(sample_SNPS)

test_NGT <- setNames(lapply(names(final_variants_of_interest),function(x){setNames(data.frame(sample_NGTs_filter[[x]][,final_variants_of_interest[[x]]]),final_variants_of_interest[[x]] )}),names(sample_SNPS))

final_NGT<-setNames(lapply(names(test_NGT),function(x){
  if(ncol(test_NGT[[x]])>1){
    j<-test_NGT[[x]][,apply(test_NGT[[x]],2,function(z){sum(z%in%c(1,2))>=2})]
    z<-j[!apply(j,1,function(y){any(y==3)}),]
    q<-z[,apply(z,2,function(h){sum(h%in%c(1,2))>=2})]
  } else{
    j<-data.frame(test_NGT[[x]][test_NGT[[x]]==1|test_NGT[[x]]==2])
    q<-data.frame(j[!j==3])      
  }
  return(q)
}),names(test_NGT))

final_variants_of_interest <- lapply(final_NGT,colnames)

final_SNP <-setNames(lapply(names(final_variants_of_interest),function(x){
  sample_SNPS_filter_cols[[x]]%>%filter(SNP_variant%in%final_variants_of_interest[[x]])}),names(sample_SNPS))

final_SNP_mat<-do.call(bind_rows,final_SNP)
final_SNP_mat$Gene <- ifelse(final_SNP_mat$Gene=="","FLT3",final_SNP_mat$Gene)
final_SNP_mat$gene <- ifelse(final_SNP_mat$gene=="","FLT3",final_SNP_mat$gene)
final_SNP_mat$Gene <- ifelse(is.na(final_SNP_mat$Gene),final_SNP_mat$gene,final_SNP_mat$Gene) 
colnames(final_SNP_mat)[5] <- "Sample.ID"

pheno_mut_melted<-inner_join(pheno,final_SNP_mat[,2:5])
pheno_mut_melted$Gene<- factor(pheno_mut_melted$Gene,levels=names(sort(table(pheno_mut_melted$Gene), decreasing=TRUE)))

for(x in 1:length(names(final_NGT))){
  colnames(final_NGT[x][[1]])<-pheno_mut_melted[pheno_mut_melted$SNP_variant%in%colnames(final_NGT[[x]])&
                                                  pheno_mut_melted$Sample.ID==names(final_NGT)[x],"protein"]
  rownames(final_NGT[x][[1]])<-paste("Cell",1:dim(final_NGT[x][[1]])[[1]],sep="_")
}


#These lines were used to regenerated a new blacklist
protein_variant_counts<-data.frame(sort(table(pheno_mut_melted[,"protein"]),decreasing=TRUE))
SNP_variant_counts<-data.frame(sort(table(pheno_mut_melted[,"SNP_variant"]),decreasing=TRUE))
Muts_per_patientx<-data.frame(sort(table(pheno_mut_melted[,"Sample.ID"]),decreasing=TRUE))

processing_NGT <- final_NGT
names(processing_NGT) <- names(final_NGT)

cell_number_per_sample_cutoff<-40
final_sample_set<-names(which(do.call(rbind,lapply(processing_NGT,dim))[,1]>cell_number_per_sample_cutoff))

samples_of_interest<-setdiff(names(processing_NGT),final_sample_set)
lapply(processing_NGT[samples_of_interest],nrow)
lapply(sample_NGTs[samples_of_interest],nrow)

cells_of_interest<-setNames(lapply(processing_NGT[final_sample_set],function(x){
  rownames(x)
}),final_sample_set)

saveRDS(processing_NGT[final_sample_set],file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/final_NGTs.rds")
saveRDS(pheno_mut_melted,file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/pheno_mut_melted.rds")
saveRDS(final_sample_set,file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/final_sample_set.rds")
saveRDS(cells_of_interest,file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/cells_of_interest.rds")

processing_NGT<-final_NGTs

clonal_sample_set <- final_sample_set[do.call(rbind,lapply(final_NGTs,dim))[,2]>1]
NGT_to_clone<-lapply(processing_NGT[clonal_sample_set],function(y){
  bulk_VAF_order <-names(sort(colSums(y),decreasing=TRUE))
  x <- data.frame(y[,bulk_VAF_order],"Clone"=apply(y[,bulk_VAF_order],1,function(z){paste(z,sep="_",collapse="_")}))
})

clonal_abundance<- lapply(NGT_to_clone,function(x){
  y <- data.frame(data.table(data.frame("Count"=as.matrix(table(x[,"Clone"])),
                                        "Clone"=names(table(x[,"Clone"]))),
                             key="Count"))
  y$Clone <- factor(y$Clone,levels=rev(c(y$Clone)))
  return(y)
})
set.seed(68864)
clonal_abundance_boot_CI<-setNames(mclapply(mc.cores=8,names(NGT_to_clone),function(sample_to_test){
  bootstrap_NGT_mat(NGT_to_clone_list=NGT_to_clone,
                    clonal_abundance=clonal_abundance,
                    sample_to_test=sample_to_test,
                    replicate=1000,
                    clone_cutoff=10)}),names(NGT_to_clone))

oops<-data.frame("Clones"=do.call(rbind,lapply(clonal_abundance_boot_CI,dim))[,1],
                 "Cells"=do.call(rbind,lapply(clonal_abundance_boot_CI,function(x){sum(x[,1])})))
oops$Sample <- rownames(oops)
clonal_sample_set_final <- data.frame(oops%>%filter(Clones>1&Cells>100))$Sample

#clonal_abundance_ADO_CI<-setNames(mclapply(mc.cores=8,(clonal_sample_set_final),function(sample_to_test){
 # bootstrap_allele_dropout_NGT_mat(NGT_to_clone_list=NGT_to_clone,
  #                                 clonal_abundance=clonal_abundance,
   #                                sample_to_test=sample_to_test,
    #                               replicate=100,
     #                              clone_cutoff=10)}),(clonal_sample_set_final))

#oops2<-data.frame("Clones"=do.call(rbind,lapply(clonal_abundance_ADO_CI,dim))[,1],
  #                "Cells"=do.call(rbind,lapply(clonal_abundance_ADO_CI,function(x){sum(x[,1])})))
#oops2$Sample <- clonal_sample_set_final
#clonal_sample_set_final_ADO <- data.frame(oops2%>%filter(Clones>1&Cells>100))$Sample

#names(clonal_abundance_ADO_CI)<-clonal_sample_set_final
#clonal_abundance_random_CI<-setNames(mclapply(mc.cores=2,names(NGT_to_clone),function(sample_to_test){
#                                             random_distribution_NGT_CI(NGT_to_clone_list=NGT_to_clone,
#                                                              clonal_abundance=clonal_abundance,
#                                                              sample_to_test=sample_to_test,
#                                                              replicate=10)}),names(NGT_to_clone))

harmonized_clonal_abundace <- setNames(mclapply(mc.cores=2,clonal_sample_set_final,function(sample_to_test){
  set<- clonal_abundance_boot_CI[[sample_to_test]]
  # set <- clonal_abundance[[sample_to_test]]%>%filter(Count>=5)
  # set$Clone <- as.character(set$Clone)
 # set<-inner_join(clonal_abundance_boot_CI[[sample_to_test]],
 #                 clonal_abundance_ADO_CI[[sample_to_test]],by="Clone")[,-5]
  colnames(set)[1]<-"Count"
  return(set)}),clonal_sample_set_final)



final_sample_summary <- setNames(lapply(names(harmonized_clonal_abundace),function(sample_to_test){
  print(sample_to_test)
  if(nrow(harmonized_clonal_abundace[[sample_to_test]])==0) {return("No clones after boostrapping")}
  dedup_NGT<-as.matrix(do.call(rbind,strsplit(as.character(harmonized_clonal_abundace[[sample_to_test]][,"Clone"]),split="_")))
  mode(dedup_NGT) <- "numeric"
  cols_to_remove<-which(colSums(dedup_NGT)==0)
  print(cols_to_remove)      
  if(nrow(dedup_NGT)==1) {
    return("Only 1 clone left")
  }else  if(2%in%c(cols_to_remove)){
    return("Removed all but 1 variant")
  }else if(length(cols_to_remove)>=1){
    NGT_to_clone_subset <- NGT_to_clone[[sample_to_test]]%>%filter(as.character(Clone)%in%as.character(harmonized_clonal_abundace[[sample_to_test]]$Clone))
    NGT_to_clone_subset <- NGT_to_clone_subset[,-cols_to_remove]
    NGT_to_clone_subset$Clone <- apply(NGT_to_clone_subset[,-ncol(NGT_to_clone_subset)],1,function(x){paste(x,sep="_",collapse="_")})
    harmonized_clonal_abundace[[sample_to_test]]$Clone<-apply(do.call(rbind,strsplit(as.character(harmonized_clonal_abundace[[sample_to_test]]$Clone),split="_") )[,-c(cols_to_remove)],1,paste,sep="_",collapse="_")
    
  } else{
    
    NGT_to_clone_subset <- NGT_to_clone[[sample_to_test]]%>%filter(as.character(Clone)%in%as.character(harmonized_clonal_abundace[[sample_to_test]]$Clone))
  }
  clonal_architecture<-melt(NGT_to_clone_subset[!duplicated(NGT_to_clone_subset),])
  colnames(clonal_architecture) <- c("Clone","Mutant","Genotype")
  clonal_architecture$Genotype <- factor(with(clonal_architecture,ifelse(Genotype==3,NA,
                                                                         ifelse(Genotype==0,"WT",
                                                                                ifelse(Genotype==1,"Heterozygous",
                                                                                       ifelse(Genotype==2,"Homozygous",Genotype))))),
                                         levels=c("WT","Heterozygous","Homozygous"))
  return(list("Clones"=harmonized_clonal_abundace[[sample_to_test]],
              "NGT"=NGT_to_clone_subset,
              "Architecture"=clonal_architecture))
}),names(harmonized_clonal_abundace))

final_sample_summary[which(do.call(rbind,lapply(final_sample_summary,length))==1 )]<-NULL

clonal_sample_set_after_boostrap <-names(final_sample_summary)

removed_genes_colSums <- setNames(lapply((clonal_sample_set_after_boostrap),function(sample_to_test){
  dedup_NGT<-as.matrix(do.call(rbind,strsplit(as.character(harmonized_clonal_abundace[[sample_to_test]][,"Clone"]),split="_")))
  mode(dedup_NGT) <- "numeric"
  cols_to_remove<-which(colSums(dedup_NGT)==0)
  if(length(cols_to_remove)>1&sample_to_test!="MA4849B"){
    colSums(NGT_to_clone[[sample_to_test]][,cols_to_remove])
  }}),(clonal_sample_set_after_boostrap))


saveRDS(final_sample_summary,file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/final_sample_summary.rds")
saveRDS(clonal_sample_set_after_boostrap,file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/clonal_sample_set_after_boostrap")
