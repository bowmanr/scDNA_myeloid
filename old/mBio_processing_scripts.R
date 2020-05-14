library(dplyr)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(grid)
library(cowplot)
library(RColorBrewer)
library(data.table)
library(ggbeeswarm)
library(purrr)
library(vegan)
library(igraph)
library(reshape2)
library(gdata)
library(TRONCO)
library(qgraph)
library(ggnet)
library(network)
library(sna)
library(ggplot2)
library(tidyverse)
library(tidygraph)
library(ReinforcementLearning)
library(pheatmap)
library(ggplotify)
library(cooccur)
library(knitr)
library(forcats)
library(tidyr)
library(corrplot)
library(ggcorrplot)
library(dplyr)
library(ROI)
library(ROI.plugin.glpk)
library(ompr)
library(ompr.roi)
library(ape)
library(phylobase)
library(ggtree)
library(miscTools)
library(parallel)
library(ggridges)
library(pals)
library(widyr)
library(circlize)
library(rcompanion)
library(ggrepel)
library(ggthemes)
library(ggnet)
library(network)
library(sna)

options(stringsAsFactors=FALSE)

random_distribution_NGT_CI <- function(NGT_to_clone_list,clonal_abundance,sample_to_test,replicates,clone_cutoff){
                        mat_of_interest <- NGT_to_clone[[sample_to_test]][,1:(ncol(NGT_to_clone[[sample_to_test]])-1)]
                        bulk_VAF_order <-names(sort(colSums(mat_of_interest),decreasing=TRUE))
                        test<-replicate(n=replicates,apply(mat_of_interest,MARGIN=c(2),function(x){
                          z<-sample(x,size=length(x),replace=FALSE)
                          return(as.numeric(z))}),simplify = "array")
                        test2<-apply(test,3,function(q){
                          x<-data.frame(q[,bulk_VAF_order],
                                        "Clone"=apply(q[,bulk_VAF_order],1, function(p){paste(p,sep="_",collapse="_")}))
                          y<-data.table(data.frame("Count"=as.matrix(table(x[,"Clone"])),
                                                   "Clone"=names(table(x[,"Clone"]))),
                                        key="Count")
                          y$Clone<-factor(y$Clone,levels=rev(as.character(y$Clone)))
                          return(data.frame(y))})
                        if(class(test2)=="list"){y <- lapply(test2,data.frame) %>% reduce(full_join, by = "Clone")}
                        if(class(test2)=="array"){y <- apply(test2,3,data.frame) %>% reduce(full_join, by = "Clone")}
                        y[is.na(y)]<-0
                        rownames(y)<-y[,"Clone"]
                        y<-as.matrix(y[,-2])
                        mode(y)<-'numeric'
                        z<-t(apply(y,1,function(p){quantile(p,probs=c(0.025,0.975))}))
                        z_mat<-data.frame(z,"Clone"=rownames(z))
                        set<-setNames(data.frame(full_join(data.frame(clonal_abundance[[sample_to_test]]),z_mat,by="Clone")),c("Count","Clone","LCI","UCI"))%>%filter(LCI>=clone_cutoff&!is.na(Count))
                        return(set)
}

bootstrap_allele_dropout_NGT_mat <- function(NGT_to_clone_list,clonal_abundance,sample_to_test,replicates,clone_cutoff){
                  mat_of_interest <- NGT_to_clone[[sample_to_test]][,1:(ncol(NGT_to_clone[[sample_to_test]])-1)]
                  bulk_VAF_order <-names(sort(colSums(mat_of_interest),decreasing=TRUE))
                  test<-replicate(n=replicates,apply(mat_of_interest,MARGIN=c(1,2),function(x){
                    if(x==1){z<-x} 
                    if(x==0){z<-sample(x=c(0,1),size=1,prob=c(0.9,0.1))}
                    if(x==2){z<-sample(x=c(2,1),size=1,prob=c(0.9,0.1))}
                    return(as.numeric(z))}),simplify = "array")
                  test2<-apply(test,3,function(q){
                    x<-data.frame(q[,bulk_VAF_order],
                                  "Clone"=apply(q[,bulk_VAF_order],1, function(p){paste(p,sep="_",collapse="_")}))
                    y<-data.table(data.frame("Count"=as.matrix(table(x[,"Clone"])),
                                             "Clone"=names(table(x[,"Clone"]))),
                                  key="Count")
                    y$Clone<-factor(y$Clone,levels=rev(as.character(y$Clone)))
                    return(data.frame(y))})
                  
                  if(class(test2)=="list"){y <- lapply(test2,data.frame) %>% reduce(full_join, by = "Clone")}
                  if(class(test2)=="array"){y <- apply(test2,3,data.frame) %>% reduce(full_join, by = "Clone")}
                  y[is.na(y)]<-0
                  rownames(y)<-y[,"Clone"]
                  y<-as.matrix(y[,-2])
                  mode(y)<-'numeric'
                  z<-t(apply(y,1,function(p){quantile(p,probs=c(0.025,0.975))}))
                  z_mat<-data.frame(z,"Clone"=rownames(z))
                  set<-setNames(data.frame(full_join(data.frame(clonal_abundance[[sample_to_test]]),z_mat,by="Clone")),c("Count","Clone","LCI","UCI"))%>%filter(LCI>=clone_cutoff&!is.na(Count))
                  return(set)
} 

bootstrap_NGT_mat <- function(NGT_to_clone_list,clonal_abundance,sample_to_test,replicates,clone_cutoff){
  test<-replicate(n=replicates,resample_fun(NGT_to_clone[[sample_to_test]]),simplify = "array")
  if(class(test)=="list"){
    y <- lapply(test,data.frame) %>% reduce(full_join, by = "Clone")
  }
  if(class(test)=="array"){
    y <- apply(test,3,data.frame) %>% reduce(full_join, by = "Clone")
  }
  y[is.na(y)]<-0
  rownames(y)<-y[,"Clone"]
  y<-as.matrix(y[,-2])
  mode(y)<-'numeric'
  z<-t(apply(y,1,function(p){
    quantile(p,probs=c(0.025,0.975))
  }))
  z_mat<-data.frame(z,"Clone"=rownames(z))
  set<-setNames(data.frame(inner_join(data.frame(clonal_abundance[[sample_to_test]]),z_mat,by="Clone")),c("Count","Clone","LCI","UCI"))%>%filter(LCI>=clone_cutoff)
  return(set)
}

resample_fun<-function(data){
  x <- data[sample(x=1:nrow(data),replace=TRUE),]
  return(as.matrix(data.table(data.frame("Count"= as.matrix(table(x[,"Clone"])),
                                         "Clone"=names(table(x[,"Clone"]))),
                              key="Count")))
}

extract_NGT_files <- function(sample,variants){
  NGT <- read.csv(paste("./",sample,"/","NGT.csv",sep="")) # - genotype call converted to categorical (0-reference, 1-heterozygous mutation, 2-homozygous mutation, 3-unknown)
 # cells_with_unknown<-NGT[apply(NGT[,variants],MARGIN=1,FUN=function(x){sum(x==3)>0}),"Cell"]
#  matrix_of_interest<-NGT[!NGT$Cell%in%cells_with_unknown,variants]
  #matrix_of_interest<- matrix_of_interest[,colSums(matrix_of_interest)>0]
  return(NGT)
}

count_unknown_cells<- function(sample,variants){
    NGT <- read.csv(paste("./",sample,"/","NGT.csv",sep="")) # - genotype call converted to categorical (0-reference, 1-heterozygous mutation, 2-homozygous mutation, 3-unknown)
    cells_with_unknown<-NGT[apply(NGT[,variants],MARGIN=1,FUN=function(x){sum(x==3)>0}),"Cell"]

    return(c(length(cells_with_unknown)/nrow(NGT)))
  }

extract_DP_files <- function(sample,variants){
  DP <- read.csv(paste("./",sample,"/","DP.csv",sep="")) # - genotype call converted to categorical (0-reference, 1-heterozygous mutation, 2-homozygous mutation, 3-unknown)
  cells_with_unknown<-DP[apply(DP[,variants],MARGIN=1,FUN=function(x){sum(x==3)>0}),"Cell"]
  matrix_of_interest<-DP[!DP$Cell%in%cells_with_unknown,variants]
  #matrix_of_interest<- matrix_of_interest[,colSums(matrix_of_interest)>0]
  return(matrix_of_interest)
}

extract_SNP_info <- function(sample){
  AF <- read.csv(paste("./",sample,"/","AF.csv",sep=""))                    # - variant allele frequency
  
  if(any(grepl("SNP_INFO.csv",list.files(paste("./",sample,"/",sep=""))))){
    SNP_INFO <- read.csv(paste("./",sample,"/","SNP_INFO.csv",sep="")) 
    rownames(SNP_INFO) <-colnames(AF)[-c(1:2)]  # - renames SNP info to have the variants as rownames
    SNP_INFO[,c("protein")] <- apply(SNP_INFO,1,function(x){
      ifelse(x["protein"]!="",x["protein"],
             ifelse(x["cdna"]!="",paste(x["gene"],x["cdna"],sep="_"),
                    ifelse(x["amplicon"]%in%c("MSK_RL_AMP54", "MSK_RL_AMP55", "MSK_RL_AMP56","MSK_RL_AMP57","MSK_RL_AMP58"),paste("FLT3_INS",x["POS"],x["ALT"],sep="_"),"")))
    })# - variant metadata
  } else {
    SNP_INFO <- read.csv(paste("./",sample,"/","Variants.csv",sep=""))  
    rownames(SNP_INFO) <-colnames(AF)[-c(1:2)]  # - renames SNP info to have the variants as rownames
    colnames(SNP_INFO)[3] <- "protein"
    SNP_INFO$cDNA <- as.character(SNP_INFO$cDNA)
    SNP_INFO$protein <- as.character(SNP_INFO$protein)
    SNP_INFO$Variant <- as.character(SNP_INFO$Variant)
    SNP_INFO$Gene <- as.character(SNP_INFO$Gene)
    SNP_INFO[,c("protein")] <- apply(SNP_INFO,1,function(x){
      ifelse(x["protein"]!="",x["protein"],
             ifelse(x["cDNA"]!="",paste(x["Gene"],x["cDNA"],sep="_"),
                    ifelse(grepl("^chr13",x["Variant"]),paste("FLT3_INS",x["Variant"],sep="_"),"")))
     }) # - variant metadata
  }
  return(SNP_INFO)
}

# run the following just once to load the function
plot_variants_and_selection <- function(sample){
  #Load in the data
  AF <- read.csv(paste("./",sample,"/","AF.csv",sep=""))                    # - variant allele frequency
  DP <- read.csv(paste("./",sample,"/","DP.csv",sep=""))                     # - read depth
  GQ <- read.csv(paste("./",sample,"/","GQ.csv",sep=""))                    # - quality scores
  NGT <- read.csv(paste("./",sample,"/","NGT.csv",sep="")) # - genotype call converted to categorical (0-reference, 1-heterozygous mutation, 2-homozygous mutation, 3-unknown)
  
  if(any(grepl("SNP_INFO.csv",list.files(paste("./",sample,"/",sep=""))))){
    SNP_INFO <- read.csv(paste("./",sample,"/","SNP_INFO.csv",sep="")) 
    rownames(SNP_INFO) <-colnames(AF)[-c(1:2)]  # - renames SNP info to have the variants as rownames
    SNP_INFO[,c("protein")] <- apply(SNP_INFO,1,function(x){
      ifelse(x["protein"]!="",x["protein"],
             ifelse(x["cdna"]!="",paste(x["gene"],x["cdna"],sep="_"),
                    ifelse(x["amplicon"]%in%c("MSK_RL_AMP54", "MSK_RL_AMP55", "MSK_RL_AMP56","MSK_RL_AMP57","MSK_RL_AMP58"),paste("FLT3_INS",x["POS"],x["ALT"],sep="_"),"")))
    })# - variant metadata
  } else {
    SNP_INFO <- read.csv(paste("./",sample,"/","Variants.csv",sep=""))  
    rownames(SNP_INFO) <-colnames(AF)[-c(1:2)]  # - renames SNP info to have the variants as rownames
    colnames(SNP_INFO)[3] <- "protein"
    SNP_INFO$cDNA <- as.character(SNP_INFO$cDNA)
    SNP_INFO$protein <- as.character(SNP_INFO$protein)
    SNP_INFO$Variant <- as.character(SNP_INFO$Variant)
    SNP_INFO$Gene <- as.character(SNP_INFO$Gene)
    SNP_INFO[,c("protein")] <- apply(SNP_INFO,1,function(x){
      ifelse(x["protein"]!="",x["protein"],
             ifelse(x["cDNA"]!="",paste(x["Gene"],x["cDNA"],sep="_"),
                    ifelse(grepl("^chr13",x["Variant"]),paste("FLT3_INS",x["Variant"],sep="_"),"")))
    })# - variant metadata
  }
  
  protein_changes_of_interest <- (SNP_INFO %>% filter(protein!="")%>%select(protein))[,1] # select column of variants
  SNP_changes_of_interest <- rownames(SNP_INFO)[SNP_INFO[,"protein"]%in%protein_changes_of_interest]
  SNP_protein_key <- data.frame("Mutant"=SNP_changes_of_interest,"Protein"=protein_changes_of_interest)
  
  # Cell and clone selection and Average depth per cell per allele
  distribution_of_unknowns_by_variant <- apply(NGT[,SNP_changes_of_interest],MARGIN=2,table)
  NGT_subset<-NGT[,c("Sample","Cell",SNP_changes_of_interest)]
  #NGT_interest$Clone <- apply(NGT_interest[,-c(1,2)],1,function(x){paste(x,sep="_",collapse="_")})
  #clonal_abundance <- sort(table(NGT_interest$Clone))
  
  DP_subset<- DP[,c("Sample","Cell",SNP_changes_of_interest)]
  DP_NGT_melt<-inner_join(melt(DP_subset[,-1],id.var="Cell"),melt(NGT_subset[,-1],id.var="Cell"),by=c("Cell","variable"))
  colnames(DP_NGT_melt) <- c("Cell","Mutant","Depth","Genotype")
  DP_NGT_melt$Genotype <- ifelse(DP_NGT_melt$Genotype==3,"Unknown",
                                 ifelse(DP_NGT_melt$Genotype==0,"WT",
                                        ifelse(DP_NGT_melt$Genotype==1,"Heterozygous",
                                               ifelse(DP_NGT_melt$Genotype==2,"Homozygous",DP_NGT_melt$Genotype))))
  DP_NGT_melt$Genotype  <- factor(DP_NGT_melt$Genotype ,levels=c("WT","Heterozygous","Homozygous","Unknown"))
  
  summarized<-data.frame(DP_NGT_melt%>%group_by(Mutant,Genotype)%>%summarise("Depth"=mean(Depth)))
  summarized_no_unknown <- summarized %>% filter(Genotype!="Unknown")
  grouped <- group_by(summarized_no_unknown, Mutant)
  output <- data.frame(summarise(grouped, mean=mean(Depth), sd=sd(Depth))) 
  output$Include <- rep("Include",length(rownames(output)))
  output_final <- inner_join(SNP_protein_key,output)
  
  depth_plot<- ggplot(summarized_no_unknown,aes(x=Mutant,y=Depth,color=Genotype))+geom_point()+coord_flip()+
    scale_x_discrete(limits = rev(levels(summarized_no_unknown$Mutant)))+theme_bw()#+
  
  write.csv(output_final,paste0("./",sample,"/",sample,"_VARIANT_SELECTION.csv"),row.names=FALSE)
  ggsave(paste0("./",sample,"/",sample,"_VARIANT_SELECTION.pdf"),depth_plot,width=10,height=20)
}

run_mbio_analysis<-function(sample,diversity){
  #Load in the data
  
  if(any(grepl("VARIANT_SELECTION_LAM",list.files(paste("./",sample,"/",sep=""))))){
    LAM_VARIANTS_MAT <- read.csv(paste("./",sample,"/",sample,"_VARIANT_SELECTION_LAM.csv",sep=""))  
    LAM_VARIANTS_SNP <- as.character((LAM_VARIANTS_MAT %>% filter(Include=="Include")%>%select(Mutant))[,1])
    LAM_VARIANTS_Protein <- as.character((LAM_VARIANTS_MAT %>% filter(Include=="Include")%>%select(Protein))[,1])
  } else {
    print("No variant selection performed!") 
    return("next")
  }
  AF <- read.csv(paste("./",sample,"/","AF.csv",sep=""))                    # - variant allele frequency
  DP <- read.csv(paste("./",sample,"/","DP.csv",sep=""))                     # - read depth
  GQ <- read.csv(paste("./",sample,"/","GQ.csv",sep=""))                    # - quality scores
  NGT <- read.csv(paste("./",sample,"/","NGT.csv",sep="")) # - genotype call converted to categorical (0-reference, 1-heterozygous mutation, 2-homozygous mutation, 3-unknown)
  
  if(any(grepl("SNP_INFO.csv",list.files(paste("./",sample,"/",sep=""))))){
    SNP_INFO <- read.csv(paste("./",sample,"/","SNP_INFO.csv",sep="")) 
    rownames(SNP_INFO) <-colnames(AF)[-c(1:2)]  # - renames SNP info to have the variants as rownames
    SNP_INFO[,c("protein")] <- apply(SNP_INFO,1,function(x){
      ifelse(x["protein"]!="",x["protein"],
             ifelse(x["cdna"]!="",paste(x["gene"],x["cdna"],sep="_"),
                    ifelse(x["amplicon"]%in%c("MSK_RL_AMP54", "MSK_RL_AMP55", "MSK_RL_AMP56","MSK_RL_AMP57","MSK_RL_AMP58"),paste("FLT3_INS",x["POS"],x["ALT"],sep="_"),"")))
    })# - variant metadata
  } else {
    SNP_INFO <- read.csv(paste("./",sample,"/","Variants.csv",sep=""))  
    rownames(SNP_INFO) <-colnames(AF)[-c(1:2)]  # - renames SNP info to have the variants as rownames
    colnames(SNP_INFO)[3] <- "protein"
    SNP_INFO$cDNA <- as.character(SNP_INFO$cDNA)
    SNP_INFO$protein <- as.character(SNP_INFO$protein)
    SNP_INFO$Variant <- as.character(SNP_INFO$Variant)
    SNP_INFO$Gene <- as.character(SNP_INFO$Gene)
    SNP_INFO[,c("protein")] <- apply(SNP_INFO,1,function(x){
      ifelse(x["protein"]!="",x["protein"],
             ifelse(x["cDNA"]!="",paste(x["Gene"],x["cDNA"],sep="_"),
                    ifelse(grepl("^chr13",x["Variant"]),paste("FLT3_INS",x["Variant"],sep="_"),"")))
    })# - variant metadata
  }
  
  SNP_changes_of_interest <- LAM_VARIANTS_SNP
  protein_changes_of_interest <- LAM_VARIANTS_Protein
  SNP_protein_key <- data.frame("Mutant"=SNP_changes_of_interest,"Protein"=protein_changes_of_interest)
  
  # Cell and clone selection
  distribution_of_unknowns_by_variant <- apply(NGT[,SNP_changes_of_interest],MARGIN=2,table)
  #print(distribution_of_unknowns_by_variant)
  cells_with_unknown<-NGT[apply(NGT[,SNP_changes_of_interest],MARGIN=1,FUN=function(x){sum(x==3)>0}),"Cell"]
  matrix_of_interest<-NGT[!NGT$Cell%in%cells_with_unknown,SNP_changes_of_interest]
  bulk_VAF_order<-SNP_INFO[colnames(matrix_of_interest)[order(colSums(matrix_of_interest),decreasing=TRUE)],"protein"]
  bulk_VAF_order <- bulk_VAF_order[!duplicated(bulk_VAF_order)]
  matrix_of_interest$Clone <- apply(matrix_of_interest,1,function(x){paste(x,sep="_",collapse="_")})
  clonal_abundance <- sort(table(matrix_of_interest$Clone))
  clones_to_include<-names(clonal_abundance)[clonal_abundance>length(rownames(matrix_of_interest))*0]
  matrix_subset_clones <- matrix_of_interest[matrix_of_interest$Clone%in%clones_to_include,]
  
  #Average depth per cell per allele
  DP_subset<- DP[-cells_with_unknown,c("Sample","Cell",rownames(SNP_INFO)[SNP_INFO$protein%in%protein_changes_of_interest])]
  NGT_subset <- NGT[NGT$Cell%in%DP_subset$Cell,colnames(DP_subset)]
  DP_NGT_melt<-inner_join(melt(DP_subset[,-1],id.var="Cell"),melt(NGT[,-1],id.var="Cell"),by=c("Cell","variable"))
  colnames(DP_NGT_melt) <- c("Cell","Mutant","Depth","Genotype")
  DP_NGT_melt$Genotype <- ifelse(DP_NGT_melt$Genotype==3,"Unknown",
                                 ifelse(DP_NGT_melt$Genotype==0,"WT",
                                        ifelse(DP_NGT_melt$Genotype==1,"Heterozygous",
                                               ifelse(DP_NGT_melt$Genotype==2,"Homozygous",DP_NGT_melt$Genotype))))
  DP_NGT_melt$Genotype  <- factor(DP_NGT_melt$Genotype ,levels=c("WT","Heterozygous","Homozygous","Unknown"))
  
  summarized<-data.frame(DP_NGT_melt%>%group_by(Mutant,Genotype)%>%summarise("Depth"=mean(Depth)))
  
  
  # Verify variants selected are unique and encode protein changes
  !any(duplicated(SNP_INFO[colnames(matrix_subset_clones),"protein"])) #Should return TRUE
  colnames(matrix_subset_clones)<- c(as.character(SNP_INFO[colnames(matrix_subset_clones)[1:length(SNP_changes_of_interest)],"protein"]),"Clone")
  
  # Aggregate data to get clone counts and architecture.
  dedup<-matrix_subset_clones[!duplicated(matrix_subset_clones),]
  clonal_architecture <- melt(dedup)
  colnames(clonal_architecture) <- c("Clone","Mutant","Genotype")
  clonal_architecture$Genotype <- ifelse(clonal_architecture$Genotype==3,NA,
                                         ifelse(clonal_architecture$Genotype==0,"WT",
                                                ifelse(clonal_architecture$Genotype==1,"Heterozygous",
                                                       ifelse(clonal_architecture$Genotype==2,"Homozygous",clonal_architecture$Genotype))))
  clonal_architecture$Genotype  <- factor(clonal_architecture$Genotype ,levels=c("WT","Heterozygous","Homozygous"))
  clonal_abundance <- data.frame("Count"=as.matrix(table(matrix_subset_clones$Clone)))
  clonal_abundance$Clone <- rownames(clonal_abundance)
  clonal_abundance <- data.table(clonal_abundance,key="Count")
  clonal_abundance$Clone <- factor(clonal_abundance$Clone,levels=rev(c(clonal_abundance$Clone)))
  clonal_architecture$Clone <- factor(as.character(clonal_architecture$Clone),levels=rev(as.character(clonal_abundance$Clone)))
  clonal_architecture$Mutant <-factor(clonal_architecture$Mutant , levels=rev(as.character(bulk_VAF_order)))
  
  # Generate clonal architecture heatmap
  gg_heatmap <- ggplot(data = clonal_architecture, aes(x = Clone, y = factor(Mutant), fill = Genotype)) + 
    geom_tile() + 
    scale_fill_manual(values=c("WT"=brewer.pal(7,"Reds")[1],
                               "Heterozygous"=brewer.pal(7,"Reds")[3],
                               "Homozygous"=brewer.pal(7,"Reds")[6],
                               "Unknown"="grey50"),"Genotype")  +
    theme_classic() +
    ylab("Mutation")+
    theme(legend.position = "bottom", legend.direction = "horizontal",
          axis.text.x = element_blank(), 
          axis.line=element_blank(),
          axis.title.x=element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin=unit(c(0,1,1,1),"cm"))
  
  # Generate clonal abundance barplot
  gg_clonal_barplot <- ggplot(data = data.frame(clonal_abundance), aes(x = factor(Clone), y = Count)) + 
    geom_bar(stat = "identity", aes(fill = Count)) + theme_gray() +
    theme_classic()+
    ylim(0,max(clonal_abundance$Count)*1.3) + 
    ylab("Cell Count")+
    geom_text(aes(label=Count), position=position_dodge(width=0.9), vjust=-0.25)+
    scale_fill_distiller(name = "Value", palette = "Reds", direction = 1) +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(),  axis.ticks.x = element_blank(), axis.line=element_blank(),
          legend.position = "none",
          plot.margin=unit(c(1,1,-0.2,1),"cm"))
  
  # Generate genotype calls by allele
  distribution_of_unknowns_by_variant <- t(apply(NGT[,SNP_changes_of_interest],MARGIN=2,function(x){table(factor(x,levels=c(0,1,2,3)))}))
  rownames(distribution_of_unknowns_by_variant) <- SNP_INFO[rownames(distribution_of_unknowns_by_variant),"protein"]
  colnames(distribution_of_unknowns_by_variant) <- c("WT","Heterozygous","Homozygous","Unknown")
  allele_melted <-melt(distribution_of_unknowns_by_variant) 
  colnames(allele_melted) <- c("Mutant","Genotype","Cells")
  allele_melted$Genotype <- factor(allele_melted$Genotype,levels=rev(c("WT","Heterozygous","Homozygous","Unknown")))
  allele_melted$Mutant <-factor(allele_melted$Mutant , levels=as.character(bulk_VAF_order))
  
  #Generate allele barplot
  gg_alleles <-  ggplot(data = data.frame(allele_melted), aes(x = factor(Mutant), y = Cells)) + 
    geom_bar(stat = "identity", aes(fill = Genotype),position="stack") + theme_gray() +
    theme_minimal()+
    ylab("Cell Count")+
    scale_fill_manual(values=c("WT"=brewer.pal(7,"Reds")[1],
                               "Heterozygous"=brewer.pal(7,"Reds")[3],
                               "Homozygous"=brewer.pal(7,"Reds")[6],
                               "Unknown"="grey50"),"Genotype") +
    theme(axis.text.x=element_text(angle=30,hjust=1,vjust=1))
  
  
  #Generate depth by alelle by mutant                       
  gg_depth <- ggplot(summarized,aes(x=Mutant,color=Genotype,y=Depth))+
    geom_point(position=position_dodge(width=0.4))+
    scale_color_manual(values=c("WT"=brewer.pal(7,"Reds")[1],
                                "Heterozygous"=brewer.pal(7,"Reds")[3],
                                "Homozygous"=brewer.pal(7,"Reds")[6],
                                "Unknown"="grey50"),"Genotype") +
    theme(axis.text.x=element_text(angle=30,hjust=1,vjust=1),
          plot.margin=unit(c(1,1,1,2),"cm"))
  
  # Compute spacing of graphs
  plots <- list(gg_depth,gg_clonal_barplot,gg_alleles,gg_heatmap)
  grobs <- list()
  widths <- list()
  for (i in 1:length(plots)){
    grobs[[i]] <- ggplotGrob(plots[[i]])
    widths[[i]] <- grobs[[i]]$widths[2:5]
  }
  maxwidth <- do.call(grid::unit.pmax, widths)
  for (i in 1:length(grobs)){
    grobs[[i]]$widths[2:5] <- as.list(maxwidth)
  }
  
  # Generate plot
  grob.title <- textGrob(sample, hjust = 0.5, vjust = 0.5, gp = gpar(fontsize = 20))
  pdf(paste(output_folder,"Graphs/",sample,".pdf",sep=""), width = 16, height = 8) # Open a new pdf file
  grid.arrange(grobs=grobs, ncol = 2, top = grob.title)
  graphics.off()
  
  
  if(diversity==TRUE){
    clonal_abundance_mat <- clonal_abundance
    clonal_abundance$Clone <- rownames(clonal_abundance)
    clonal_abundance_mat$Clone <- rownames(clonal_abundance_mat)
    clonal_abundance_count <- data.table(clonal_abundance,key="Count")
    shannon <- diversity(clonal_abundance_mat[,"Count"])
    write.csv2(shannon,paste(output_folder,"Diversity_score/",sample,".csv",sep=""))
  }
}


stats_output<-function(sample){
  #Load in the data
  
  if(any(grepl("VARIANT_SELECTION_LAM",list.files(paste("./",sample,"/",sep=""))))){
    LAM_VARIANTS_MAT <- read.csv(paste("./",sample,"/",sample,"_VARIANT_SELECTION_LAM.csv",sep=""))  
    LAM_VARIANTS_SNP <- as.character((LAM_VARIANTS_MAT %>% filter(Include=="Include")%>%select(Mutant))[,1])
    LAM_VARIANTS_Protein <- as.character((LAM_VARIANTS_MAT %>% filter(Include=="Include")%>%select(Protein))[,1])
  } else {
    print("No variant selection performed!") 
    return("next")
  }
  AF <- read.csv(paste("./",sample,"/","AF.csv",sep=""))                    # - variant allele frequency
  DP <- read.csv(paste("./",sample,"/","DP.csv",sep=""))                     # - read depth
  GQ <- read.csv(paste("./",sample,"/","GQ.csv",sep=""))                    # - quality scores
  NGT <- read.csv(paste("./",sample,"/","NGT.csv",sep="")) # - genotype call converted to categorical (0-reference, 1-heterozygous mutation, 2-homozygous mutation, 3-unknown)
  
  if(any(grepl("SNP_INFO.csv",list.files(paste("./",sample,"/",sep=""))))){
    SNP_INFO <- read.csv(paste("./",sample,"/","SNP_INFO.csv",sep="")) 
    rownames(SNP_INFO) <-colnames(AF)[-c(1:2)]  # - renames SNP info to have the variants as rownames
    SNP_INFO[,c("protein")] <- apply(SNP_INFO,1,function(x){
      ifelse(x["protein"]!="",x["protein"],
             ifelse(x["cdna"]!="",paste(x["gene"],x["cdna"],sep="_"),
                    ifelse(x["amplicon"]%in%c("MSK_RL_AMP54", "MSK_RL_AMP55", "MSK_RL_AMP56","MSK_RL_AMP57","MSK_RL_AMP58"),paste("FLT3_INS",x["POS"],x["ALT"],sep="_"),"")))
    })# - variant metadata
  } else {
    SNP_INFO <- read.csv(paste("./",sample,"/","Variants.csv",sep=""))  
    rownames(SNP_INFO) <-colnames(AF)[-c(1:2)]  # - renames SNP info to have the variants as rownames
    colnames(SNP_INFO)[3] <- "protein"
    SNP_INFO$cDNA <- as.character(SNP_INFO$cDNA)
    SNP_INFO$protein <- as.character(SNP_INFO$protein)
    SNP_INFO$Variant <- as.character(SNP_INFO$Variant)
    SNP_INFO$Gene <- as.character(SNP_INFO$Gene)
    SNP_INFO[,c("protein")] <- apply(SNP_INFO,1,function(x){
      ifelse(x["protein"]!="",x["protein"],
             ifelse(x["cDNA"]!="",paste(x["Gene"],x["cDNA"],sep="_"),
                    ifelse(grepl("^chr13",x["Variant"]),paste("FLT3_INS",x["Variant"],sep="_"),"")))
    })# - variant metadata
  }
  
  SNP_changes_of_interest <- LAM_VARIANTS_SNP
  protein_changes_of_interest <- LAM_VARIANTS_Protein
  SNP_protein_key <- data.frame("Mutant"=SNP_changes_of_interest,"Protein"=protein_changes_of_interest)
  
  # Cell and clone selection
  distribution_of_unknowns_by_variant <- apply(NGT[,SNP_changes_of_interest],MARGIN=2,table)
  #print(distribution_of_unknowns_by_variant)
  cells_with_unknown<-NGT[apply(NGT[,SNP_changes_of_interest],MARGIN=1,FUN=function(x){sum(x==3)>0}),"Cell"]
  matrix_of_interest<-NGT[!NGT$Cell%in%cells_with_unknown,SNP_changes_of_interest]
  bulk_VAF_order<-SNP_INFO[colnames(matrix_of_interest)[order(colSums(matrix_of_interest),decreasing=TRUE)],"protein"]
  bulk_VAF_order <- bulk_VAF_order[!duplicated(bulk_VAF_order)]
  matrix_of_interest$Clone <- apply(matrix_of_interest,1,function(x){paste(x,sep="_",collapse="_")})

  # Aggregate data to get clone counts and architecture.
  dedup<-matrix_of_interest[!duplicated(matrix_of_interest),]
  clonal_abundance <- data.frame("Count"=as.matrix(table(matrix_of_interest$Clone)))
  clonal_abundance$Clone <-  factor(rownames(clonal_abundance),levels=rev(c(rownames(clonal_abundance))))
  clonal_abundance$Count <- as.numeric(clonal_abundance$Count)
  clonal_abundance$Dominance <- apply(clonal_abundance,1,function(x){
    (as.numeric(x["Count"])/sum(clonal_abundance["Count"]) ) / (max(as.numeric(clonal_abundance[,"Count"]))/sum(clonal_abundance["Count"]))
  })
  clonal_abundance$Mutation_count <- apply(clonal_abundance,1,function(x){
    sum(as.numeric(strsplit( as.character(x[2]),split="_")[[1]]))
  })
  shannon <- diversity(clonal_abundance[,"Count"])
  print(shannon)
  matrix_output <- dedup%>%select(-Clone)
  matrix_output_t <- t(matrix_output)
  rownames(matrix_output_t) <- 0:(dim(matrix_output_t)[[1]]-1)
  write.table(shannon,paste0(output_folder,"Diversity_score/",sample,".txt"),sep="\t",row.names=FALSE)
  write.table(matrix_output_t,paste0(output_folder,"Genotype_matrix/",sample,".txt"),sep=" ",row.names=TRUE,quote=FALSE,col.names = FALSE)
  write.table(clonal_abundance,paste0(output_folder,"Summary_mat/",sample,".txt"),sep="\t",row.names=FALSE)
}


permute_data <- function(sample,nrun){
  #Load in the data
  
  if(any(grepl("VARIANT_SELECTION_LAM",list.files(paste("./",sample,"/",sep=""))))){
    LAM_VARIANTS_MAT <- read.csv(paste("./",sample,"/",sample,"_VARIANT_SELECTION_LAM.csv",sep=""))  
    LAM_VARIANTS_SNP <- as.character((LAM_VARIANTS_MAT %>% filter(Include=="Include")%>%select(Mutant))[,1])
    LAM_VARIANTS_Protein <- as.character((LAM_VARIANTS_MAT %>% filter(Include=="Include")%>%select(Protein))[,1])
  } else {
    print("No variant selection performed!") 
    return("next")
  }
  AF <- read.csv(paste("./",sample,"/","AF.csv",sep=""))                    # - variant allele frequency
  DP <- read.csv(paste("./",sample,"/","DP.csv",sep=""))                     # - read depth
  GQ <- read.csv(paste("./",sample,"/","GQ.csv",sep=""))                    # - quality scores
  NGT <- read.csv(paste("./",sample,"/","NGT.csv",sep="")) # - genotype call converted to categorical (0-reference, 1-heterozygous mutation, 2-homozygous mutation, 3-unknown)
  
  if(any(grepl("SNP_INFO.csv",list.files(paste("./",sample,"/",sep=""))))){
    SNP_INFO <- read.csv(paste("./",sample,"/","SNP_INFO.csv",sep="")) 
    rownames(SNP_INFO) <-colnames(AF)[-c(1:2)]  # - renames SNP info to have the variants as rownames
    SNP_INFO[,c("protein")] <- apply(SNP_INFO,1,function(x){
      ifelse(x["protein"]!="",x["protein"],
             ifelse(x["cdna"]!="",paste(x["gene"],x["cdna"],sep="_"),
                    ifelse(x["amplicon"]%in%c("MSK_RL_AMP54", "MSK_RL_AMP55", "MSK_RL_AMP56","MSK_RL_AMP57","MSK_RL_AMP58"),paste("FLT3_INS",x["POS"],x["ALT"],sep="_"),"")))
    })# - variant metadata
  } else {
    SNP_INFO <- read.csv(paste("./",sample,"/","Variants.csv",sep=""))  
    rownames(SNP_INFO) <-colnames(AF)[-c(1:2)]  # - renames SNP info to have the variants as rownames
    colnames(SNP_INFO)[3] <- "protein"
    SNP_INFO$cDNA <- as.character(SNP_INFO$cDNA)
    SNP_INFO$protein <- as.character(SNP_INFO$protein)
    SNP_INFO$Variant <- as.character(SNP_INFO$Variant)
    SNP_INFO$Gene <- as.character(SNP_INFO$Gene)
    SNP_INFO[,c("protein")] <- apply(SNP_INFO,1,function(x){
      ifelse(x["protein"]!="",x["protein"],
             ifelse(x["cDNA"]!="",paste(x["Gene"],x["cDNA"],sep="_"),
                    ifelse(grepl("^chr13",x["Variant"]),paste("FLT3_INS",x["Variant"],sep="_"),"")))
    })# - variant metadata
  }
  
  SNP_changes_of_interest <- LAM_VARIANTS_SNP
  protein_changes_of_interest <- LAM_VARIANTS_Protein
  SNP_protein_key <- data.frame("Mutant"=SNP_changes_of_interest,"Protein"=protein_changes_of_interest)
  
  # Cell and clone selection
  distribution_of_unknowns_by_variant <- apply(NGT[,SNP_changes_of_interest],MARGIN=2,table)
  #print(distribution_of_unknowns_by_variant)
  cells_with_unknown<-NGT[apply(NGT[,SNP_changes_of_interest],MARGIN=1,FUN=function(x){sum(x==3)>0}),"Cell"]
  matrix_of_interest<-NGT[!NGT$Cell%in%cells_with_unknown,SNP_changes_of_interest]
  bulk_VAF_order<-SNP_INFO[colnames(matrix_of_interest)[order(colSums(matrix_of_interest),decreasing=TRUE)],"protein"]
  bulk_VAF_order <- bulk_VAF_order[!duplicated(bulk_VAF_order)]
  matrix_of_interest$Clone <- apply(matrix_of_interest,1,function(x){paste(x,sep="_",collapse="_")})
  
  # Aggregate data to get clone counts and architecture.
  dedup<-matrix_of_interest[!duplicated(matrix_of_interest),]
  clonal_abundance <- data.frame("Count"=as.matrix(table(matrix_of_interest$Clone)))
  clonal_abundance$Clone <-  factor(rownames(clonal_abundance),levels=rev(c(rownames(clonal_abundance))))
  clonal_abundance$Count <- as.numeric(clonal_abundance$Count)
  
  clonal_architecture <- melt(dedup)
  colnames(clonal_architecture) <- c("Clone","Mutant","Genotype")
  clonal_architecture$Genotype <- ifelse(clonal_architecture$Genotype==3,NA,
                                         ifelse(clonal_architecture$Genotype==0,"WT",
                                                ifelse(clonal_architecture$Genotype==1,"Heterozygous",
                                                       ifelse(clonal_architecture$Genotype==2,"Homozygous",clonal_architecture$Genotype))))
  clonal_architecture$Genotype  <- factor(clonal_architecture$Genotype ,levels=c("WT","Heterozygous","Homozygous"))

  results_holder <- unlist(list(1:nrun))
  randomization_list<-lapply(results_holder,function(x){
    apply(matrix_of_interest[,1:length(colnames(matrix_of_interest))-1],2,function(x){
      y <- sample(x,length(x),replace=FALSE)
    })})
  
  #randomization_list<-plyr::alply(randomization_array,3)
  # names(randomization_list) <- NULL
  
  empirical_distribution<-suppressMessages(lapply(randomization_list, function(randomization){
    randomization$Clone <- apply(randomization,1,function(x){paste(x,sep="_",collapse="_")})
    loop_dedup<-matrix_of_interest[!duplicated(randomization),]
    loop_clonal_architecture <- melt(loop_dedup)
    colnames(loop_clonal_architecture) <- c("Clone","Mutant","Genotype")
    loop_clonal_architecture$Genotype <- ifelse(loop_clonal_architecture$Genotype==3,NA,
                                                ifelse(loop_clonal_architecture$Genotype==0,"WT",
                                                       ifelse(loop_clonal_architecture$Genotype==1,"Heterozygous",
                                                              ifelse(loop_clonal_architecture$Genotype==2,"Homozygous",clonal_architecture$Genotype))))
    loop_clonal_architecture$Genotype  <- factor(loop_clonal_architecture$Genotype ,levels=c("WT","Heterozygous","Homozygous"))
    loop_clonal_abundance <- data.frame("Count"=as.matrix(table(randomization$Clone)))
    loop_clonal_abundance$Clone <- rownames(loop_clonal_abundance)
    return(data.frame(loop_clonal_abundance))
  }))
  
  
  empirical_distribution_mat<-empirical_distribution%>%purrr::reduce(dplyr::full_join,by="Clone")
  colnames(empirical_distribution_mat)<-1:length(colnames(empirical_distribution_mat))
  rownames(empirical_distribution_mat) <- empirical_distribution_mat[,2]
  empirical_distribution_mat <- empirical_distribution_mat[,-2]
  empirical_distribution_mat[is.na(empirical_distribution_mat)]<-0
  
  clonal_abundance$Clone <- as.character(clonal_abundance$Clone)
  clonal_abundance$Count <- as.numeric(clonal_abundance$Count)
  clonal_abundance$pvalue <-apply(clonal_abundance,1,function(x){
  1-  sum(as.numeric(x[["Count"]])>empirical_distribution_mat[x[["Clone"]],])/nrun
    })
  
  
  #print(dedup)
  colnames(clonal_abundance)[2] <- paste(colnames(dedup)[1:length(SNP_changes_of_interest)],sep="_",collapse="__")
  clonal_abundance <- clonal_abundance[order(clonal_abundance[,"Count"],decreasing = TRUE),]
  clonal_abundance$Dominance <- apply(clonal_abundance,1,function(x){
    (as.numeric(x["Count"])/sum(clonal_abundance["Count"]) ) / (max(as.numeric(clonal_abundance[,"Count"]))/sum(clonal_abundance["Count"]))
      })
  clonal_abundance$Poisson <- apply(clonal_abundance,1,function(x){
        poisson.test(as.numeric(x["Count"]),mean(as.numeric(clonal_abundance$Count)),alternative="greater")$p.value
      })
  clonal_abundance$Mutation_count <- apply(clonal_abundance,1,function(x){
    sum(as.numeric(strsplit( as.character(x[2]),split="_")[[1]]))
  })
  
  write.table(clonal_abundance,paste0(output_folder,"Summary_mat/",sample,".txt"),sep="\t",row.names=FALSE)
  
}



replot_clonal_selection <- function(sample){
  
  if(any(grepl("VARIANT_SELECTION_LAM",list.files(paste("./",sample,"/",sep=""))))){
    LAM_VARIANTS_MAT <- read.csv(paste("./",sample,"/",sample,"_VARIANT_SELECTION_LAM.csv",sep=""))  
    LAM_VARIANTS_SNP <- as.character((LAM_VARIANTS_MAT %>% filter(Include=="Include")%>%select(Mutant))[,1])
    LAM_VARIANTS_Protein <- as.character((LAM_VARIANTS_MAT %>% filter(Include=="Include")%>%select(Protein))[,1])
  } else {
    print("No variant selection performed!") 
    return("next")
  }
  DP <- read.csv(paste("./",sample,"/","DP.csv",sep=""))                     # - read depth
  NGT <- read.csv(paste("./",sample,"/","NGT.csv",sep="")) # - genotype call converted to categorical (0-reference, 1-heterozygous mutation, 2-homozygous mutation, 3-unknown)
  
  if(any(grepl("SNP_INFO.csv",list.files(paste("./",sample,"/",sep=""))))){
    SNP_INFO <- read.csv(paste("./",sample,"/","SNP_INFO.csv",sep="")) 
    rownames(SNP_INFO) <-colnames(NGT)[-c(1:2)]  # - renames SNP info to have the variants as rownames
    SNP_INFO[,c("protein")] <- apply(SNP_INFO,1,function(x){
      ifelse(x["protein"]!="",x["protein"],
             ifelse(x["cdna"]!="",paste(x["gene"],x["cdna"],sep="_"),
                    ifelse(x["amplicon"]%in%c("MSK_RL_AMP54", "MSK_RL_AMP55", "MSK_RL_AMP56","MSK_RL_AMP57","MSK_RL_AMP58"),paste("FLT3_INS",x["POS"],x["ALT"],sep="_"),"")))
    })# - variant metadata
  } else {
    SNP_INFO <- read.csv(paste("./",sample,"/","Variants.csv",sep=""))  
    rownames(SNP_INFO) <-colnames(NGT)[-c(1:2)]  # - renames SNP info to have the variants as rownames
    colnames(SNP_INFO)[3] <- "protein"
    SNP_INFO$cDNA <- as.character(SNP_INFO$cDNA)
    SNP_INFO$protein <- as.character(SNP_INFO$protein)
    SNP_INFO$Variant <- as.character(SNP_INFO$Variant)
    SNP_INFO$Gene <- as.character(SNP_INFO$Gene)
    SNP_INFO[,c("protein")] <- apply(SNP_INFO,1,function(x){
      ifelse(x["protein"]!="",x["protein"],
             ifelse(x["cDNA"]!="",paste(x["Gene"],x["cDNA"],sep="_"),
                    ifelse(grepl("^chr13",x["Variant"]),paste("FLT3_INS",x["Variant"],sep="_"),"")))
    })# - variant metadata
  }
  
  SNP_changes_of_interest <- LAM_VARIANTS_SNP
  protein_changes_of_interest <- LAM_VARIANTS_Protein
  SNP_protein_key <- data.frame("Mutant"=SNP_changes_of_interest,"Protein"=protein_changes_of_interest)
  rownames(SNP_protein_key) <- SNP_protein_key$Mutant
  # Cell and clone selection
  distribution_of_unknowns_by_variant <- apply(NGT[,SNP_changes_of_interest],MARGIN=2,table)
  #print(distribution_of_unknowns_by_variant)
  cells_with_unknown<-NGT[apply(NGT[,SNP_changes_of_interest],MARGIN=1,FUN=function(x){sum(x==3)>0}),"Cell"]
  matrix_of_interest<-NGT[!NGT$Cell%in%cells_with_unknown,SNP_changes_of_interest]
  bulk_VAF_order<-SNP_INFO[colnames(matrix_of_interest)[order(colSums(matrix_of_interest),decreasing=TRUE)],"protein"]
  bulk_VAF_order <- bulk_VAF_order[!duplicated(bulk_VAF_order)]
  matrix_of_interest$Clone <- apply(matrix_of_interest,1,function(x){paste(x,sep="_",collapse="_")})
  
  dedup<-matrix_of_interest[!duplicated(matrix_of_interest),]
  colnames(dedup) <- SNP_protein_key[colnames(dedup),"Protein"]
  clonal_architecture <- melt(dedup)
  colnames(clonal_architecture) <- c("Clone","Mutant","Genotype")
  clonal_architecture$Genotype <- ifelse(clonal_architecture$Genotype==3,NA,
                                         ifelse(clonal_architecture$Genotype==0,"WT",
                                                ifelse(clonal_architecture$Genotype==1,"Heterozygous",
                                                       ifelse(clonal_architecture$Genotype==2,"Homozygous",clonal_architecture$Genotype))))
  clonal_architecture$Genotype  <- factor(clonal_architecture$Genotype ,levels=c("WT","Heterozygous","Homozygous"))
  clonal_abundance <- data.frame("Count"=as.matrix(table(matrix_of_interest$Clone)))
  clonal_abundance$Clone <- rownames(clonal_abundance)
  clonal_abundance <- data.table(clonal_abundance,key="Count")
  clonal_abundance$Clone <- factor(clonal_abundance$Clone,levels=rev(c(clonal_abundance$Clone)))
  clonal_architecture$Clone <- factor(as.character(clonal_architecture$Clone),levels=rev(as.character(clonal_abundance$Clone)))
  clonal_architecture$Mutant <-factor(clonal_architecture$Mutant , levels=rev(as.character(bulk_VAF_order)))
  
  clones_import <- read.delim(paste0("/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/08_08_19/Summary_mat/",sample,".txt"),sep="\t")
  colnames(clones_import)[2] <- "Clone"
  select_clones_full <- clones_import %>% filter(pvalue<0.1|Poisson<0.05)%>%select(Clone)
  select_clones_dominant <- clones_import %>% filter(pvalue<0.1|Poisson<0.05)%>%filter(Dominance>0.05)%>%select(Clone)

  
  
  if(length(select_clones_dominant[,1])>3){
    clonal_architecture_subset <- clonal_architecture%>%filter(Clone%in%select_clones_dominant[,1])
    clonal_abundance_subset <- clonal_abundance%>%filter(Clone%in%select_clones_dominant[,1])
    name_var <- "VAF Subset"
  } else if(length(select_clones_full[,1])>20){
    print("Complex clonal picture, skipping")
    return(next)
  } else if(length(select_clones_full[,1])<2){
    print("Simple clonal picture, plotting all clones")
    clonal_architecture_subset <- clonal_architecture
    clonal_abundance_subset <- clonal_abundance
    name_var <- "All Clones"
  } else {
    clonal_architecture_subset <- clonal_architecture%>%filter(Clone%in%select_clones_full[,1])
    clonal_abundance_subset <- clonal_abundance%>%filter(Clone%in%select_clones_full[,1])
    name_var <- "Full"
  }
  # Generate clonal architecture heatmap
  gg_heatmap <- ggplot(data = clonal_architecture_subset, aes(x = Clone, y = factor(Mutant), fill = Genotype)) + 
    geom_tile() + 
    scale_fill_manual(values=c("WT"=brewer.pal(7,"Reds")[1],
                               "Heterozygous"=brewer.pal(7,"Reds")[3],
                               "Homozygous"=brewer.pal(7,"Reds")[6],
                               "Unknown"="grey50"),"Genotype")  +
    theme_classic() +
    ylab("Mutation")+
    theme(legend.position = "bottom", legend.direction = "horizontal",
          axis.text.x = element_blank(), 
          axis.line=element_blank(),
          axis.title.x=element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin=unit(c(0,1,1,1),"cm"))
  
  # Generate clonal abundance barplot
  gg_clonal_barplot <- ggplot(data = data.frame(clonal_abundance_subset), aes(x = factor(Clone), y = Count)) + 
    geom_bar(stat = "identity", aes(fill = Count)) + theme_gray() +
    theme_classic()+
    ylim(0,max(clonal_abundance$Count)*1.3) + 
    ylab("Cell Count")+
    geom_text(aes(label=Count), position=position_dodge(width=0.9), vjust=-0.25)+
    scale_fill_distiller(name = "Value", palette = "Reds", direction = 1) +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(),  axis.ticks.x = element_blank(), axis.line=element_blank(),
          legend.position = "none",
          plot.margin=unit(c(1,1,-0.2,1),"cm"))

  # Compute spacing of graphs
  plots <- list(gg_clonal_barplot,gg_heatmap)
  grobs <- list()
  widths <- list()
  for (i in 1:length(plots)){
    grobs[[i]] <- ggplotGrob(plots[[i]])
    widths[[i]] <- grobs[[i]]$widths[2:5]
  }
  maxwidth <- do.call(grid::unit.pmax, widths)
  for (i in 1:length(grobs)){
    grobs[[i]]$widths[2:5] <- as.list(maxwidth)
  }
  
  # Generate plot
  grob.title <- textGrob(paste(sample,"-",name_var), hjust = 0.5, vjust = 0.5, gp = gpar(fontsize = 20))
  pdf(paste(output_folder,"Graphs_significant/",sample,".pdf",sep=""), width = 16, height = 8) # Open a new pdf file
  grid.arrange(grobs=grobs, ncol = 1, top = grob.title)
  graphics.off()

}


plot_clonal_heatmap_and_barplot <-function(sample,output_folder,clonal_abundance,clonal_architecture,clone_cell_count,shrink)
{
  #clonal_abundance_subset <-data.frame( clonal_abundance[[sample]])%>%filter(Count>=5)
  #clonal_architecture_subset <- data.frame(clonal_architecture[[sample]])%>%
   #                                         filter(Clone%in%as.character(clonal_abundance_subset$Clone))
  clonal_abundance_subset<-final_sample_summary[[sample]]$Clones#clonal_abundance
  clonal_abundance_subset$Clone<-factor(final_sample_summary[[sample]]$Clones[,"Clone"],levels=rev(c(final_sample_summary[[sample]]$Clones[,"Clone"])))#clonal_abundance
  clonal_architecture_subset<-final_sample_summary[[sample]]$Architecture%>%
                                      filter(Clone%in%as.character(clonal_abundance_subset$Clone))
  
  clonal_architecture_subset$Clone <- factor(clonal_architecture_subset$Clone, levels=levels(clonal_abundance_subset$Clone))
  gg_heatmap <- ggplot(data = clonal_architecture_subset, aes(x = Clone, y = factor(Mutant,levels=rev(levels(factor(Mutant)))), fill = Genotype)) + 
    geom_tile() +# scale_y_discrete(limits = rev(levels(Mutant)))+

    scale_fill_manual(values=c("WT"=brewer.pal(7,"Reds")[1],
                               "Heterozygous"=brewer.pal(7,"Reds")[3],
                               "Homozygous"=brewer.pal(7,"Reds")[6],
                               "Unknown"="grey50"),"Genotype")  +
    theme_classic(base_size=6) +
    ylab("Mutation")+
    theme(legend.position = "bottom", legend.direction = "horizontal",
          axis.text.x = element_blank(), 
          axis.line=element_blank(),
          axis.title.x=element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin=unit(c(0,1,1,1),"cm"))
  
  # Generate clonal abundance barplot
  gg_clonal_barplot <- ggplot(data = data.frame(clonal_abundance_subset), aes(x = Clone, y = Count)) + 
    geom_bar(stat = "identity", aes(fill = Count)) + theme_gray() +
    theme_classic(base_size=6)+
    ylim(0,max(clonal_abundance_subset$Count)*1.3) + 
    ylab("Cell Count")+
    geom_text(aes(label=Count), position=position_dodge(width=0.9), vjust=-0.25,size=1)+
    scale_fill_distiller(name = "Value", palette = "Reds", direction = 1) +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(),  axis.ticks.x = element_blank(), axis.line=element_blank(),
          legend.position = "none",
          plot.margin=unit(c(1,1,-0.2,1),"cm"))
  
  grob.title <- textGrob(paste(sample,"-"), hjust= 0.5, vjust = 0.5, gp = gpar(fontsize = 12))
  final_plot<-plot_grid(grob.title,gg_clonal_barplot,gg_heatmap,ncol=1,align="hv",axis="l",rel_heights = c(0.1,1,0.75))
 
   save_plot(paste(output_folder,sample,".pdf",sep=""),final_plot, ncol=1) # Open a new pdf file
}

generate_and_plot_cooccurence <- function(mut_mat){
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
    geom_tile(aes(fill = factor(value)), color='grey90') +
    scale_fill_manual(values=c("-1"="firebrick3","0"="white","1"="steelblue2"),"Correlation",
                      labels=c("Mutually Exclusive","Not Significant","Mutually Inclusive"))+
    theme_classic(base_size=10)+xlab("")+ylab("")+
    theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1),
          axis.line = element_blank(),
          legend.position = c(0.8,1), 
          legend.justification = c(1, 1),
          legend.direction = "vertical")+
    theme(legend.key.size = unit(0.5,"line"))
  return(list("plot"=grob_corrplot,
         "data"=cooccur_mat))
}

diversity_from_bulk_VAF <- function(sample){
    if(any(grepl("VARIANT_SELECTION_LAM",list.files(paste("./",sample,"/",sep=""))))){
      LAM_VARIANTS_MAT <- read.csv(paste("./",sample,"/",sample,"_VARIANT_SELECTION_LAM.csv",sep=""))  
      LAM_VARIANTS_SNP <- as.character((LAM_VARIANTS_MAT %>% filter(Include=="Include")%>%select(Mutant))[,1])
      LAM_VARIANTS_Protein <- as.character((LAM_VARIANTS_MAT %>% filter(Include=="Include")%>%select(Protein))[,1])
    } else {
      print("No variant selection performed!") 
      return("next")
    }
    #AF <- read.csv(paste("./",sample,"/","AF.csv",sep=""))                    # - variant allele frequency
    #DP <- read.csv(paste("./",sample,"/","DP.csv",sep=""))                     # - read depth
    #GQ <- read.csv(paste("./",sample,"/","GQ.csv",sep=""))                    # - quality scores
    NGT <- read.csv(paste("./",sample,"/","NGT.csv",sep="")) # - genotype call converted to categorical (0-reference, 1-heterozygous mutation, 2-homozygous mutation, 3-unknown)
    
    if(any(grepl("SNP_INFO.csv",list.files(paste("./",sample,"/",sep=""))))){
      SNP_INFO <- read.csv(paste("./",sample,"/","SNP_INFO.csv",sep="")) 
      rownames(SNP_INFO) <-colnames(NGT)[-c(1:2)]  # - renames SNP info to have the variants as rownames
      SNP_INFO[,c("protein")] <- apply(SNP_INFO,1,function(x){
        ifelse(x["protein"]!="",x["protein"],
               ifelse(x["cdna"]!="",paste(x["gene"],x["cdna"],sep="_"),
                      ifelse(x["amplicon"]%in%c("MSK_RL_AMP54", "MSK_RL_AMP55", "MSK_RL_AMP56","MSK_RL_AMP57","MSK_RL_AMP58"),paste("FLT3_INS",x["POS"],x["ALT"],sep="_"),"")))
      })# - variant metadata
    } else {
      SNP_INFO <- read.csv(paste("./",sample,"/","Variants.csv",sep=""))  
      rownames(SNP_INFO) <-colnames(NGT)[-c(1:2)]  # - renames SNP info to have the variants as rownames
      colnames(SNP_INFO)[3] <- "protein"
      SNP_INFO$cDNA <- as.character(SNP_INFO$cDNA)
      SNP_INFO$protein <- as.character(SNP_INFO$protein)
      SNP_INFO$Variant <- as.character(SNP_INFO$Variant)
      SNP_INFO$Gene <- as.character(SNP_INFO$Gene)
      SNP_INFO[,c("protein")] <- apply(SNP_INFO,1,function(x){
        ifelse(x["protein"]!="",x["protein"],
               ifelse(x["cDNA"]!="",paste(x["Gene"],x["cDNA"],sep="_"),
                      ifelse(grepl("^chr13",x["Variant"]),paste("FLT3_INS",x["Variant"],sep="_"),"")))
      })# - variant metadata
    }
    
    SNP_changes_of_interest <- LAM_VARIANTS_SNP
    protein_changes_of_interest <- LAM_VARIANTS_Protein
    SNP_protein_key <- data.frame("Mutant"=SNP_changes_of_interest,"Protein"=protein_changes_of_interest)
    
    # Cell and clone selection
    distribution_of_unknowns_by_variant <- apply(NGT[,SNP_changes_of_interest],MARGIN=2,table)
    #print(distribution_of_unknowns_by_variant)
    cells_with_unknown<-NGT[apply(NGT[,SNP_changes_of_interest],MARGIN=1,FUN=function(x){sum(x==3)>0}),"Cell"]
    matrix_of_interest<-NGT[!NGT$Cell%in%cells_with_unknown,SNP_changes_of_interest]
    bulk_VAF_order<-SNP_INFO[colnames(matrix_of_interest)[order(colSums(matrix_of_interest),decreasing=TRUE)],"protein"]
    bulk_VAF_order <- bulk_VAF_order[!duplicated(bulk_VAF_order)]
    matrix_of_interest$Clone <- apply(matrix_of_interest,1,function(x){paste(x,sep="_",collapse="_")})
    
    output<-diversity(colSums(matrix_of_interest[colnames(matrix_of_interest)!="Clone"])/length(rownames(matrix_of_interest)))
    write.table(output,paste0(output_folder,"Diversity_score_faux_VAF/",sample,".txt"),sep="\t")
  
}

## Old Functions
plot_variants_AF_and_DP <- function(sample){
  
  #Load in the data
  AF <- read.csv(paste("./",sample,"/","AF.csv",sep=""))                    # - variant allele frequency
  DP <- read.csv(paste("./",sample,"/","DP.csv",sep=""))                     # - read depth
  GQ <- read.csv(paste("./",sample,"/","GQ.csv",sep=""))                    # - quality scores
  NGT <- read.csv(paste("./",sample,"/","NGT.csv",sep=""))                 # - genotype call converted to categorical (0-reference, 1-heterozygous mutation, 2-homozygous mutation, 3-unknown)
  SNP_INFO <- read.csv(paste("./",sample,"/","SNP_INFO.csv",sep=""))         # - variant metadata
  rownames(SNP_INFO) <-colnames(AF)[-c(1:2)]  # - renames SNP info to have the variants as rownames
  
  AF_NGT_Melt<-inner_join(AF %>% unite("Sample_cell",c("Sample","Cell")) %>%melt("Sample_cell"),
                          NGT %>% unite("Sample_cell",c("Sample","Cell")) %>%melt("Sample_cell"),
                          by=c("Sample_cell","variable"))
  colnames(AF_NGT_Melt) <- c("Cell","Variant","AF","NGT")
  gg_AF_NGT<-  ggplot(AF_NGT_Melt,aes(x=as.factor(NGT),y=AF,color=NGT))+ geom_quasirandom()+facet_wrap(~Variant,ncol=3)
  
  DP_NGT_Melt<-inner_join(DP %>% unite("Sample_cell",c("Sample","Cell")) %>%melt("Sample_cell"),
                          NGT %>% unite("Sample_cell",c("Sample","Cell")) %>%melt("Sample_cell"),
                          by=c("Sample_cell","variable"))
  colnames(DP_NGT_Melt) <- c("Cell","Variant","DP","NGT")
  gg_DP_NGT<- ggplot(DP_NGT_Melt,aes(x=as.factor(NGT),y=DP,fill=NGT))+geom_boxplot()+facet_wrap(~Variant,ncol=3,scale="free_y")
  
  
  ggsave(paste0("./",sample,"/",sample,"_DP_NGT.pdf"),gg_DP_NGT,width=10,height=40)
  ggsave(paste0("./",sample,"/",sample,"_AF_NGT.pdf"),gg_AF_NGT,width=10,height=40)
  
}


plot_COMPASS_data <-function(sample,output_folder,melted_protein_mat,clonal_abundance,clonal_architecture,clone_cell_count,shrink)
  {
  melted_protein_mat <-melted_protein_mat[[sample]]
    clonal_abundance_subset <-data.frame( clonal_abundance[[sample]])#%>%filter(Count>=5)
    clonal_architecture_subset <- data.frame(clonal_architecture[[sample]])#%>%
    #  filter(Clone%in%as.character(clonal_abundance_subset$Clone))
    genes_to_display <-unique(as.character(clonal_architecture_subset[clonal_architecture_subset$Genotype%in%c("Heterozygous","Homozygous"),"Mutant"]))
    
    if(shrink==TRUE){
      clonal_architecture_subset <-clonal_architecture_subset%>%filter(Mutant%in%genes_to_display)
    }
    clonal_architecture_subset$Clone <- factor(clonal_architecture_subset$Clone, levels=levels(clonal_abundance_subset$Clone))
    melted_protein_mat_subset<-melted_protein_mat%>%filter(Clone%in%as.character(clonal_architecture_subset$Clone))
    melted_protein_mat_subset$Clone <-factor(melted_protein_mat_subset$Clone, levels=levels(clonal_abundance_subset$Clone))
    gg_protein_heatmap<-ggplot(melted_protein_mat_subset, aes(y=variable, x=(Clone))) + geom_tile(aes(fill = value),colour = "white") + 
      scale_fill_distiller(palette = "PRGn")+
      theme_minimal(base_size=6) +
      theme( axis.text.x = element_blank(), axis.title.x = element_blank(), 
             legend.position = "right",legend.direction = "vertical",
             axis.ticks.x = element_blank(),
             plot.margin=unit(c(0,0,0,0),"cm"))
    
     gg_heatmap <- ggplot(data = clonal_architecture_subset, aes(x = Clone, y = factor(Mutant,levels=rev(levels(factor(Mutant)))), fill = Genotype)) + 
      geom_tile() +# scale_y_discrete(limits = rev(levels(Mutant)))+
      scale_fill_manual(values=c("WT"=brewer.pal(7,"Reds")[1],
                                 "Heterozygous"=brewer.pal(7,"Reds")[3],
                                 "Homozygous"=brewer.pal(7,"Reds")[6],
                                 "Unknown"="grey50"),"Genotype")  +
      theme_classic(base_size=6) +
      ylab("Mutation")+
      theme(legend.position = "right", legend.direction = "vertical",
            axis.text.x = element_blank(), 
            axis.line=element_blank(),
            axis.title.x=element_blank(),
            axis.ticks.x = element_blank(),
            plot.margin=unit(c(0,0,0,0),"cm"))
    
    clone_colors <- alphabet()[1:length(unique(clonal_abundance_subset$Clone))]
    names(clone_colors) <-levels(clonal_abundance_subset$Clone)
    # Generate clonal abundance barplot
    gg_clonal_barplot <- ggplot(data = data.frame(clonal_abundance_subset), aes(x = Clone, y = Count)) + 
      geom_bar(stat = "identity", aes(fill = Clone)) + theme_gray() +
      theme_classic(base_size=6)+
      ylim(0,max(clonal_abundance_subset$Count)*1.3) + 
      ylab("Cell Count")+
      geom_text(aes(label=Count), position=position_dodge(width=0.9), vjust=-0.25,size=1)+
      scale_fill_manual(values=clone_colors)+
      theme(axis.title.x = element_blank(), axis.text.x = element_blank(),  axis.ticks.x = element_blank(), axis.line=element_blank(),
            legend.position = "none",
            plot.margin=unit(c(0,0,0,0),"cm"))
    
    grob.title <- textGrob(paste(sample,"-"), hjust= 0.5, vjust = 0.5, gp = gpar(fontsize = 12))
    final_plot<-plot_grid(gg_clonal_barplot,addSmallLegend(gg_heatmap),addSmallLegend(gg_protein_heatmap),ncol=1,align="hv",axis="lrtb",rel_heights = c(1,(length(genes_to_display)*0.15),1))
    
    save_plot(paste(output_folder,sample,".pdf",sep=""),final_plot, ncol=1) # Open a new pdf file
  }

addSmallLegend <- function(myPlot, pointSize = 3, textSize = 8, spaceLegend = 0.5) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_text(size = textSize), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
}


create_reward_matrix_deprecated<-function(Known_mat,weights){
  set.seed(68864)
  names(weights) <- apply(Known_mat,2,function(x){paste(x,sep="_",collapse="_")})
  num_type <- 2
  num_mutations <- nrow(Known_mat); 
  num_clones <- ncol(Known_mat)
  num_states <- num_type^num_mutations
  
  states<-data.frame(expand.grid(rep(list(0:num_type), num_mutations)))
  state_interactions<-expand.grid(apply(states[,1:num_mutations],1,function(x){paste(x,collapse="_",sep="_")}),
                                  apply(states[,1:num_mutations],1,function(x){paste(x,collapse="_",sep="_")}))
  
  states$Reward_states<-ifelse(apply(states,1,function(x){
    any(apply(Known_mat,2,function(y){ all(x==t(y))}))
  }),2,NA)
  states$Clone <- apply(states[,1:num_mutations],1,function(x){paste(x,collapse="_",sep="_")})
  
  state_interactions$possible<-ifelse(apply(state_interactions,1,function(x){
    A<-as.numeric(do.call(cbind,strsplit(as.character(x[1]),split="_")))
    B<-as.numeric(do.call(cbind,strsplit(as.character(x[2]),split="_")))
    sum(abs(A-B))<=1
  }),2,NA)
  
  dat<-acast(Var1~Var2,value.var="possible",data=data.frame(state_interactions,
                                                            "value"=-1))
  diag(dat)<-0
  lowerTriangle(dat)<-NA
  
  test_dat<-do.call(cbind,lapply(colnames(dat),function(Clone){
    if(Clone%in%names(weights)){
      output<-ifelse(dat[,Clone]==2,weights[Clone],dat[,Clone])
    } else{
      output<-dat[,Clone]
    }
    return(output)
  }))
  colnames(test_dat) <-rownames(test_dat)
  graph<-graph.adjacency(test_dat,mode="directed",weighted=TRUE)
  graph_mat <- get.data.frame(graph)%>% drop_na()
  
  subgraphs<-make_ego_graph(graph_from_data_frame(graph_mat,directed=TRUE), order=1,nodes=names(weights), mode="all")
  subgraph_bind<-do.call(rbind,lapply(subgraphs,get.data.frame)) %>% distinct(to,from,weight, .keep_all = TRUE)
  subgraph_subset<-subgraph_bind%>%filter(!to%in%setdiff(setdiff(subgraph_bind$to,subgraph_bind$from),names(weights)))
  return(subgraph_subset)
}


create_reward_matrix<-function(Known_mat,weights){
  
  set.seed(68864)
  names(weights) <- apply(Known_mat,2,function(x){paste(x,sep="_",collapse="_")})
  num_type <- 2
  num_mutations <- nrow(Known_mat); 
  mutant_names<-rownames(Known_mat)
  num_clones <- ncol(Known_mat)
  num_states <- num_type^num_mutations
  
  possible_mut_list<- unlist(apply(Known_mat,1,function(x){list(0:max(unique(x))) }),recursive = FALSE)
  
  states<-data.frame(expand.grid(possible_mut_list))
  state_interactions<-data.frame(expand.grid(apply(states[,1:num_mutations],1,function(x){paste(x,collapse="_",sep="_")}),
                                             apply(states[,1:num_mutations],1,function(x){paste(x,collapse="_",sep="_")})))
  
  state_interactions$possible<-ifelse(apply(state_interactions,1,function(x){
    A<-as.numeric(do.call(cbind,strsplit(as.character(x[1]),split="_")))
    B<-as.numeric(do.call(cbind,strsplit(as.character(x[2]),split="_")))
    sum(abs(A-B))<=1
  }),0,NA)
  
  state_interactions$action<-apply(state_interactions,1,function(x){
    A<-as.numeric(do.call(cbind,strsplit(as.character(x[1]),split="_")))
    B<-as.numeric(do.call(cbind,strsplit(as.character(x[2]),split="_")))
    if(!is.na(x["possible"])){
      if(sum(abs(B-A))==0){
        return("stay")
      } else{
        return(mutant_names[which((B-A)==1)])
      }
    }
  })
  
  dat<-setNames(state_interactions%>%filter(action%in%c(mutant_names,"stay")),
                c("State","NextState","Reward","Action"))[,c(1,4,2,3)]
  
  dat$Reward <- as.numeric(apply(dat,1,function(x){
    ifelse(x$NextState%in%names(weights),weights[x$NextState],x$Reward)
  }))
  dat$Reward <- as.numeric(apply(dat,1,function(x){
    ifelse(x$Action%in%"stay",0,x$Reward)
  }))
  dat$State <- as.character(dat$State)
  dat$NextState <- as.character(dat$NextState)
  dat$Action <- as.character(dat$Action)
  
  control <- list(alpha = 0.8, gamma = 0.9)
  model <- ReinforcementLearning(data = dat, s = "State", a = "Action", r = "Reward",  s_new = "NextState",  iter =  1,control=control)
  x<- model$Q
  rownames(x) <- substring(rownames(x),2)
  Q_mat <- setNames(melt(x),c("State","Action","Q"))
  set<-inner_join(dat,Q_mat,by=c("State","Action"))
  set$Valid <- TRUE
  return(set)
  }


plot_clonal_heatmap_and_barplot_with_SD <-function(sample,output_folder,clonal_abundance,clonal_architecture,clone_cell_count,shrink)
{
  clonal_abundance_subset <-data.frame( clonal_abundance[[sample]])#%>%filter(Count>=5)
  clones_to_use <- intersect(as.character(clonal_abundance_subset$Clone),as.character(clonal_architecture[[sample]]$Clone))
  clonal_architecture_subset <- data.frame(clonal_architecture[[sample]])%>%
    filter(data.frame(clonal_architecture[[sample]])$Clone%in%clones_to_use)
  clonal_abundance_subset <- data.frame(clonal_abundance_subset)%>%
    filter(Clone%in%clones_to_use)
  genes_to_display <-unique(as.character(clonal_architecture_subset[clonal_architecture_subset$Genotype%in%c("Heterozygous","Homozygous"),"Mutant"]))
  
  if(shrink==TRUE){
    clonal_architecture_subset <-clonal_architecture_subset%>%filter(Mutant%in%genes_to_display)
  }
  clonal_architecture_subset$Clone <- factor(clonal_architecture_subset$Clone, levels=rev(clonal_abundance_subset$Clone))
  clonal_abundance_subset$Clone <- factor(clonal_abundance_subset$Clone, levels=levels(clonal_architecture_subset$Clone))
  gg_heatmap <- ggplot(data = clonal_architecture_subset, aes(x = Clone, y = factor(Mutant,levels=rev(levels(factor(Mutant)))), fill = Genotype)) + 
    geom_tile() +# scale_y_discrete(limits = rev(levels(Mutant)))+
    
    scale_fill_manual(values=c("WT"=brewer.pal(7,"Reds")[1],
                               "Heterozygous"=brewer.pal(7,"Reds")[3],
                               "Homozygous"=brewer.pal(7,"Reds")[6],
                               "Unknown"="grey50"),"Genotype")  +
    theme_classic(base_size=6) +
    ylab("Mutation")+
    theme(legend.position = "bottom", legend.direction = "horizontal",
          axis.text.x = element_blank(), 
          axis.line=element_blank(),
          axis.title.x=element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin=unit(c(0,1,1,1),"cm"))
  
  # Generate clonal abundance barplot
  gg_clonal_barplot <- ggplot(data = data.frame(clonal_abundance_subset), aes(x = Clone, y = Count)) + 
    geom_bar(stat = "identity", aes(fill = Count)) + theme_gray() +
    theme_classic(base_size=6)+
    ylim(0,max(clonal_abundance_subset$Count)*1.3) + 
    ylab("Cell Count")+
    geom_errorbar(aes(ymin = X2.5., ymax = X97.5.), width = 0.2)+
    geom_text(aes(label=Count), position=position_dodge(width=0.9), vjust=-0.25,size=1)+
    scale_fill_distiller(name = "Value", palette = "Reds", direction = 1) +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(),  axis.ticks.x = element_blank(), axis.line=element_blank(),
          legend.position = "none",
          plot.margin=unit(c(1,1,-0.2,1),"cm"))
  
  grob.title <- textGrob(paste(sample,"-"), hjust= 0.5, vjust = 0.5, gp = gpar(fontsize = 12))
  final_plot<-plot_grid(grob.title,gg_clonal_barplot,gg_heatmap,ncol=1,align="hv",axis="l",rel_heights = c(0.1,1,0.75))
  
  save_plot(paste(output_folder,sample,".pdf",sep=""),final_plot, ncol=1) # Open a new pdf file
}

query_initiating_mutations<-function(graph_results){
  start_index<-paste(rep(0,length(strsplit(graph_results$State[1],split="_")[[1]])),sep="_",collapse="_")
  possible_starting_actions<-graph_results%>%filter(State==start_index&Action!="stay")%>%pull(Action)
  final_results<-list()
  initating_action_count<-0
  for(initating_action in possible_starting_actions){
    print(initating_action)
    set <- graph_results
    initating_action_count<-initating_action_count+1
    storage_results<- list()
    branches<-0
    state_to_kill <- set%>%filter(State==start_index&Action==initating_action)%>%pull(NextState)
    start_killed <- sum(set%>%filter(State==state_to_kill)%>%pull(Valid))
    while(start_killed>0){
      #print(branches)
     # print(start_killed)
      branches <- branches +1
      number_of_mutations<-0
      state_log<- list()
      optimal_reward<-list()
      action_log<-list()
      current_state<- start_index
      indicator<-TRUE
      nextState<-0
      while(current_state!=nextState)  {
       # print(number_of_mutations)
        number_of_mutations <- number_of_mutations+1
        if(number_of_mutations==1){
          state_log[[number_of_mutations]] <- start_index
        }
        current_state  <- state_log[[number_of_mutations]]
        nextState_indicator<- FALSE
        
        while(nextState_indicator==FALSE){
          
          if(number_of_mutations==1){
            max_potential_action_index<-  set%>%
              filter(State==current_state&Action==initating_action)
          } else {
            max_potential_action_index <- set%>%
              filter(State==current_state&Valid==TRUE)%>%
              filter(Q==max(Q))%>%sample_n(1)
          }
          if(nrow(max_potential_action_index)==0){
            break
          }
          max_potential_action <- max_potential_action_index%>%pull(NextState)
          next_valid_action <- any(set%>%filter(State==max_potential_action&Action!="stay")%>%pull(Valid))  
          if(next_valid_action==TRUE){
            nextState <-max_potential_action
            current_action <-  max_potential_action_index%>%pull(Action)
            nextState_indicator==TRUE
            break
          } else{
            set[set$State%in%max_potential_action_index["State"]&
                  set$Action%in%max_potential_action_index["Action"],"Valid"] <- FALSE  
          }
        }
        if(nrow(set%>%filter(State==current_state&Action==current_action))==0){
          optimal_reward[[number_of_mutations]] <-NA
        } else {
          optimal_reward[[number_of_mutations]] <- set%>%
            filter(State==current_state&Action==current_action)%>%
            pull(Reward) 
        }
        state_log[[number_of_mutations+1]]<- nextState
        action_log[[number_of_mutations]] <- current_action 
        if(current_action==nextState){
          indicator==FALSE
          state_log[[number_of_mutations+1]]<-NULL
          break
        }
      }
      optimal_reward[[number_of_mutations+1]] <- NA
      action_log[[number_of_mutations+1]] <- NA
      storage_results[[branches]] <-data.frame("states"=do.call(rbind,state_log),#[1:(length(state_log)-1)]),
                                               "actions"=do.call(rbind,action_log),
                                               "reward"=do.call(rbind,optimal_reward),
                                               "nextState"=do.call(rbind,c(state_log[2:length(state_log)],NA)) )
      storage_results[[branches]] <- storage_results[[branches]]%>%
        filter(states!=nextState)
      storage_results[[branches]]$cumulative_reward <- cumsum(storage_results[[branches]]$reward)
      
      #storage_results[[branches]] <-storage_results[[branches]][1:which.max(storage_results[[branches]]$cumulative_reward), ]
      set[set$State%in%current_state&set$Action%in%current_action,"Valid"] <- FALSE
      start_killed <- sum(set%>%filter(State==state_to_kill)%>%pull(Valid))
    }
    final_results[[initating_action_count]]<-storage_results[!duplicated(storage_results)]
  }
  names(final_results)<-possible_starting_actions
  return(final_results)
}

