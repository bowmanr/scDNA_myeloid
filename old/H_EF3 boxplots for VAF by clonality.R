source("/Volumes/LevineLab-1/Levine Lab/Bobby/Collaborations/MissionBio/Scripts/mBio_processing_scripts.R")


final_sample_summary<-readRDS(file="/Volumes/LevineLab-1/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/final_sample_summary.rds")
clonal_sample_set_after_boostrap<-readRDS(file="/Volumes/LevineLab-1/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/clonal_sample_set_after_boostrap")


clone_size_by_gene<- lapply(final_sample_summary[clonal_sample_set_after_boostrap],function(x){
  clones<-x$Clones[,2]
  weights<-x$Clones[,1]/sum(x$Clones[,1])
  names(weights) <- clones
  dedup<-x$NGT[!duplicated(x$NGT)&x$NGT[,"Clone"]%in%clones,]
  weights <- weights[dedup$Clone]
  z<- apply(dedup[,1:(ncol(dedup)-1)],2,function(y){
    if(!all(y==0)){max(weights[y!=0])} 
  })
  do.call(cbind,list("Max_Clone_Size"=z,
                     "Clonality"=lapply(z,function(t){
                       ifelse(t==max(weights),"Dominant","Subclone")
                     }),
                     "VAF"=colSums(x$NGT[,1:(ncol(x$NGT)-1)])/(nrow(x$NGT)*2)))
  
})


clone_size_by_gene <- lapply(names(clone_size_by_gene),function(x){
  data.frame(clone_size_by_gene[[x]],"Sample"=x)
})

clone_size_mat<-setNames(data.frame(do.call(rbind,clone_size_by_gene)),c("Clone_size","Clonality","VAF","Sample"))
clone_size_mat<- data.frame(clone_size_mat,
                            "Variant"=rownames(clone_size_mat),
                            "Gene"=as.character(do.call(rbind,strsplit(rownames(clone_size_mat),split="\\."))[,1]))%>%
  filter(Clonality%in%c("Dominant","Subclone"))
clone_size_mat$Clone_size <- as.numeric(unlist(clone_size_mat$Clone_size))
clone_size_mat$VAF <- as.numeric(unlist(clone_size_mat$VAF))
clone_size_mat$Clonality <- factor(clone_size_mat$Clonality,levels=c("Subclone","Dominant"))
clone_size_mat$Gene <- as.character(clone_size_mat$Gene)
clone_size_mat$Gene <- factor(do.call(rbind,strsplit(
  ifelse(clone_size_mat$Gene=="FLT3_INS_chr13","FLT3",
         ifelse(clone_size_mat$Gene=="FLT3_INS_","FLT3",clone_size_mat$Gene)),
  split="_c"))[,1])

tally_set<-table(clone_size_mat$Gene,as.character(clone_size_mat$Clonality))
clone_size_mat$Gene <- factor(clone_size_mat$Gene,levels=rev(rownames(tally_set)[order(tally_set[,2]/(tally_set[,1]+tally_set[,2]))]))

clone_size_mat$Group <- ifelse(grepl("^CH",clone_size_mat$Sample),"CH",
                               ifelse(grepl("^MA",clone_size_mat$Sample),"Post MPN",
                                      ifelse(grepl("^P",clone_size_mat$Sample),"TP53",
                                             ifelse(grepl("^CB",clone_size_mat$Sample),"CBF",
                                                    ifelse(grepl("^E",clone_size_mat$Sample),"Epigenetic",
                                                           ifelse(grepl("^S",clone_size_mat$Sample),"Epigenetic",
                                                                  ifelse(grepl("^CP",clone_size_mat$Sample),"Epigenetic",
                                                                         ifelse(grepl("^F",clone_size_mat$Sample),"FLT3-ITD",
                                                                                ifelse(grepl("R",clone_size_mat$Sample),"RAS","Error")))))))))


clone_size_mat$Group <-factor(clone_size_mat$Group, levels=c("CH","Post MPN","CBF","Epigenetic","TP53","RAS","FLT3-ITD"))

clone_size_mat[clone_size_mat$Sample%in%c("S5346","S5941","R1551A","CP0417d"),"Group"] <- "FLT3-ITD"
clone_size_mat[clone_size_mat$Sample%in%c("17-005","17-006","Meyer_1"),"Group"]<- "Epigenetic"
clone_size_mat[clone_size_mat$Sample%in%c("17-021","17-020"),"Group"]<- "RAS"
clone_size_mat[clone_size_mat$Sample%in%c("CP0417d","CP8045d","CP0289d","CP7964d","CP1703d","CP5853d"),"Group"]<- "RAS"
clone_size_mat[clone_size_mat$Sample%in%c("CP6408d","CP5764d"),"Group"]<- "FLT3-ITD"

clone_size_mat%>%filter(Group=="Error")


mutants_in_each_sample<-do.call(rbind,lapply(final_sample_summary,function(x){
  y<-colnames(x$NGT)
  z <- list()
  z$TET2 <- ifelse(any(grepl("TET2",y)),1,0)
  z$DNMT3A <- ifelse(any(grepl("DNMT3A",y)),1,0)
  z$DNMT3A_R882 <- ifelse(any(grepl("DNMT3A.p.R882",y)),1,0)
  z$ASXL1 <- ifelse(any(grepl("ASXL1",y)),1,0)
  z$IDH <- ifelse(any(grepl("IDH",y)),1,0)
  z$RAS <- ifelse(any(grepl("RAS",y)),1,0)
  z$FLT3_ITD <-ifelse(any(grepl("FLT3_INS",y)),1,0)
  z$FLT3 <- ifelse(any(grepl("FLT3.p",y)),1,0)
  z$JAK2 <- ifelse(any(grepl("JAK2",y)),1,0)
  data.frame(t(do.call(rbind,z)))
}))
rownames(mutants_in_each_sample) <- names(final_sample_summary)
mutants_in_each_sample$Sample <- rownames(mutants_in_each_sample)
mutants_in_each_sample$DTAI <- with(mutants_in_each_sample,ifelse((TET2==1|DNMT3A==1|IDH==1|ASXL1==1)&(RAS==0&FLT3==0&FLT3_ITD==0),1,0 ))
mutants_in_each_sample$DTAI_RAS <- with(mutants_in_each_sample,ifelse((TET2==1|DNMT3A==1|IDH==1|ASXL1==1)&(RAS==1&FLT3==0&FLT3_ITD==0),1,0 ))
mutants_in_each_sample$DTAI_FLT3 <- with(mutants_in_each_sample,ifelse((TET2==1|DNMT3A==1|IDH==1|ASXL1==1)&(RAS==0&(FLT3==1|FLT3_ITD==1)),1,0 ))
mutants_in_each_sample$DTAI_FLT3_RAS <- with(mutants_in_each_sample,ifelse((TET2==1|DNMT3A==1|IDH==1|ASXL1==1)&(RAS==1&(FLT3==1|FLT3_ITD==1)),1,0 ))
mutants_in_each_sample$Signaling <- with(mutants_in_each_sample,ifelse((TET2==0&DNMT3A==0&IDH==0&ASXL1==0)&(RAS==1|FLT3==1|FLT3_ITD==1|JAK2==1),1,0 ))


color_red<-brewer.pal(5,"Reds")[5]

test<-merge(clone_size_mat,mutants_in_each_sample,by="Sample")



test$Group2 <- ifelse(test$Group=="CH","CH",
                      ifelse(test$DTAI==1&test$Group!="CH","DTAI",
                             ifelse(test$DTAI_RAS ==1,"DTAI-RAS",
                                    ifelse(test$DTAI_FLT3 ==1,"DTAI-FLT3",
                                           ifelse(test$Signaling==1,"JAK/RAS/FLT3",
                                                  ifelse(test$DTAI_FLT3_RAS==1,"DTAI-FLT3-RAS","Other"))))))

test[test$Sample%in%c("MA1715B","MA0092A","MA1715B","MA4244A","MA6300A","MA6363B","MA6363C","MA9521A","MA9521B"),"Group2"]<- "MPN"

test$Group <-factor(test$Group, levels=c("CH","Post MPN","CBF","Epigenetic","TP53","RAS","FLT3-ITD"))
test$Group2 <-factor(test$Group2, levels=c("CH","MPN","JAK/RAS/FLT3","DTAI","DTAI-RAS","DTAI-FLT3","DTAI-FLT3-RAS","Other"))
test%>%filter(as.character(Group2)=="Other")

test%>%filter(Group2=="DTAI")%>%filter(Dominant_clone_size<0.6)%>%filter(Number_of_clones<15)%>%filter(Number_of_mutations_in_dominant_clone>2)%>%pull(Sample)

test[test$Sample%in%c("S5346","S5941","CP0417d","R1551"),"Group"] <- "FLT3-ITD"
test[test$Sample%in%c("17-005","17-006","Meyer_1"),"Group"]<- "Epigenetic"
test[test$Sample%in%c("17-021","17-020"),"Group"]<- "RAS"
test[test$Sample%in%c("CP0417d","CP8045d","CP0289d","CP7964d","CP1703d","CP5853d"),"Group"]<- "RAS"
test[test$Sample%in%c("CP6408d","CP5764d"),"Group"]<- "FLT3-ITD"


test%>%filter(Group=="Error")

gene_order <- c("DNMT3A","TET2","ASXL1","IDH1","IDH2","JAK2","NRAS","KRAS","FLT3","NPM1")


test%>%filter(Clonality=="Dominant"&VAF<.2)

summarized_test<-data.frame(test%>%filter(Gene%in%gene_order&!Group2%in%c("CH","MPN"))%>%group_by(Gene,Clonality)%>%summarise(mean=mean(VAF),
                                                                                             sd = sd(VAF),
                                                                                             sem = sd(VAF)/sqrt(length(VAF))))


clonality_VAF<-ggplot(test%>%filter(Gene%in%gene_order&!Group2%in%c("CH","MPN")),aes(x=Clonality,y=VAF,color=Clonality))+
                    facet_wrap(~factor(Gene,levels=gene_order),scale="free_x",ncol=5)+
                    ggbeeswarm::geom_beeswarm()+xlab("")+
                    geom_errorbar(data=summarized_test,aes(x=Clonality,y=mean,ymin=mean-sem,ymax=mean+sem),color="black")+
                    scale_color_manual(values=c("Dominant"=color_red,"Subclone"="grey50"))+
                    theme_classic()+guides(fill=FALSE)+
                    theme(axis.ticks.x = element_blank(),
                          axis.text.x = element_blank())+
                    scale_y_continuous(limits=c(0,1.1),breaks=c(0,.25,.5,.75,1),labels=c("0","0.25","0.50","0.75","1.0"))+ylab("Computed VAF")

ggsave(clonality_VAF, width=6,height=4.25,
       file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/Figures/EF3D-clonality_VAF.pdf")


gene_order2 <- c("DNMT3A","TET2","ASXL1","IDH1","NRAS","KRAS","FLT3","NPM1")

EF3D_clonality_VAF_pvalues<-data.frame(test%>%filter(Gene%in%gene_order2&!Group2%in%c("CH","MPN"))%>%group_by(Gene)%>%summarise(p_value= t.test(VAF~Clonality)$p.value))
write.table(EF3D_clonality_VAF_pvalues,
            file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/Stats/EF3D_clonality_VAF.txt",sep="\t",row.names=FALSE,quote=FALSE)

