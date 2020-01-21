source("/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Scripts/mBio_processing_scripts.R")
final_sample_summary<-readRDS(file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/final_sample_summary.rds")
clonal_sample_set_after_boostrap<-readRDS(file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/clonal_sample_set_after_boostrap")

WT_cells<-data.frame("WT"=do.call(rbind,setNames(lapply(final_sample_summary[clonal_sample_set_after_boostrap],function(x){
  clones<-x$Clones[,2]
  weights<-x$Clones[,1]/sum(x$Clones[,1])
  names(weights) <- clones
  dedup<-x$NGT[!duplicated(x$NGT)&x$NGT[,"Clone"]%in%clones,]
  weights <- weights[dedup$Clone]
  
  if(any(rowSums(dedup[,1:(ncol(dedup)-1)])==0)){
    return(weights[rowSums(dedup[,1:(ncol(dedup)-1)])==0 ])
  } else{
    return(0)
  }
}),clonal_sample_set_after_boostrap)),"Sample"=clonal_sample_set_after_boostrap)

colnames(WT_cells)[1]<-"WT"

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

test<-merge(mutants_in_each_sample,WT_cells,by="Sample")

test$Group <- ifelse(grepl("^CH",test$Sample),"CH",
                     ifelse(grepl("^MA",test$Sample),"Post MPN",
                            ifelse(grepl("^P",test$Sample),"TP53",
                                   ifelse(grepl("^CB",test$Sample),"CBF",
                                          ifelse(grepl("^E",test$Sample),"Epigenetic",
                                                 ifelse(grepl("^S",test$Sample),"Epigenetic",
                                                        ifelse(grepl("^CP",test$Sample),"Epigenetic",
                                                               ifelse(grepl("^F",test$Sample),"FLT3-ITD",
                                                                      ifelse(grepl("R",test$Sample),"RAS","Error")))))))))


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
test[test$Sample%in%c("CB6743"),"Group2"] <- "JAK/RAS/FLT3"

test%>%filter(Group2=="DTAI")%>%filter(Dominant_clone_size<0.6)%>%filter(Number_of_clones<15)%>%filter(Number_of_mutations_in_dominant_clone>2)%>%pull(Sample)

test[test$Sample%in%c("S5346","S5941","CP0417d","R1551"),"Group"] <- "FLT3-ITD"
test[test$Sample%in%c("17-005","17-006","Meyer_1"),"Group"]<- "Epigenetic"
test[test$Sample%in%c("17-021","17-020"),"Group"]<- "RAS"
test[test$Sample%in%c("CP0417d","CP8045d","CP0289d","CP7964d","CP1703d","CP5853d"),"Group"]<- "RAS"
test[test$Sample%in%c("CP6408d","CP5764d"),"Group"]<- "FLT3-ITD"


test%>%filter(Group=="Error")

gg_number_of_WT_cells<-ggplot(test,aes(y=WT,x=Group2,fill=Group2))+
  geom_boxplot(outlier.shape = NA)+  
  geom_jitter(width = 0.1,size=0.5)+
  theme_classic(base_size = 8)+
  ylab("Fraction Non-mutant")+
  xlab("")+
  theme(axis.text.x = element_text(angle=30,hjust=1)) +
  scale_fill_brewer(type="seq",palette = "Reds",aesthetics = "fill",guide=FALSE)

ggsave(gg_number_of_WT_cells, width=2.5,height=2.25,
       file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/Figures/F1X-number_WT_cells.pdf")



pvalues_Number_of_WT_cells<-test%>%{melt(pairwise.t.test(.$WT,g=.$Group2,data=.,p.adjust.method="fdr")$p.value)}%>%filter(!is.na(value))%>%filter(value<0.1)

write.table(pvalues_Number_of_WT_cells ,file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/Stats/F1X-number_of_WT_cells.txt",sep="\t",row.names = FALSE)
