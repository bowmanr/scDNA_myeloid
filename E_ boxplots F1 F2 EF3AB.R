source("/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Scripts/mBio_processing_scripts.R")


final_sample_summary<-readRDS(file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/final_sample_summary.rds")
clonal_sample_set_after_boostrap<-readRDS(file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/clonal_sample_set_after_boostrap")

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

clonal_level_info<-data.frame(do.call(rbind,lapply(names(final_sample_summary),function(y){
  x <- final_sample_summary[[y]]$Clones
  data.frame("Shannon"=vegan::diversity(x[,1],index="shannon"),
             "Number_of_clones"=length(x[,1]),
             "Number_of_mutations"=dim(do.call(rbind,strsplit(as.character(x[,2]),split="_")))[[2]],
             "Number_of_mutations_in_dominant_clone"=sum(as.numeric(do.call(rbind,strsplit(as.character(x[dim(x)[1],2]),split="_")))),
             "Dominant_clone_size"=max(x[,1]/sum(x[,1])) ,
             "Second_clone_size"=(x[,1]/sum(x[,1]))[2],
             "Sample"=y,
             "Number_of_large_clones"=sum((x[,1]/sum(x[,1]))>0.1))
})))

gene_level_info <- data.frame(do.call(rbind,lapply(names(final_sample_summary),function(y){
  x <- final_sample_summary[[y]]$NGT
  data.frame("Max_VAF"=max(colSums(x[,!c(colnames(x)=="Clone")])/(length(rownames(x))*2)),
             "Sample"=y,
             "Gene_Shannon"=vegan::diversity(colSums(x[,!c(colnames(x)=="Clone")]),index="shannon"))
})))

test<-merge(mutants_in_each_sample,merge(clonal_level_info,gene_level_info,by="Sample"),by="Sample")

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


test%>%filter(Group=="Other")


gg_number_of_mutations<-ggplot(test%>%group_by(Group2)%>%summarise(mean=mean(Number_of_mutations),
                                    sd = sd(Number_of_mutations),
                                    sem = sd(Number_of_mutations)/sqrt(length(Number_of_mutations))),
      aes(x=Group2,y=mean,fill=Group2))+
      geom_bar(stat="identity",color="black")+
       geom_errorbar(aes(ymin = mean-sem, ymax = mean+sem),width=0.5,lwd=0.5)+
  theme_classic(base_size = 8)+
      ylab("Number of mutations")+xlab("")+ggtitle("")+
  scale_y_continuous(limits = c(0,9), expand = c(0, 0)) +
  theme(axis.text.x = element_text(angle=30,hjust=1)) +
  scale_fill_brewer(type="seq",palette = "Reds",aesthetics = "fill",guide=FALSE)


ggsave(gg_number_of_mutations, width=2.5,height=2.25,
       file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/Figures/F1C-number_of_mutations.pdf")


gg_number_of_clones<-ggplot(test,aes(y=Number_of_clones,x=Group2,fill=Group2))+
                    geom_boxplot(outlier.shape = NA)+  
                      geom_jitter(width = 0.1,size=0.5)+
  theme_classic(base_size = 8)+
  ylab("Number of clones")+
                    xlab("")+
                    theme(axis.text.x = element_text(angle=30,hjust=1)) +
                    scale_fill_brewer(type="seq",palette = "Reds",aesthetics = "fill",guide=FALSE)
  
ggsave(gg_number_of_clones, width=2.5,height=2.25,
       file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/Figures/F1E-number_of_clones.pdf")



gg_Number_of_mutations_in_Dclone<-ggplot(test%>%group_by(Group2)%>%summarise(mean=mean(Number_of_mutations_in_dominant_clone),
                                                                             sd = sd(Number_of_mutations_in_dominant_clone),
                                                                             sem = sd(Number_of_mutations_in_dominant_clone)/sqrt(length(Number_of_mutations_in_dominant_clone))),
                                         aes(x=Group2,y=mean,fill=Group2))+
  geom_bar(stat="identity",color="black")+
  geom_errorbar(aes(ymin = mean-sem, ymax = mean+sem),width=0.5,lwd=0.5)+
  theme_classic(base_size = 8)+
  ylab("Number of mutations \n in dominant clone")+xlab("")+ggtitle("")+
  scale_y_continuous(limits = c(0,4.5), expand = c(0, 0)) +
  theme(axis.text.x = element_text(angle=30,hjust=1)) +
  scale_fill_brewer(type="seq",palette = "Reds",aesthetics = "fill",guide=FALSE)

ggsave(gg_Number_of_mutations_in_Dclone, width=2.5,height=2.25,
       file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/Figures/EF3-number_of_mutations_in_dominant_clone.pdf")



gg_shannon<-ggplot(test,aes(y=Shannon,x=Group2,fill=Group2))+
  geom_boxplot(outlier.shape = NA)+  
  geom_jitter(width = 0.1,size=0.5)+
  theme_classic(base_size = 8)+
  ylab("Shannon diveristy index")+
  xlab("")+
  theme(axis.text.x = element_text(angle=30,hjust=1)) +
  scale_fill_brewer(type="seq",palette = "Reds",aesthetics = "fill",guide=FALSE)

ggsave(gg_shannon, width=2.5,height=2.25,
       file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/Figures/F2A-shannon.pdf")




gg_dominant_clone_size<-ggplot(test,aes(y=Dominant_clone_size,x=Group2,fill=Group2))+
  geom_boxplot(outlier.shape = NA)+  
  geom_jitter(width = 0.1,size=0.5)+
  theme_classic(base_size = 8)+
  ylab("Fraction of sample \n in dominant clone")+
  xlab("")+
  theme(axis.text.x = element_text(angle=30,hjust=1)) +
  scale_fill_brewer(type="seq",palette = "Reds",aesthetics = "fill",guide=FALSE)

ggsave(gg_dominant_clone_size, width=2.5,height=2.25,
       file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/Figures/F2B-percent_dominant_clone.pdf")



ggsave(plot_grid(gg_shannon,gg_dominant_clone_size,ncol=1,align="v",axis="lr"), width=2.5,height=4.25,
       file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/Figures/F2AB-shannon-percent_dominant_clone.pdf")


pvalues_Dominant_clone_size<-test%>%{melt(pairwise.t.test(.$Dominant_clone_size,g=.$Group2,data=.,p.adjust.method="fdr")$p.value)}%>%filter(!is.na(value))%>%filter(value<0.1)
pvalues_Number_of_clones<-test%>%{melt(pairwise.t.test(.$Number_of_clones,g=.$Group2,data=.,p.adjust.method="fdr")$p.value)}%>%filter(!is.na(value))%>%filter(value<0.1)
pvalues_Number_of_mutations<-test%>%{melt(pairwise.t.test(.$Number_of_mutations,g=.$Group2,data=.,p.adjust.method="fdr")$p.value)}%>%filter(!is.na(value))%>%filter(value<0.1)
pvalues_Shannon<-test%>%{melt(pairwise.t.test(.$Shannon,g=.$Group2,data=.,p.adjust.method="fdr")$p.value)}%>%filter(!is.na(value))%>%filter(value<0.1)
pvalues_Number_of_mutations_in_dominant_clone<-test%>%{melt(pairwise.t.test(.$Number_of_mutations_in_dominant_clone,g=.$Group2,data=.,p.adjust.method="fdr")$p.value)}%>%filter(!is.na(value))%>%filter(value<0.1)

write.table(pvalues_Number_of_mutations_in_dominant_clone ,file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/Stats/EF3-number_of_mutations_in_dominant_clone.txt",sep="\t",row.names = FALSE)
write.table(pvalues_Dominant_clone_size ,file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/Stats/F2B-percent_dominant_clone.txt",sep="\t",row.names = FALSE)
write.table(pvalues_Shannon ,file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/Stats/F2A-shannon.txt",sep="\t",row.names = FALSE)
write.table(pvalues_Number_of_mutations ,file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/Stats/F1C-number_of_mutations.txt",sep="\t",row.names = FALSE)
write.table(pvalues_Number_of_clones ,file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/Stats/F1E-number_of_clones.txt",sep="\t",row.names = FALSE)



clone_size_by_genetic_density<- do.call(rbind,lapply(final_sample_summary[clonal_sample_set_after_boostrap],function(x){
  possible_clones_subset <-x$Clones%>%filter(Clone%in% x$Clones[,"Clone"] )
  clones<-possible_clones_subset[,"Clone"]
  dedup<-x$NGT[!duplicated(x$NGT)&x$NGT[,"Clone"]%in%clones,]
  set_mat<-full_join(possible_clones_subset[,1:2],dedup)
  counts <-set_mat[,"Count"]
  weights<-set_mat[,"Count"]/sum(set_mat[,"Count"])
  genetic_complexity <- rowSums(set_mat[,-c(1:2)])
  return(data.frame("Clone_size"=weights,
                    "Genetic_density"=genetic_complexity))
  
}))


gg_clone_size_by_genetic_density<-ggplot(clone_size_by_genetic_density,aes(y=Clone_size,x=factor(Genetic_density),fill=factor(Genetic_density)))+
  geom_jitter(width = 0.1,size=0.5)+
  geom_boxplot(outlier.shape = NA)+  
  
  theme_bw(base_size = 8)+
  ylab("Fraction of sample in clone")+
  xlab("Number of mutant alleles")+
  scale_fill_brewer(type="seq",palette = "Greens",aesthetics = "fill",guide=FALSE)

ggsave(plot_grid(gg_Number_of_mutations_in_Dclone,gg_clone_size_by_genetic_density,align="hv",axis="tb",ncol=2), 
            width=5.5,height=2.25,
       file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/Figures/EF3AB-dominant_clone_number-clone_size_by_genetic_density.pdf")

