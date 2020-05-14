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


color_red<-brewer.pal(5,"Reds")[5]

ggA<-ggplot(tally(clone_size_mat%>%group_by(Gene,Group,Clonality)),
            aes(x=factor(Gene),fill=Clonality,y=n,label=n)) +guides(fill=FALSE,color=FALSE)+
  scale_y_continuous( expand = c(0, 0.0))+

  geom_bar(stat="identity",position="fill")+xlab("")+coord_flip()+
scale_fill_manual(values=c("Dominant"=color_red,"Subclone"="grey80"))+
  ylab("Fraction of mutant samples \n with mutation in dominant clone")+theme_bw(base_size=8)+theme(legend.position = "bottom")
  


ggB<-ggplot(clone_size_mat, aes(y = Clone_size, x =Gene, fill = Gene)) +
  geom_boxplot(alpha = 0.5,outlier.shape = NA)+geom_point(aes(color=Clonality,group=Clonality),position = position_jitterdodge(),size=0.3)+
  scale_fill_manual(values=tol.rainbow(n=length(levels(clone_size_mat$Gene))))+
  scale_color_manual(values=c("Dominant"=color_red,"Subclone"="grey20"))+
  coord_flip()+theme_bw(base_size=8)+guides(fill=FALSE,color=FALSE)+
  theme(axis.text.y = element_blank(),
        axis.ticks.y  = element_blank(),
        axis.title.y = element_blank())+
  scale_y_continuous(limits = c(0,1), expand = c(0, 0.05)) +
  ylab("Fraction of cells \n in largest mutant clone")+
  theme(legend.position = "bottom")



spacer <- plot_grid(NULL)


ggsave(plot_grid(plot_grid(ggA),plot_grid(NULL),plot_grid(ggB),align="h",axis="tb",ncol=3,rel_widths=c(1,0.05,1)), width=4.5,height=4.25,
       file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/Figures/F2C-clone_size_by_gene.pdf")






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






test<-merge(mutants_in_each_sample,clone_size_mat,by="Sample")

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


test[test$Sample%in%c("S5346","S5941","CP0417d","R1551"),"Group"] <- "FLT3-ITD"
test[test$Sample%in%c("17-005","17-006","Meyer_1"),"Group"]<- "Epigenetic"
test[test$Sample%in%c("17-021","17-020"),"Group"]<- "RAS"
test[test$Sample%in%c("CP0417d","CP8045d","CP0289d","CP7964d","CP1703d","CP5853d"),"Group"]<- "RAS"
test[test$Sample%in%c("CP6408d","CP5764d"),"Group"]<- "FLT3-ITD"


test%>%filter(Group=="Error")




genes_of_interest <- c("ASXL1","TET2","DNMT3A","DNMT3A.p.R882","IDH1","IDH2")

#test$Gene <- ifelse(grepl("IDH",as.character(test$Gene)),"IDH",as.character(test$Gene))
test$Gene <- ifelse(grepl("DNMT3A.p.R882",as.character(test$Variant)),"DNMT3A.p.R882",as.character(test$Gene))

mutation_by_dominance<-ggplot(tally(test%>%filter(Gene%in%genes_of_interest)%>%
                    group_by(Gene,Group2,Clonality)),
            aes(x=factor(Group2),#,levels=c("ASXL1","IDH1","IDH2","DNMT3A","DNMT3A R882","TET2","NPM1")),
                fill=Clonality,y=n)) +facet_wrap(~factor(Gene,levels=c("DNMT3A","TET2","ASXL1","DNMT3A.p.R882","IDH1","IDH2")),ncol=3)+
  geom_bar(stat="identity",position="stack",alpha=1)+xlab("")+
  scale_fill_manual(values=c("Dominant"=color_red,"Subclone"="grey80"))+
  ylab("Number of samples")+theme_bw(base_size=10)+theme(legend.position = "right")+  theme(axis.text.x =element_text(angle=30,hjust=1))


ggsave(mutation_by_dominance, width=6,height=3,
       file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/Figures/EF3C-mutation_by_patient_group_clonality.pdf")


