source("/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Scripts/mBio_processing_scripts.R")
library(UpSetR)

final_sample_summary<-readRDS(file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/final_sample_summary.rds")
clonal_sample_set_after_boostrap<-readRDS(file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/clonal_sample_set_after_boostrap")

## Are epigenetic mutations found in the same cell, if so, are these cells the dominant clone?
## Are there Tet2 IDH in the same sample, and if so same cell?
## Are there signaling in the same sample, and if so same cell?


DTAI_genes <- c("ASXL1","DNMT3A","TET2","IDH1","IDH2")

genes_in_each_sample<-setNames(mclapply(mc.cores = 4,final_sample_summary[clonal_sample_set_after_boostrap],function(x){
  aggregate_by<-do.call(rbind,strsplit(colnames(x$NGT),split="[_\\.]")) [,1][-length(colnames(x$NGT))]
}), clonal_sample_set_after_boostrap)

list_of_mutants<-setNames(lapply(DTAI_genes,function(x){
  names(genes_in_each_sample)[do.call(c,setNames(lapply(genes_in_each_sample,function(y){
    x%in%y
  }),names(genes_in_each_sample)))]
}),DTAI_genes)

elements <- unique(unlist(list_of_mutants))
data <- unlist(lapply(list_of_mutants, function(x) {
  x <- as.vector(match(elements, x))
}))
data[is.na(data)] <- as.integer(0)
data[data != 0] <- as.integer(1)
data <- data.frame(matrix(data, ncol = length(list_of_mutants), byrow = F))
data <- data[which(rowSums(data) != 0), ]
names(data) <- names(list_of_mutants)
rownames(data) <-elements
data$Group <- ifelse(grepl("CH",rownames(data)),"CH","AML")


DTAI_in_dominant_clone<-lapply(final_sample_summary[clonal_sample_set_after_boostrap],function(sample){
  variants<-colnames(sample$NGT)[1:(length(colnames(sample$NGT))-1)]
  genes <-do.call(rbind,strsplit(variants,split="[_\\.]"))[,1]
  clones<-sample$Clones[,"Clone"]
  WT_clone_to_exclude <-!apply(do.call(rbind,strsplit(clones,split="_")),1,function(x){all(x=="0")})
  dominant_clone <-sample$Clones%>%filter(!apply(do.call(rbind,strsplit(Clone,split="_")),1,function(x){all(x=="0")})) %>%filter(Count==max(Count))%>%pull(Clone)
  dominant_clone_genes<-genes[do.call(c,strsplit(dominant_clone,split="_"))!=0]
  unique(intersect(dominant_clone_genes,DTAI_genes))
})

list_of_mutants_in_dominant_clone<-setNames(lapply(DTAI_genes,function(x){
  names(DTAI_in_dominant_clone)[do.call(c,setNames(lapply(DTAI_in_dominant_clone,function(y){
    x%in%y
  }),names(DTAI_in_dominant_clone)))]
}),DTAI_genes)



elements_dominant <- unique(unlist(list_of_mutants_in_dominant_clone))
data_dom <- unlist(lapply(list_of_mutants_in_dominant_clone, function(x) {
  x <- as.vector(match(elements_dominant, x))
}))
data_dom[is.na(data_dom)] <- as.integer(0)
data_dom[data_dom != 0] <- as.integer(1)
data_dom <- data.frame(matrix(data_dom, ncol = length(list_of_mutants_in_dominant_clone), byrow = F))
data_dom <- data_dom[which(rowSums(data_dom) != 0), ]
names(data_dom) <- names(list_of_mutants)
rownames(data_dom) <-elements_dominant
data_dom$Group <- ifelse(grepl("CH",rownames(data_dom)),"CH","AML")


rownames(data_dom)%in%rownames(data)

data_dom$Match<-ifelse(sapply(rownames(data_dom),function(sample){
  all(data_dom[sample,]==data[sample,])
}),"Match","Absent")
data_dom$Sample<-rownames(data_dom)
data$Sample<-rownames(data)

test_set<-full_join(data,data_dom[,7:8],by="Sample")
test_set$Match[is.na(test_set$Match)]<-"Absent"

test_set[test_set$Sample%in%c("MA1715B","MA6300A","MA9521A","MA9521B"),"Group"]<- "MPN"
test_set[test_set$Sample%in%c("E4840","E4838new"),"Group"]<- "CMML"
test_set[test_set$Sample%in%c("MA4244A","MA2725","MA9521B","MA0092A","R2715","MA1715B","MA6363B"),"Group"]<- "MF"

Myfunc <- function(row,feature) {
  data <- (row[feature]=="Match")
}

AML<-upset(test_set%>%filter(Group=="AML"), sets=DTAI_genes,order.by = c("degree"),
           main.bar.color = "grey60",decreasing=FALSE,
           mainbar.y.label = "Number of samples",
           sets.x.label = "Number of \n mutant samples",
           text.scale=1.25,
           shade.alpha = 0.75,
           show.numbers=FALSE,
           mb.ratio = c(0.6, 0.4),
           queries=list(list(query = Myfunc, params = list("Match"), color = brewer.pal(5,"Reds")[5], active = TRUE )))


tiff(width=4,height=3.75,res=300, units = "in",file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/Figures/F2D-upset.tiff") # or other device

AML
dev.off()



IDH1_IDH2<-test_set%>%filter(Group=="AML")%>%filter(IDH1==1&DNMT3A==0&ASXL1==0&TET2==0&IDH2==1)%>%pull(Sample)
TET2_IDH2<-test_set%>%filter(Group=="AML")%>%filter(IDH1==0&DNMT3A==0&ASXL1==0&TET2==1&IDH2==1)%>%pull(Sample)
TET2_IDH1<-test_set%>%filter(Group=="AML")%>%filter(IDH1==1&DNMT3A==0&ASXL1==0&TET2==1&IDH2==0)%>%pull(Sample)
DNMT3A_TET2_IDH2<-test_set%>%filter(Group=="AML")%>%filter(IDH1==0&DNMT3A==1&ASXL1==0&TET2==1&IDH2==1)%>%pull(Sample)
DNMT3A_TET2_IDH1<-test_set%>%filter(Group=="AML")%>%filter(IDH1==1&DNMT3A==1&ASXL1==0&TET2==1&IDH2==0)%>%pull(Sample)


# 38 AML samples with only 1 mutation and 42 with 2 or more mutations
test_set%>%filter(Group=="AML")%>%filter((ASXL1+DNMT3A+TET2+IDH1+IDH2)==1)%>%distinct(Sample)%>%summarize(Count=n())
test_set%>%filter(Group=="AML")%>%filter((ASXL1+DNMT3A+TET2+IDH1+IDH2)>=2)%>%distinct(Sample)%>%summarize(Count=n())

test_set%>%filter(Group=="AML")%>%filter((DNMT3A+TET2+IDH1+IDH2+ASXL1)>=2)%>%
                                          filter(DNMT3A==1&TET2==1&IDH1==0&IDH2==0)%>%
                                        distinct(Sample)%>%summarize(Count=n())/42*100

test_set%>%filter(Group=="AML")%>%filter((DNMT3A+TET2+IDH1+IDH2+ASXL1)>=2)%>%
  filter(DNMT3A==1&TET2==0&IDH1==0&IDH2==0)%>%
  distinct(Sample)%>%summarize(Count=n())/42*100

test_set%>%filter(Group=="AML")%>%filter((DNMT3A+TET2+IDH1+IDH2+ASXL1)>=2)%>%
  filter(DNMT3A==1&TET2==0&IDH1==0&IDH2==1)%>%
  distinct(Sample)%>%summarize(Count=n())/42*100

test_set%>%filter(Group=="AML")%>%filter((DNMT3A+TET2+IDH1+IDH2+ASXL1)>=2)%>%
  filter(Match=="Match")%>%
  distinct(Sample)%>%summarize(Count=n())/42*100



final_sample_summary[["CP5764d"]]$Clones
colnames(final_sample_summary[["CP5764d"]]$NGT)



TET2_c.3501.21_3501.20dupGT
ASXL1<-test_set%>%filter(Group=="AML")%>%filter(IDH1==0&DNMT3A==0&ASXL1==1&TET2==0&IDH2==0)%>%pull(Sample)

TET2_IDH1<-test_set%>%filter(Group=="AML")%>%filter(IDH1==1&DNMT3A==0&ASXL1==0&TET2==1&IDH2==0)%>%pull(Sample)
IDH1_IDH2<-test_set%>%filter(Group=="AML")%>%filter(IDH1==1&DNMT3A==0&ASXL1==0&TET2==0&IDH2==1)%>%pull(Sample)
TET2_IDH2<-test_set%>%filter(Group=="AML")%>%filter(IDH1==0&DNMT3A==0&ASXL1==0&TET2==1&IDH2==1)%>%pull(Sample)

DNMT3A_IDH1_subclone<-test_set%>%filter(Group=="AML")%>%filter(IDH1==1&DNMT3A==1&ASXL1==0&TET2==0&IDH2==0)%>%filter(Match=="Absent")%>%pull(Sample)
DNMT3A_IDH1_dominant<-test_set%>%filter(Group=="AML")%>%filter(IDH1==1&DNMT3A==1&ASXL1==0&TET2==0&IDH2==0)%>%filter(Match=="Match")%>%pull(Sample)

DNMT3A_IDH2_subclone<-test_set%>%filter(Group=="AML")%>%filter(IDH1==0&DNMT3A==1&ASXL1==0&TET2==0&IDH2==1)%>%filter(Match=="Absent")%>%pull(Sample)
DNMT3A_IDH2_dominant<-test_set%>%filter(Group=="AML")%>%filter(IDH1==0&DNMT3A==1&ASXL1==0&TET2==0&IDH2==1)%>%filter(Match=="Match")%>%pull(Sample)

DNMT3A_TET2_subclone<-test_set%>%filter(Group=="AML")%>%filter(IDH1==0&DNMT3A==1&ASXL1==0&TET2==1&IDH2==0)%>%filter(Match=="Absent")%>%pull(Sample)
DNMT3A_TET2_dominant<-test_set%>%filter(Group=="AML")%>%filter(IDH1==0&DNMT3A==1&ASXL1==0&TET2==0&IDH2==0)%>%filter(Match=="Match")%>%pull(Sample)


