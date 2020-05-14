source("/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Scripts/mBio_processing_scripts.R")


final_sample_summary<-readRDS(file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/final_sample_summary.rds")
clonal_sample_set_after_boostrap<-readRDS(file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/clonal_sample_set_after_boostrap")

mutants_in_each_sample<-do.call(rbind,lapply(final_sample_summary,function(x){
  y<-colnames(x$NGT)
  z <- list()
  z$TET2 <- ifelse(any(grepl("TET2",y)),1,0)
  z$DNMT3A <- ifelse(any(grepl("DNMT3A",y)),1,0)
  z$ASXL1 <- ifelse(any(grepl("ASXL1",y)),1,0)
  z$IDH1 <- ifelse(any(grepl("IDH1",y)),1,0)
  z$IDH2 <- ifelse(any(grepl("IDH2",y)),1,0)
  z$KRAS <- ifelse(any(grepl("KRAS",y)),1,0)
  z$NRAS <- ifelse(any(grepl("NRAS",y)),1,0)
  z$PTPN11 <- ifelse(any(grepl("PTPN11",y)),1,0)
  z$FLT3 <- ifelse(any(grepl("FLT3",y)),1,0)
  z$JAK2 <- ifelse(any(grepl("JAK2",y)),1,0)
  data.frame(t(do.call(rbind,z)))
}))


rownames(mutants_in_each_sample) <- names(final_sample_summary)
mutants_in_each_sample$Sample <- rownames(mutants_in_each_sample)
mutants_in_each_sample$Group <- ifelse(grepl("CH",mutants_in_each_sample$Sample),"CH","AML")
mutants_in_each_sample[mutants_in_each_sample$Sample%in%c("MA1715B","MA6300A","MA9521A","MA9521B"),"Group"]<- "MPN"
mutants_in_each_sample[mutants_in_each_sample$Sample%in%c("E4840","E4838new"),"Group"]<- "CMML"
mutants_in_each_sample[mutants_in_each_sample$Sample%in%c("MA4244A","MA2725","MA9521B","MA0092A","R2715","MA1715B","MA6363B"),"Group"]<- "MF"
mutants_in_each_sample$signal2 <- ifelse(rowSums(mutants_in_each_sample[,6:10])>=2,1,0)


patient_samples <- list(
  #"DNMT3A TET2"=mutants_in_each_sample%>%filter(Group=="AML")%>%filter(TET2==1&DNMT3A==1&IDH1==0&IDH2==0&ASXL1==0)%>%pull(Sample),
  "DNMT3A IDH1"=mutants_in_each_sample%>%filter(Group=="AML")%>%filter(TET2==0&DNMT3A==1&IDH1==1&IDH2==0&ASXL1==0)%>%pull(Sample),
  "DNMT3A IDH2"=mutants_in_each_sample%>%filter(Group=="AML")%>%filter(TET2==0&DNMT3A==1&IDH1==0&IDH2==1&ASXL1==0)%>%pull(Sample),
  "DNMT3A_samples"=mutants_in_each_sample%>%filter(Group=="AML")%>%filter(TET2==0&DNMT3A==1&IDH1==0&IDH2==0&ASXL1==0)%>%pull(Sample),
 # "TET2_samples"=mutants_in_each_sample%>%filter(Group=="AML")%>%filter(TET2==1&DNMT3A==0&IDH1==0&IDH2==0&ASXL1==0)%>%pull(Sample),
  "IDH1_samples"=mutants_in_each_sample%>%filter(Group=="AML")%>%filter(TET2==0&DNMT3A==0&IDH1==1&IDH2==0&ASXL1==0)%>%pull(Sample),
  "IDH2_samples"=mutants_in_each_sample%>%filter(Group=="AML")%>%filter(TET2==0&DNMT3A==0&IDH1==0&IDH2==1&ASXL1==0)%>%pull(Sample))


comutant_status<-do.call(rbind,lapply(names(patient_samples), function(gene_pair){
  x <-patient_samples[[gene_pair]]
  data.frame("Total"=mutants_in_each_sample%>%filter(Sample%in%x)%>%distinct(Sample)%>%summarise(Count=n())%>%pull(Count),
             "FLT3"=mutants_in_each_sample%>%filter(Sample%in%x&signal2==0)%>%tally(FLT3)%>%pull(n),
             "PTPN11"=mutants_in_each_sample%>%filter(Sample%in%x&signal2==0)%>%tally(PTPN11)%>%pull(n),
             "JAK2"=mutants_in_each_sample%>%filter(Sample%in%x&signal2==0)%>%tally(JAK2)%>%pull(n),
             "KRAS"=mutants_in_each_sample%>%filter(Sample%in%x&signal2==0)%>%tally(KRAS)%>%pull(n),
             "NRAS"=mutants_in_each_sample%>%filter(Sample%in%x&signal2==0)%>%tally(NRAS)%>%pull(n),
             "Multiple mutants"=mutants_in_each_sample%>%filter(Sample%in%x)%>%tally(signal2)%>%pull(n))
}))
  
rownames(comutant_status)<-  names(patient_samples)
  
comutant_status$None <-comutant_status[,"Total"]-rowSums(comutant_status[,-1])
comutant_status$Group <- rownames(comutant_status)
final<-comutant_status[,-1]%>%melt()
final$Group <- factor(do.call(rbind,strsplit(final$Group,split="_"))[,1],levels=rev(c("DNMT3A","TET2","ASXL1","IDH1","IDH2","DNMT3A TET2","DNMT3A IDH1","DNMT3A IDH2")))
final$variable <- gsub("\\."," ",final$variable)
final$variable <- factor(final$variable,levels=rev(c("JAK2","PTPN11","NRAS","KRAS","FLT3","Multiple mutants","None")))

color_set<-rev(brewer.pal(9,"Set1")[c(1,8,2,3,4,5,9)])

gg_co_mutants<-ggplot(final, aes(x=Group,fill=variable,y=value))+geom_bar(stat="identity",position="fill")+
                theme_classic(base_size=8)+
                ylab("Fraction of samples")+
                xlab("")+
                coord_flip()+
               # theme(axis.text.x = element_text(angle=30,hjust=1))+
                scale_y_continuous(expand=c(0,0))+
                scale_fill_manual(values=color_set,"Signaling mutation",guide = guide_legend(reverse = TRUE))


ggsave(addSmallLegend(gg_co_mutants),
       width=3.5,height=1.5,
       file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/Figures/F2G-comutant_distribution.pdf")

multi_mutant <- unlist(patient_samples[1:2])
patients_of_interest<-mutants_in_each_sample%>%filter(signal2==1&Sample%in%multi_mutant)%>%pull(Sample)


output_folder <- "/Users/bowmanr/Desktop/mBio_test/"  

mutually_exclusive<-setNames(lapply(patients_of_interest,function(y){
  x<-(final_sample_summary[[y]]$NGT)
  x<-x[,!grepl("Clone",colnames(x))]
  x[is.na(x)] <-0
  x[x>0] <-1
  x<- as.matrix(x)
  mode(x) <-"numeric"
  grob_corrplot_single_cell_loop<-generate_and_plot_cooccurence(x)
}),patients_of_interest)

lapply(patients_of_interest,function(sample){
data.frame(mutually_exclusive[[sample]]$data)[,10:12]#%>%filter(score!=0)
})


mutually_exclusive<-setNames(lapply(patients_of_interest,function(y){
  x<-(final_sample_summary[[y]]$NGT)
  x<-x[,!grepl("Clone",colnames(x))]
  x[is.na(x)] <-0
  x[x>0] <-1
  epi<-x[,grepl("DNMT3A",colnames(x))|grepl("TET2",colnames(x))|grepl("IDH1",colnames(x))|grepl("IDH2",colnames(x))]
  signal<-x[,grepl("PTPN11",colnames(x))|grepl("JAK2",colnames(x))|grepl("FLT3",colnames(x))|grepl("NRAS",colnames(x))|grepl("KRAS",colnames(x))]
  
  data.frame("Epigenetic"=sum(apply(epi,1,function(z){sum(z==1)>=2}))/nrow(epi),
             "Signaling"=sum(apply(signal,1,function(z){sum(z==1)>=2}))/nrow(signal))
         
}),patients_of_interest)

gg_fraction_comutated_cells<-ggplot(melt(do.call(rbind,mutually_exclusive)),aes(x=variable,y=value,fill=variable))+
                    geom_boxplot()+
                    geom_jitter(width=0.1)+
                    theme_classic(base_size=8)+
                    scale_fill_brewer(type="qual",palette = "Set1","Mutation pairs")+
                    xlab("")+ylab("Fraction of co-mutated cells")

melt(do.call(rbind,mutually_exclusive))%>%summarise(p=t.test(value~variable)$p.value)


ggsave((gg_fraction_comutated_cells),
       width=2.5,height=2,
       file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/Figures/EF3F-comutant_cell_fraction.pdf")

