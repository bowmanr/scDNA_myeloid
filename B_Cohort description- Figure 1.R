source("/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Scripts/mBio_processing_scripts.R")

final_NGTs<-readRDS(file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/final_NGTs.rds")
pheno_mut_melted<-readRDS(file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/pheno_mut_melted.rds")
final_sample_set<-readRDS(file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/final_sample_set.rds")


final_mut_melt <-pheno_mut_melted%>%filter(Sample.ID%in%final_sample_set)
final_mut_melt$Gene<- factor(final_mut_melt$Gene,levels=names(sort(table(final_mut_melt$Gene), decreasing=TRUE)))

## Total number of mutations identified
gg_mut_count<-ggplot(final_mut_melt,aes(x=Gene))+geom_bar(stat="count")+
  theme_classic(base_size = 8)+ylab("Count")+ggtitle("Number of mutations")+
  theme(axis.text.x = element_text(angle=45, hjust=1,vjust=1))

## Melt to get a tally of how many mutations per patient
melted_mut_mat <- data.frame(count(final_mut_melt, Gene, Sample.ID))
melted_mut_mat$Gene<- factor(melted_mut_mat$Gene,levels=names(sort(table(melted_mut_mat$Gene),decreasing=TRUE)))
gg_mut_patient<-ggplot(melted_mut_mat,aes(x=Gene))+geom_bar(stat="count")+theme_classic(base_size = 8)+ylab("Count")+ggtitle("Number of patients with mutation")+theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1))

## Total nunber of mutant genes per patient
gg_mutated_genes_per_patient<-melted_mut_mat%>%group_by(Sample.ID)%>%tally%>%ggplot(aes(x=n))+geom_bar()+theme_classic(base_size = 8)+ylab("Count")+ggtitle("Mutant genes per patient")+xlab("Number of genes")

## Total number of mutations per patient
melted_variant_mat <- data.frame(count(final_mut_melt, protein, Sample.ID))
melted_variant_mat$protein<- factor(melted_variant_mat$protein, levels=names(sort(table(melted_variant_mat$protein), decreasing=TRUE)))
gg_mutations_per_patient<-melted_variant_mat%>%group_by(Sample.ID)%>%tally%>%ggplot(aes(x=n))+geom_bar()+theme_classic(base_size = 10)+ylab("Count")+ggtitle("Variants per patient")+xlab("Number of variants")

### Cells per sample
gg_cells_per_patient<-data.frame("Cells"=do.call(rbind,lapply(final_NGTs,nrow)))%>%ggplot(aes(x=Cells))+geom_histogram(binwidth = 100)+theme_classic(base_size = 10)+ylab("Count")+ggtitle("Informative Cells per Sample") #+scale_x_continuous(breaks=c(0,500,1000,2000,4000,8000))

### create matrix for oncoprint
mut_mat <- table(melted_mut_mat$Sample.ID,melted_mut_mat$Gene)
mut_mat <- mut_mat[,colSums(mut_mat)!=0]
tronco_mat <-import.genotypes(mut_mat, event.type='myVariant')
grob_oncoprint<-grid.grabExpr(oncoprint(tronco_mat))

### Prepare matrix for co occurence map
grob_corrplot<-generate_and_plot_cooccurence(mut_mat)

#EF1A
ggsave(gg_mut_count+ theme(plot.title = element_text(hjust = 0.5)), width=3.75,height=2.25,
        file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/Figures/EF1A-mutation_count.pdf")

#EF1B
ggsave(gg_mut_patient+ theme(plot.title = element_text(hjust = 0.5)), width=3.75,height=2.25,
       file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/Figures/EF1B-number_of_patients_with_mutation.pdf")

#EF1C
ggsave(gg_mutated_genes_per_patient+ theme(plot.title = element_text(hjust = 0.5)), width=3.75,height=2.25,
       file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/Figures/EF1C-mutant_genes_per_patient.pdf")

#EF1D
ggsave(gg_mutations_per_patient+ theme(plot.title = element_text(hjust = 0.5)), width=3.75,height=2.25,
       file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/Figures/EF1D-variants_per_patient.pdf")




final_mut_melt$Variant_Class<-ifelse(grepl("fs\\*",final_mut_melt$protein),"Indel",
                                 ifelse(grepl("INS_",final_mut_melt$protein),"Indel",
                                        ifelse(grepl("ins",final_mut_melt$protein),"Indel",
                                               ifelse(grepl("ext",final_mut_melt$protein),"Indel",
                                                      ifelse(grepl("del",final_mut_melt$protein),"Indel",
                                                           ifelse(grepl("\\*$",final_mut_melt$protein),"Nonsense","Missense"))))))
final_mut_melt$Gene <- as.character(final_mut_melt$Gene)
final_mut_melt$Variant_Class <- as.character(final_mut_melt$Variant_Class)
final_mut_melt$Sample.ID <- as.character(final_mut_melt$Sample.ID)


mut_mat_wide<-data.frame(final_mut_melt%>%group_by(Sample.ID)%>%pivot_wider(id_cols=Sample.ID, names_from = Gene, values_from = Variant_Class,values_fn = list(Variant_Class = list)),stringsAsFactors = FALSE)
mut_mat_wide[is.na(mut_mat_wide)] = " "

mut_mat_wide_storage <- list()
for(i in 2:ncol(mut_mat_wide)){
  mut_mat_wide_storage[[i]]<- do.call(rbind,lapply(mut_mat_wide[,i],function(x){paste(x,sep=";",collapse=";")}))
}
mut_mat_wide_storage[[1]] <- mut_mat_wide[,1]
final_mat<-do.call(cbind,mut_mat_wide_storage)
colnames(final_mat) <- colnames(mut_mat_wide)
rownames(final_mat) <- final_mat[,1]
final_mat <- t(final_mat[,-1])



col = c("Indel" = "darkorchid2", "Nonsense" = "black", "Missense" = "darkgreen")
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w, h-unit(0.25, "mm"), 
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  # big blue
  Indel = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Indel"], col = NA))
  },
  # bug red
  Nonsense = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Nonsense"], col = NA))
  },
  # small green
  Missense = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, 
              gp = gpar(fill = col["Missense"], col = NA))
  }
)


library(ComplexHeatmap)

sample_key<-read.delim(file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/From_LAM/OS_and_Cheatsheet/Mission_Bio_MRN_RunID_Samples.txt",sep="\t")
sample_clinical<-read.delim(file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/From_LAM/OS_and_Cheatsheet/sample_clinical_information.txt",sep="\t")

merged_set<-full_join(sample_key,sample_clinical,by="Molecular_accession_num")

missing_patient_info<-setdiff(colnames(final_mat),merged_set%>%filter(Sample.Name%in%colnames(final_mat))%>%distinct(Sample.Name)%>%pull(Sample.Name))
complete_patient_info<-intersect(colnames(final_mat),merged_set%>%filter(Sample.Name%in%colnames(final_mat))%>%distinct(Sample.Name)%>%pull(Sample.Name))


set<-rbind(data.frame(merged_set%>%filter(Sample.Name%in%complete_patient_info)%>%select(Dx,Disease.Status,Sample.Name)),
      data.frame("Dx"=ifelse(grepl("^CH",missing_patient_info),"CH","AML"),
                 "Disease.Status"=NA,
                "Sample.Name"=missing_patient_info))

set[set$Sample.Name%in%c("S5522"),"Dx"] <- "tAML"
set[set$Sample.Name%in%c("F7469","F7469dnew","MA4244B","MA6300B","MA6363C"),"Dx"] <- "sAML"
set[set$Sample.Name%in%c("RF2768","R6501","CB0346","CB4833","E3083new","E6502","E6522","E8261","F5348","R2278new"),"Dx"] <- "AML"

set[set$Sample.Name%in%c("S5046","S5331","S4410","S6225","S5522","CB4833","E3083new","E6502"),"Disease.Status"] <- "Relapse"
set[set$Sample.Name%in%c("S3326","S5346","S5480","S5941","Meyer_1","17-005","17-006","17-020","17-021"),"Disease.Status"] <- "Newly Diagnosed"
set[set$Sample.Name%in%c("RF2768","R6501","CB0346","F5348","R2278new"),"Disease.Status"] <- "Newly Diagnosed"
set[set$Sample.Name%in%c("E6522","F7469","F7469dnew","MA6300B","MA6363C"),"Disease.Status"] <- "Newly Transformed"
set[set$Sample.Name%in%c("E8261","MA4244B"),"Disease.Status"] <- "Refractory"


set[set$Sample.Name%in%c("MA0092A","MA1715B","MA4244A","MA4849B","MA6300A","MA6363B","MA9521A","MA9521B"),"Dx"] <- "MPN"
set[set$Sample.Name%in%c("MA0092A","MA1715B","MA4244A","MA4849B","MA6300A","MA6363B","MA9521A","MA9521B"),"Disease.Status"] <- "Chronic Stage"

set_dedup<-set[!duplicated(as.matrix(set%>%filter(Sample.Name%in%colnames(final_mat)))),]

set_dedup<-read.delim(file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/From_LAM/OS_and_Cheatsheet/BB_merged_sample_clinical_information.txt",sep="\t",row.names=4)[,-1]

colnames(set_dedup) <-c("Diagnosis","Status")
set_dedup<- set_dedup[colnames(final_mat),]
set_dedup[is.na(set_dedup[,2]),"Status"] <- "Other"

color_set <- list("Diagnosis" = setNames(okabe(n=length(unique(set_dedup[,1]))),c( "CH","MDS","MPN", "AML", "sAML" ,"tAML","Other") ),
                       "Status" = setNames(tol(n=(length(unique(set_dedup[,2])))),c("Newly Diagnosed","Relapse/Refractory","Newly Transformed","Persistent","Chronic Stage", "Other")))                       
                                                    
top_annotation=HeatmapAnnotation(cbar = anno_oncoprint_barplot(),
                  df=set_dedup[,1:2],
                  col =list("Diagnosis"=color_set[[1]],
                            "Status"=color_set[[2]])
                  )



column_title = ""
heatmap_legend_param = list(title = "Alternations", 
                            at = c("Indel", "Nonsense", "Missense"), 
                            labels = c("Indel", "Nonsense", "Missense"))
oncoPrint(final_mat,
          alter_fun = alter_fun, 
          col = col, 
          top_annotation=top_annotation,
         #heatmap_legend_param = list(direction = "horizontal"),
          #column_title = column_title# 
          heatmap_legend_param = heatmap_legend_param)
