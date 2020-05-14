source("/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Scripts/mBio_processing_scripts.R")

final_sample_summary<-readRDS(file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/final_sample_summary.rds")
clonal_sample_set_after_boostrap<-readRDS(file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/clonal_sample_set_after_boostrap")

graph_results<-readRDS(file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/MDP_results.rds")
final_results<-readRDS(file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/MDP_trajectory_for_each_gene.rds")

A<- "NPM1"
B<-"FLT3"

double_mutant_samples<-names(graph_results)[do.call(c,lapply(names(graph_results),function(x){
  sample <- graph_results[[x]]
  any(grepl(A,unique(sample$Action)))&any(grepl(B,unique(sample$Action)))
}))]



double_mutant_order<-setNames(lapply(double_mutant_samples,function(x){
  print(x)
  sample <- final_sample_summary[[x]]
  mutations <- colnames(sample$NGT)
  clones<- sample$Clone
  A_het_clones <-sample$Architecture%>%filter(Genotype=="Heterozygous")%>%filter(grepl(A,Mutant))%>%distinct(Clone)%>%pull(Clone)
  B_het_clones <-sample$Architecture%>%filter(Genotype=="Heterozygous")%>%filter(grepl(B,Mutant))%>%distinct(Clone)%>%pull(Clone)
  
  AB_het_clones<-intersect(A_het_clones,B_het_clones)
  AB_clone_sizes <- graph_results[[x]]%>%filter(NextState%in%AB_het_clones)%>%filter(grepl(A,Action)|grepl(B,Action))%>%filter(Reward==max(Reward))%>%distinct(Reward)%>%pull(Reward)
  antecdenent_options<-graph_results[[x]]%>%filter(NextState%in%AB_het_clones)%>%filter(grepl(A,Action)|grepl(B,Action))%>%filter(Reward==max(Reward))
  
  antecendent_clones<-graph_results[[x]]%>%filter(NextState%in%antecdenent_options$State)%>%filter(Reward!=0)
  
  set<-setNames(merge(antecdenent_options[,1:4],antecendent_clones[,3:4],by.x="State",by.y="NextState"),
                c("Antecedent","Mutation","Max_het_state","Max_het_size","Antecednet_size"))
  
  Mut_A <- set%>%filter(grepl(A,set$Mutation))%>%distinct(Mutation,Antecednet_size)%>%filter(Antecednet_size==max(Antecednet_size))%>%pull(Mutation)
  Mut_B <- set%>%filter(grepl(B,set$Mutation))%>%distinct(Mutation,Antecednet_size)%>%filter(Antecednet_size==max(Antecednet_size))%>%pull(Mutation)
  
  set <- set%>%filter(Mutation%in%c(Mut_A,Mut_B))
  
  convergence <- ifelse((any(grepl(A,set$Mutation)) &any(grepl(B,set$Mutation))),"Yes","No")
  
  
  print(convergence)
  
  if(convergence=="Yes"){
    sizes <- setNames(data.frame("Double_mutant"=unique(set$Max_het_size),
                                 B=unique(set%>%filter(grepl(A,set$Mutation))%>%pull(Antecednet_size)),
                                 A=unique(set%>%filter(grepl(B,set$Mutation))%>%pull(Antecednet_size))),
                      c("Double_mutant",B,A))
  } else if(any(grepl(B,set$Mutation))) {
    sizes <- setNames(data.frame("Double_mutant"=unique(set$Max_het_size),
                                 A=unique(set%>%filter(grepl(B,set$Mutation))%>%pull(Antecednet_size))),
                      c("Double_mutant",A))
  } else if(any(grepl(A,set$Mutation))) {
    sizes <- setNames(data.frame("Double_mutant"=unique(set$Max_het_size),
                                 B=unique(set%>%filter(grepl(A,set$Mutation))%>%pull(Antecednet_size))),
                      c("Double_mutant",B))
  } else if(nrow(set)<1){
    return(data.frame("Double_mutant"=0))
  }
  return(sizes)
}),double_mutant_samples)


double_mutant_order_mat<-plyr::rbind.fill(setNames(lapply(double_mutant_order,function(x){data.frame(x)}),c(double_mutant_samples)))
rownames(double_mutant_order_mat) <-    double_mutant_samples
double_mutant_order_mat$Sample <- double_mutant_samples

double_mutant_order_mat<- double_mutant_order_mat%>%filter(!Double_mutant==0)

double_mutant_order_mat$Convergence <- ifelse((!is.na(double_mutant_order_mat[,A])&!is.na(double_mutant_order_mat[,B]))&  
                                                double_mutant_order_mat[,A]>0&
                                                double_mutant_order_mat[,B]>0,"Convergent",
                                              ifelse(is.na(double_mutant_order_mat[,A])&double_mutant_order_mat[,B]>0,B,A))
melted_set <-setNames(melt(double_mutant_order_mat),c("Sample","Convergence","Clone","Size"))

melted_set$Clone <-factor(melted_set$Clone,levels=c(A,"Double_mutant",B))                                     
melted_set$Clone <-plyr::revalue(melted_set$Clone,c("Double_mutant"="Double mutant"))                                     
melted_set$Convergence <-factor(melted_set$Convergence,levels=c("Convergent",A,B))  
gene_name_A <- melted_set%>%filter(Convergence==A)

colors_set <- c(brewer.pal(5,"Reds")[5],
            brewer.pal(5,"Blues")[5],
            brewer.pal(5,"Greens")[5])
names(colors_set) <- c("Convergent",A,B)
alpha_set <- c(0.5,1,1)
names(alpha_set) <- names(colors_set)
ggA<-ggplot(melted_set,aes(x=Clone,y=Size/100,group=Sample,color=Convergence))+geom_point()+
  geom_line()+ggtitle(paste(A,B))+
  theme_bw(base_size=6)+theme(plot.title = element_text(hjust=0.5))+
  ylab("Fraction of sample in mutant clone")+
  xlab("")+
  scale_alpha_manual(values=alpha_set,guide=FALSE)+
  scale_color_manual(values=colors_set,drop=FALSE)

ggsave(ggA, width=2.75,height=1.7,
       file=paste0("/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/Figures/F3F1-convergence ",A,B,".pdf"))


melted_set%>%filter(Clone%in%c(A,B))%>%summarize(p=t.test(Size~Clone)$p.value)
melted_set%>%filter(Clone%in%c(A,"Double_mutant"))%>%summarize(p=t.test(Size~Clone)$p.value)
melted_set%>%filter(Clone%in%c(B,"Double_mutant"))%>%summarize(p=t.test(Size~Clone)$p.value)


ggB<-ggplot(melted_set%>%distinct(Sample,Convergence),
            aes(x=Convergence,fill=Convergence))+
  geom_bar()+theme_classic()+
  scale_y_continuous(expand=c(0,0))


ggC<-ggplot(melted_set%>%filter(Clone%in%c(A,B)),aes(x=Clone,y=Size/100,fill=Clone))+
  geom_boxplot()+
  theme_bw(base_size=10)+
  ylab("Fraction of sample in mutant clone")+
  xlab("")#+

plot_grid(ggA,ggB,ggC,align="hv",axis="ltrb",ncol=3)
