
```{r,eval=TRUE}

patients_of_interest <- intersect(multi_DTAI,multi_signaling)

mutually_exclusive<-setNames(lapply(patients_of_interest,function(y){
  print(y)
  x<-(final_sample_summary[[y]]$NGT)
  x<-x[,!grepl("Clone",colnames(x))]
  
  if(!any(grepl("DNMT3A",colnames(x)))){
    return(NULL)
  } else{
    x[is.na(x)] <-0
    x[x>0] <-1
    epi<-x[,grepl("DNMT3A",colnames(x))|grepl("IDH1",colnames(x))|grepl("IDH2",colnames(x))]
    signal<-x[,grepl("PTPN11",colnames(x))|grepl("JAK2",colnames(x))|grepl("FLT3",colnames(x))|grepl("NRAS",colnames(x))|grepl("KRAS",colnames(x))]
    
    if(class(epi)=="numeric"){
      return(NULL)
    } else if(ncol(epi)==0){
      return(NULL)
    } else {  
      data.frame("Sample"=y,
                 "Epigenetic"=sum(apply(epi,1,function(z){sum(z==1)>=2}))/nrow(epi),
                 "Signaling"=sum(apply(signal,1,function(z){sum(z==1)>=2}))/nrow(signal))
    } }
  
}),patients_of_interest)

gg_fraction_comutated_cells<-ggplot(melt(do.call(rbind,mutually_exclusive)),aes(x=variable,y=value,fill=variable))+
  geom_boxplot()+
  #   geom_path()+
  geom_jitter(width=0.1)+
  theme_classic(base_size=8)+
  scale_fill_brewer(type="qual",palette = "Set1","Mutation pairs")+
  xlab("")+ylab("Fraction of co-mutated cells")

gg_fraction_comutated_cells
```