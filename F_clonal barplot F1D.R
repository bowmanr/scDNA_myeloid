source("/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Scripts/mBio_processing_scripts.R")

final_sample_summary<-readRDS(file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/final_sample_summary.rds")
clonal_sample_set_after_boostrap<-readRDS(file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/clonal_sample_set_after_boostrap")

sample <-"MSK45"

  clonal_abundance_subset <-final_sample_summary[[sample]]$Clones 
  clonal_architecture <-final_sample_summary[[sample]]$Architecture 

    clones_to_use <- intersect(as.character(clonal_abundance_subset$Clone),as.character(clonal_architecture$Clone))
    clonal_architecture_subset <- data.frame(clonal_architecture)%>%
                                             filter(data.frame(clonal_architecture)$Clone%in%clones_to_use)
    clonal_abundance_subset <- data.frame(clonal_abundance_subset)%>%
                                            filter(Clone%in%clones_to_use)
    genes_to_display <-unique(as.character(clonal_architecture_subset[clonal_architecture_subset$Genotype%in%c("Heterozygous","Homozygous"),"Mutant"]))
    
    clonal_architecture_subset$Clone <- factor(clonal_architecture_subset$Clone, levels=rev(clonal_abundance_subset$Clone))
    clonal_abundance_subset$Clone <- factor(clonal_abundance_subset$Clone, levels=levels(clonal_architecture_subset$Clone))
    
    
    gg_heatmap <- ggplot(data = clonal_architecture_subset, aes(x = Clone, y = factor(Mutant,levels=rev(levels(factor(Mutant)))), fill = Genotype)) + 
      geom_tile() +# scale_y_discrete(limits = rev(levels(Mutant)))+
      
      scale_fill_manual(values=c("WT"=brewer.pal(7,"Reds")[1],
                                 "Heterozygous"=brewer.pal(7,"Reds")[3],
                                 "Homozygous"=brewer.pal(7,"Reds")[6],
                                 "Unknown"="grey50"),"Genotype")  +
      theme_classic(base_size=7) +
      ylab("Mutation")+
      theme(legend.position = "none", legend.direction = "vertical",
            axis.text.x = element_blank(), 
            axis.line=element_blank(),
            axis.title.x=element_blank(),
            axis.ticks.x = element_blank(),
            plot.margin=unit(c(0,0,0,0),"cm"))
    
    # Generate clonal abundance barplot
    gg_clonal_barplot <- ggplot(data = data.frame(clonal_abundance_subset), aes(x = Clone, y = Count)) + 
      geom_bar(stat = "identity", aes(fill = Count)) + theme_gray() +
      theme_classic(base_size=7)+
      scale_y_continuous(expand=c(0.01,0))+
      #ylim() + 
      ylab("Cell Count")+
      geom_errorbar(aes(ymin = LCI, ymax = UCI), width = 0.2)+
    #  geom_text(aes(label=Count), position=position_dodge(width=0.9), vjust=-0.25,size=1)+
      scale_fill_distiller(name = "Value", palette = "Reds", direction = 1) +
      theme(axis.title.x = element_blank(), 
            axis.text.x = element_blank(), 
            axis.ticks.x = element_blank(),
            axis.line.x =element_blank(),
            legend.position = "none",
            plot.margin=unit(c(0,0,0,0),"cm"))
    
    grob.title <- textGrob(paste("",""), hjust= 0.5, vjust = 0.5, gp = gpar(fontsize = 12))
    final_plot<-plot_grid(grob.title,gg_clonal_barplot,gg_heatmap,ncol=1,align="v",axis="l",rel_heights = c(0,1,0.75))
    
    ggsave(final_plot, width=4.5,height=3.5,
           file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/Figures/F1D-clonal_barplot.pdf")
    
    