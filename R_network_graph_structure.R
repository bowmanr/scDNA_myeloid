source("/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Scripts/mBio_processing_scripts.R")

final_sample_summary<-readRDS(file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/final_sample_summary.rds")
clonal_sample_set_after_boostrap<-readRDS(file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/clonal_sample_set_after_boostrap")

graph_results<-readRDS(file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/MDP_results.rds")
final_results<-readRDS(file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/MDP_trajectory_for_each_gene.rds")


plyr::compact(lapply(final_results,function(x){names(x)[c(grep("DNMT3A",names(x)),grep("TET",names(x)) )]}))


DNMT3A IDH NPM1 final sample 
E2145

NPM1 IDH FLT3 RAS
R5853


quartz()
barplot_for_observed_states("E2145")
quartz()
plot_optimal_graph_for_trajectory("E2145")


final_sample_summary[["S5046"]]$Clones

colnames(final_sample_summary[["R5853"]]$NGT)



plot_optimal_graph_for_trajectory<-function(sample){
        all_results <-final_results[[sample]]

        all_results_filtered<-lapply(all_results,function(storage_results){
        storage_results[lapply(storage_results,function(x){sum(x$reward,na.rm = TRUE)})==0]<-NULL
        storage_results<-lapply(storage_results,function(x){x[1:which.max(x$cumulative_reward),]})
        final<-do.call(rbind,storage_results)[,c(1,4,3)]
        nodes_to_remove <- setdiff(setdiff(final$nextState,final$states),final_sample_summary[[sample]]$Clones$Clone)
        final <- final%>%filter(!nextState%in%nodes_to_remove)
        
        nodes_to_remove <- setdiff(setdiff(final$nextState,final$states),final_sample_summary[[sample]]$Clones$Clone)
        final <- final%>%filter(!nextState%in%nodes_to_remove)
        final[!duplicated(final),]
      })

      all_results_filtered<-setNames(lapply(names(all_results_filtered),function(x){
        data.frame(all_results_filtered[[x]],
                   "mutation"=x,
                   "edge"=paste(all_results_filtered[[x]]$states,all_results_filtered[[x]]$nextState,sep="->"))
      }),names(all_results))

      optimal<-names(which.max(do.call(c,lapply(all_results_filtered,function(x){
        sum(x$reward)
        }))))
      print(optimal)
      final<-do.call(rbind,all_results_filtered)%>%filter(mutation==optimal)
      final%>%group_by(mutation)%>%distinct(nextState,reward)%>%summarize(total=sum(reward))

      graph<-graph_from_data_frame(final,directed=T)
      weights<-final_sample_summary[[sample]]$Clones$Count/sum(final_sample_summary[[sample]]$Clones$Count)
      names(weights) <-final_sample_summary[[sample]]$Clones$Clone
      weight_subset<-weights[names(weights)%in%names(V(graph))]
      nodes_to_add_names<-setdiff(names(V(graph)),names(weights))
      nodes_to_add <- rep(0.1,length(nodes_to_add_names))
      names(nodes_to_add)<-nodes_to_add_names
      weight_final <- c(weight_subset,nodes_to_add)[names(V(graph))]

      mutations<-names(all_results_filtered)


      clone_colors <-ifelse(names(V(graph))%in%final_sample_summary[[sample]]$Clones$Clone,brewer.pal(5,"Reds")[5],"grey80")

      #par(mar=c(0,0,0,0))
      #pdf(width=8,height=8,file="/Volumes/LevineLab/Levine Lab/Bobby/Collaborations/MissionBio/Analysis/2020/January/Figures/F2E-network_DTAI_test.pdf") # or other device
  #    quartz()
     (plot(graph,layout=layout_as_tree,
                  vertex.color=ifelse(names(V(graph))%in%final_sample_summary[[sample]]$Clones$Clone,brewer.pal(5,"Reds")[5],"grey80"),
                  vertex.frame.color=ifelse(names(V(graph))%in%final_sample_summary[[sample]]$Clones$Clone ,brewer.pal(5,"Reds")[5],"grey80"),
                  vertex.size=weight_final*100,#main=c(optimal),annotate.plot=TRUE,
                  vertex.label=NA,
                   edge.color="grey30",#ifelse(edge_attr(graph)$mutation%in%mutations[1],brewer.pal(5,"Oranges")[4],
                               #   ifelse(edge_attr(graph)$mutation%in%mutations[2],brewer.pal(5,"Blues")[4],brewer.pal(5,"Greens")[4])), 
                  arrow.size=0.01,arrow.width=0.01))
      
       }


barplot_for_observed_states <- function(sample){
        all_results <-final_results[[sample]]

        all_results_filtered<-lapply(all_results,function(storage_results){
        storage_results[lapply(storage_results,function(x){sum(x$reward,na.rm = TRUE)})==0]<-NULL
        storage_results<-lapply(storage_results,function(x){x[1:which.max(x$cumulative_reward),]})
        final<-do.call(rbind,storage_results)[,c(1,4,3)]
        nodes_to_remove <- setdiff(setdiff(final$nextState,final$states),final_sample_summary[[sample]]$Clones$Clone)
        final <- final%>%filter(!nextState%in%nodes_to_remove)
        
        nodes_to_remove <- setdiff(setdiff(final$nextState,final$states),final_sample_summary[[sample]]$Clones$Clone)
        final <- final%>%filter(!nextState%in%nodes_to_remove)
        final[!duplicated(final),]
      })

      all_results_filtered<-setNames(lapply(names(all_results_filtered),function(x){
        data.frame(all_results_filtered[[x]],
                   "mutation"=x,
                   "edge"=paste(all_results_filtered[[x]]$states,all_results_filtered[[x]]$nextState,sep="->"))
      }),names(all_results))


      final<-do.call(rbind,all_results_filtered)


      final$observed <- ifelse(final$nextState%in%final_sample_summary[[sample]]$Clones$Clone,"Observed","Unobserved")


      out<-data.frame("Mutation"=final%>%distinct(mutation)%>%pull(mutation),
                            "Reward"=final%>%group_by(mutation)%>%summarize(max_reward=sum(reward))%>%pull(max_reward),
                            "Observed"=final%>%group_by(mutation)%>%filter(observed=="Observed")%>%distinct(nextState)%>%summarize(Observed=n())%>%pull(Observed),
                            "Clones"=final%>%group_by(mutation)%>%distinct(nextState)%>%summarize(Observed=n())%>%pull(Observed))
      out$Unobserved <- out$Clones-out$Observed
      out_melt <- melt(out[,-4],id.vars=c("Mutation","Reward"))
      out_melt$variable <- factor(out_melt$variable,levels=c("Unobserved","Observed"))

      out_melt$Mutation <-factor(out_melt$Mutation,levels=out$Mutation[order(out$Unobserved/out$Observed)])


      ggA<-ggplot(out_melt,aes(x=Mutation,y=value,fill=variable))+
                  geom_bar(stat="identity")+
                  scale_fill_manual(values=c("Observed"=brewer.pal(5,"Reds")[5],
                                              "Unobserved"="grey80"))+theme_classic(base_size=8)+
                  scale_y_continuous(expand=c(0,0))+
                  theme(axis.text.x=element_text(angle=30,hjust=1))

     ggB<-ggplot(out_melt%>%filter(variable=="Observed"),aes(x=Mutation,y=Reward))+
                  geom_bar(stat="identity")+
                  scale_fill_manual(values=c("Observed"=brewer.pal(5,"Reds")[5],
                                              "Unobserved"="grey80"))+theme_classic(base_size=8)+
                  scale_y_continuous(expand=c(0,0))+
                  theme(axis.text.x=element_text(angle=30,hjust=1))


    return(plot_grid(ggA,ggB,ncol=2))
}
         