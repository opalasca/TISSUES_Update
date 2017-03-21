plot_multi_fold_enrichment<-function(data,out_file,title,labels){
  png(out_file,height=1000,width=900)
  g <- gridExtra::borderGrob(type=9, colour="black", lwd=2)
  plot<- ggplot(data, aes(mean_score, fold_enrichment,color=set)) + 
    colScale+
    facet_grid(set~.,space="free",scales="free_x") +
    #ggtitle(title) +
    xlab(label="Mean score")+
    ylab(label="Fold enrichment compared with UniProt")+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())+
    annotation_custom(g)+
    scale_x_log10(breaks=c(0.1,1,10, 100,1000, 5000))+
    annotation_logticks(sides="b", short = unit(0.02, "cm"), 
                        mid = unit(0.05, "cm"), long = unit(0.1, "cm"),color="black")+
    theme(strip.background=element_rect(fill='white',color='white'))+
    theme(strip.text.y = element_blank())+
    theme(axis.text.x = element_text(family= "Helvetica",colour = "black"))+
    theme(axis.text.y = element_text(family= "Helvetica",colour = "black"))+
    theme(axis.text = element_text(face="bold",family= "Helvetica",size=rel(1.2)))+
    theme(axis.title = element_text(family= "Helvetica",size=rel(2)))+
    theme(legend.title = element_blank())+
    theme(legend.position="bottom")+
    theme(legend.text = element_text(family= "Helvetica",size = rel(1.6)))+
    theme(legend.key.width=unit(1,"line"))+
    theme(legend.key = element_rect(fill = "white",colour = "white")) + 
    guides(colour = guide_legend(override.aes = list(size=4)))+
    theme(legend.background = element_blank())
 
  for(j in 1:length(labels)){
    dataset<-labels[j]
    susbset_data<-subset(data,set = dataset)
    #method="gam",formula=y~s(x,bs='tp')
    plot<- plot + geom_point(data=susbset_data,aes(x=mean_score,y=fold_enrichment), size=3)
       #scale_size_area()+	
  }
  print(plot)
  dev.off()
}
  

get_mRNA_data<-function(levels,labels){
  dataset<- levels[1]
  file<-paste(dataset,"uniprot_fold_enrichment_analysis.tsv",sep="_")
  
  all_data<-read.csv(file=file,header=FALSE,sep="\t")
  names(all_data)<-c("mean_stars","fold_enrichment","mean_score")
  all_data$set<-dataset
  
  for (i in 2:length(levels) ) {
    dataset<- levels[i]
    file<-paste(dataset,"uniprot_fold_enrichment_analysis.tsv",sep="_")
    data<-read.csv(file=file,header=FALSE,sep="\t")
    names(data)<-c("mean_stars","fold_enrichment","mean_score")
    data$set<-dataset
    all_data<-rbind(all_data,data)
  }
  all_data$set<-factor(all_data$set,levels=levels,labels=labels)
  return(all_data)
}

###############################
#           Analysis         #
##############################
#Colors
col<-get_colors(is_label = T)
colScale <- scale_color_manual(name = "",values = col,labels=c("Human Exon Array","Human GNF","Human RNA-seq atlas","Human HPA RNA-seq","Mouse GNF", "Mouse GNF V3", "Mouse RNA-seq Encode", "Mouse RNA-seq MIT", "Rat Array", "Rat RNA-seq MIT", "Rat RNA-seq BodyMap", "Pig Array", "Pig RNA-seq Aarhus", "Pig RNA-seq WUR"))

#Get mRNA data
levels<-c("exon","gnf","rna","hpa_rna","mouse_gnf","mouse_gnfv3","mouse_rnaseq_encode","mouse_rnaseq_mit", "rat_array", "rat_rnaseq_mit", "rat_rnaseq_bodymap", "pig_array", "pig_rnaseq_aarhus", "pig_rnaseq_wur")
labels=c("Human Exon Array","Human GNF","Human RNA-seq atlas","Human HPA RNA-seq","Mouse GNF", "Mouse GNF V3", "Mouse RNA-seq Encode", "Mouse RNA-seq MIT", "Rat Array", "Rat RNA-seq MIT", "Rat RNA-seq BodyMap", "Pig Array", "Pig RNA-seq Aarhus", "Pig RNA-seq WUR")

data<-get_mRNA_data(levels,labels)  

#Plot fold-enrichment
title <- "Datasets fold enrichment per score"
out_file <- "../figures/datasets_fold_enrichment.png"
plot_multi_fold_enrichment(data,out_file,title,labels)

