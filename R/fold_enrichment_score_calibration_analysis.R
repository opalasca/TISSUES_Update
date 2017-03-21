
datasets<-c("tm","mouse_gnf","mouse_gnfv3","mouse_rnaseq_encode","mouse_rnaseq_mit", "rat_array", "rat_rnaseq_mit", "rat_rnaseq_bodymap", "pig_array", "pig_rnaseq_aarhus", "pig_rnaseq_wur")
dataset_labels<- c("Text-mining", "Mouse GNF","Mouse GNF V3","Mouse RNA-seq Encode", "Mouse RNA-seq MIT", "Rat Array", "Rat RNA-seq BodyMap", "Rat RNA-seq MIT" ,"Pig Array", "Pig RNA-seq Aarhus", "Pig RNA-seq WUR") 

# datasets<-c("tm","mouse_rnaseq_encode","mouse_rnaseq_mit", "rat_array")
# dataset_labels<- c("Text-mining", "Mouse RNA-seq Encode", "Mouse RNA-seq MIT", "Rat Array") 

title <- "Datasets fold enrichment per confidence score"

out_file <- "../figures/datasets_score_calibration.png"
png(out_file,height=750,width=750)#,useDingbats=FALSE)

col<-get_colors()
colScale <- scale_fill_manual(name = "",values = col,labels=c("Text-mining","Mouse GNF", "Mouse GNF V3", "Mouse RNA-seq Encode", "Mouse RNA-seq MIT", "Rat Array", "Rat RNA-seq MIT", "Rat RNA-seq BodyMap", "Pig Array", "Pig RNA-seq Aarhus", "Pig RNA-seq WUR"))

plot<- ggplot()+colScale
y <- 1.5 
j <- 0.2
for (i in 1:length(datasets)){
  dataset<- datasets[i]
  label<- dataset_labels[i]
  file<-paste(dataset,"uniprot_fold_enrichment_analysis.tsv",sep="_")
  color<-col[names(col)==dataset][1]
  data<-read.csv(file=file,header=FALSE,sep="\t")
  names(data)<-c("mean_stars","fold_enrichment","score")
  x <- 2.8
  y <- y-j 
  data2.labels <- data.frame(x = c(x),y = c(y),label = c(label))
  plot<-plot + geom_point(data=data,aes(x=mean_stars,y=fold_enrichment),size=3,color=color, shape = 20) +
    geom_text(data= data2.labels, aes(x=x,y=y,label= label), color=color,size = 4) #,fontface="bold")
}

plot <- plot  +
  ggtitle(title) +
  scale_colour_discrete()+
  guides(colour = guide_legend(override.aes = list(size=3))) +
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank())+
  theme(plot.title = element_text(size = rel(2), face="bold", colour = "#003366"))+
  theme(strip.background=element_rect(fill='white',color='white'))+
  theme(strip.text.y = element_blank())+
  theme(axis.text=element_text(size=rel(1)))+
  xlab(label="Mean confidence score")+
  ylab(label="Fold enrichment compared with the gold standard")+
  theme(axis.title=element_text(size=rel(2)))+
  theme(legend.text = element_text(size = rel(2)))+
  theme(axis.text.x = element_text(colour = "black"))+
  theme(axis.text.y = element_text(colour = "black"))+
  theme(legend.title = element_text(size = 6,face="bold"))

print(plot)
dev.off()
