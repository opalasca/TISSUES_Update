levels<-c("mouse","mouse_gnf","mouse_gnfv3","mouse_rnaseq_encode","mouse_rnaseq_mit","rat", "rat_array", "rat_rnaseq_mit", "rat_rnaseq_bodymap", "pig", "pig_array", "pig_rnaseq_aarhus", "pig_rnaseq_wur", "human", "exon","gnf","rna", "hpa_rna" )
labels<-c( "Mouse","Mouse GNF", "Mouse GNF V3", "Mouse RNA-seq Encode", "Mouse RNA-seq MIT", "Rat","Rat Array", "Rat RNA-seq MIT", "Rat RNA-seq BodyMap", "Pig","Pig Array", "Pig RNA-seq Aarhus", "Pig RNA-seq WUR","Human","Human Exon Array","Human GNF","Human RNA-seq atlas","Human HPA RNA-seq")

tis<-c("heart","intestine","kidney","liver","nervous system","muscle","lung","spleen","lymph node","pancreas","skin",
       "stomach","thyroid gland","bone marrow","saliva","adrenal gland","urine","blood","eye","gall bladder","bone")
tis_labels<-c("Heart","Intestine","Kidney","Liver","Nervous system","Muscle","Lung","Spleen","Lymph node","Pancreas","Skin",
       "Stomach","Thyroid gland","Bone marrow","Saliva","Adrenal gland","Urine","Blood","Eye","Gall bladder","Bone")

data<-read.table(file="datasets_major_tissues.tsv",sep="\t",header=F)
names(data)<-c("dataset","tissues","num_prots")
data<-data[data$dataset %in% levels,]
#data<-data[data$tissues %in% tis,]

#Update the number of proteins
data$dataset<-factor(data$dataset,levels=rev(levels),labels= rev(labels))
data$tissues<-factor(data$tissues,levels=tis,labels=tis_labels)
datasets<-unique(data$dataset)

#Colors
col<-get_colors(is_label = T)

fillcolScale <- scale_fill_manual(name = "",values = col,labels<-c( "Mouse","Mouse GNF", "Mouse GNF V3", "Mouse RNA-seq Encode", "Mouse RNA-seq MIT", "Rat","Rat Array", "Rat RNA-seq MIT", "Rat RNA-seq BodyMap", "Pig","Pig Array", "Pig RNA-seq Aarhus", "Pig RNA-seq WUR","Human","Human Exon Array","Human GNF","Human RNA-seq atlas","Human HPA RNA-seq"))
colScale <- scale_color_manual(name = "",values = col,labels<-c( "Mouse","Mouse GNF", "Mouse GNF V3", "Mouse RNA-seq Encode", "Mouse RNA-seq MIT", "Rat","Rat Array", "Rat RNA-seq MIT", "Rat RNA-seq BodyMap", "Pig","Pig Array", "Pig RNA-seq Aarhus", "Pig RNA-seq WUR","Human","Human Exon Array","Human GNF","Human RNA-seq atlas","Human HPA RNA-seq"))


#plot<-ggplot(data=data,aes(dataset,tissues,colour=dataset))+geom_point(shape=16,size=12)+
plot<-ggplot(data=data,aes(dataset,tissues,colour=dataset))+geom_point(shape=16,size=rel(13))+
  coord_flip()+
  scale_x_discrete()+
  #scale_x_reverse()+
  
  colScale+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  #xlab(label="Tissues and number of proteins per dataset")+
  theme(axis.text.x = element_text(angle = 90, hjust=0.5,vjust=0.5))+
  theme(axis.text.y = element_text(angle = 0,hjust=0.5,vjust=0.5))+
  theme(legend.position="none")+
  theme(line = element_blank(),
        axis.title.y=element_text(color="white"))+
  theme(axis.text.x = element_text(family= "Helvetica",colour = "black"))+
  theme(axis.text.y = element_text(family= "Helvetica",colour = "black"))+
  theme(axis.text = element_text(face="bold",family= "Helvetica",size=rel(1.2)))+
  theme(axis.title.x = element_blank())+
  theme(axis.title.y = element_blank())

png("../figures/Summary_figure_tissues.png",height=750,width=1000)
print(plot)
dev.off()

pdf("../figures/Summary_figure_tissues.pdf",height=10.5,width=10.5*4/3)
print(plot)
dev.off()
