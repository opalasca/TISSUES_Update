calc_correlation_common_btw_all <- function(datasets, names, common_proteins, common_tissues, title, tissues_suff, width, height, barwidth, barheight){
  cat(sprintf("%f\n", (length(common_proteins))))
  pearson_matrix<-matrix(nrow=length(datasets),ncol=length(datasets))
  for (i in 1:length(datasets)){
    for (j in 1:length(datasets)){
      data1<-datasets[[i]]
      data2<-datasets[[j]]
      #cat(sprintf("%s-%s %f  ", (names(datasets))[i], (names(datasets))[j], length(common_proteins)))
      data1<-data1[data1$human %in% common_proteins,]
      data2<-data2[data2$human %in% common_proteins,]
      data1<-data1[data1$tissue %in% common_tissues,]
      data2<-data2[data2$tissue %in% common_tissues,]
      data1$pairs<-paste(data1$human,data1$tissue,sep="-")
      data2$pairs<-paste(data2$human,data2$tissue,sep="-")
      merged<-merge(data1,data2,by.x="pairs",by.y="pairs", all=TRUE)
      merged[is.na(merged)]<-0
      #cat(sprintf("%s-%s %f  ", (names(datasets))[i], (names(datasets))[j], dim(merged[])))
      pearson_matrix[i,j]<-calc_pairwise_corr(merged$stars.x, merged$stars.y)
    }
  }
  rownames(pearson_matrix)<-names
  colnames(pearson_matrix)<-names
  pearson_hm<-paste("../figures/pearson_common_btw_all_",tissues_suff,sep="")
  plot<-draw_heatmap(pearson_matrix,title,pearson_hm,min_val,width,height,barwidth,barheight)
  return(plot)
}


calc_correlation_common_btw_pairs <- function(datasets, dataset_names, names, title, suff, width, height, barwidth, barheight){
 
  pearson_matrix<-matrix(nrow=length(datasets),ncol=length(datasets))
  for (i in 1:length(datasets)){
    for (j in 1:length(datasets)){
      name1<-dataset_names[i]
      name2<-dataset_names[j]
      data1<-datasets[[i]]
      data2<-datasets[[j]]
      if  (((name1 %in% human) & (name2 %in% human)) || ((name1 %in% mouse) & (name2 %in% mouse)) || ((name1 %in% rat) & (name2 %in% rat)) || ((name1 %in% pig) & (name2 %in% pig))){}  
      else if  (((name1 %in% human)) || ((name2 %in% human))) {
        data1$proteins<-data1$human
        data2$proteins<-data2$human
      }
      else  if  (((name1 %in% mouse)) || ((name2 %in% mouse))) {
        data1$proteins<-data1$mouse
        data2$proteins<-data2$mouse
      } 
      else  if  (((name1 %in% rat)) || ((name2 %in% rat))) {
        data1$proteins<-data1$rat
        data2$proteins<-data2$rat
      }
      common_proteins<-intersect(x=data1$proteins,y=data2$proteins)
      common_proteins<-common_proteins[common_proteins!='na' & starts_with(common_proteins,'ENS')]
      common_tissues<-intersect(x=data1$tissue,y=data2$tissue)
      #cat(sprintf("%s-%s %f %f   ", (names(datasets))[i], (names(datasets))[j], length(common_proteins), length(common_tissues)))
      data1<-data1[data1$proteins %in% common_proteins, ]
      data2<-data2[data2$proteins %in% common_proteins, ]
      data1<-data1[data1$tissue %in% common_tissues,]
      data2<-data2[data2$tissue %in% common_tissues,]
      data1$pairs<-paste(data1$proteins,data1$tissue,sep="-")
      data2$pairs<-paste(data2$proteins,data2$tissue,sep="-")
      merged<-merge(data1,data2,by.x="pairs",by.y="pairs", all=TRUE)
      merged[is.na(merged)]<-0
      pearson_matrix[i,j]<-calc_pairwise_corr(merged$stars.x, merged$stars.y)
    }
  }
  rownames(pearson_matrix)<-names
  colnames(pearson_matrix)<-names
  pearson_hm<-paste("../figures/heatmap_pearson_btw_pairs")
  titlepc<-paste("Pearson's correlation based on proteins and \ntissues common between each two\n of the ten datasets")
  plot<-draw_heatmap(pearson_matrix,title,pearson_hm,min_val,width,height,barwidth,barheight)
  return(plot)
}

calc_correlation_tissues <- function(datasets, dataset_names, names, names_concat, title, suff, width, height, barwidth, barheight){
  
  pearson_matrix<-matrix(nrow=length(names_concat),ncol=length(names_concat))
  pi_last=0; pj_last=0
  pi=0; pj=0
  for (t1 in 1:length(tissues)){
      pj_last=0    
      pi_last=pi_last+pi
      for (t2 in 1:length(tissues)){
        pi=0
        pj_last=pj_last+pj
      for (i in 1:length(datasets)){
        if ((dataset_names[[i]]=="rat_array" & tissues[t1]=="liver") || (dataset_names[[i]]=="pig_rnaseq_wur" & tissues[t1]=="heart") || (dataset_names[[i]]=="pig_rnaseq_wur" & tissues[t1]=="kidney") || (dataset_names[[i]]=="exon" & tissues[t1]=="spleen") || (dataset_names[[i]]=="exon" & tissues[t1]=="lung") || (dataset_names[[i]]=="rat_array" & tissues[t1]=="lung")) {pi=pi+1; next}
        pj=0
        for (j in 1:length(datasets)){
          if ((dataset_names[[j]]=="rat_array" & tissues[t2]=="liver") || (dataset_names[[j]]=="pig_rnaseq_wur" & tissues[t2]=="heart") || (dataset_names[[j]]=="pig_rnaseq_wur" & tissues[t2]=="kidney")|| (dataset_names[[j]]=="exon" & tissues[t2]=="spleen") || (dataset_names[[j]]=="exon" & tissues[t2]=="lung") || (dataset_names[[j]]=="rat_array" & tissues[t2]=="lung")) {pj=pj+1; next}
          name1<-dataset_names[i]
          name2<-dataset_names[j]
          data1<-datasets[[i]]
          data2<-datasets[[j]]
          if  (((name1 %in% human) & (name2 %in% human)) || ((name1 %in% mouse) & (name2 %in% mouse)) || ((name1 %in% rat) & (name2 %in% rat)) || ((name1 %in% pig) & (name2 %in% pig))){}  
          else if  (((name1 %in% human)) || ((name2 %in% human))) {
            data1$proteins<-data1$human
            data2$proteins<-data2$human
          }
          else  if  (((name1 %in% mouse)) || ((name2 %in% mouse))) {
            data1$proteins<-data1$mouse
            data2$proteins<-data2$mouse
          } 
          else  if  (((name1 %in% rat)) || ((name2 %in% rat))) {
            data1$proteins<-data1$rat
            data2$proteins<-data2$rat
          }
          common_proteins<-intersect(x=data1$proteins,y=data2$proteins)
          common_proteins<-common_proteins[common_proteins!='na' & starts_with(common_proteins,'ENS')]
          #cat(sprintf("%s-%s %s-%s %f   ", (names(datasets))[i], (names(datasets))[j], tissues[t1], tissues[t2], length(common_proteins), length(common_tissues)))
          data1<-data1[data1$proteins %in% common_proteins, ]
          data2<-data2[data2$proteins %in% common_proteins, ]
          data1<-data1[data1$tissue==tissues[t1],]
          data2<-data2[data2$tissue==tissues[t2],]
          data1$pairs<-paste(data1$proteins,data1$tissue,sep="-")
          data2$pairs<-paste(data2$proteins,data2$tissue,sep="-")
          merged<-merge(data1,data2,by.x="pairs",by.y="pairs", all=TRUE)
          merged[is.na(merged)]<-0
          pc<--1
          if (nrow(merged)!=0){
            pc<-calc_pairwise_corr(merged$stars.x, merged$stars.y)
          }
          pearson_matrix[i+(t1-1)*length(datasets)-pi-pi_last, j+(t2-1)*length(datasets)-pj-pj_last]<-pc
          #cat(sprintf("[%d,%d]=%f  ",i+(t1-1)*length(datasets)-pi-pi_last, j+(t2-1)*length(datasets)-pj-pj_last, pc))
        }
      }
    }
  }
  rownames(pearson_matrix)<-names_concat
  colnames(pearson_matrix)<-names_concat
  pearson_hm<-paste("../figures/heatmap_pearson_tissues_by_pairs")
  plot<-draw_heatmap(pearson_matrix,title,pearson_hm,-1,width,height,barwidth,barheight)
  return(plot)
}


#Read data
all_data<-get_data("pairs_major_tissues_orthologs_filtered")
#all_data<-get_data("pairs_major_tissues_orthologs")
min_val=0
protein_tissues<-all_data

#datasets by organism/technology  
exon_data=NULL;gnf_data=NULL;mouse_gnf_data=NULL;mouse_gnfv3_data=NULL; rat_array_data=NULL; pig_array_data=NULL;rna_data=NULL;hparna_data=NULL;mouse_rnaseq_encode_data=NULL;mouse_rnaseq_mit_data=NULL; rat_rnaseq_bodymap_data=NULL; rat_rnaseq_mit_data=NULL; pig_rnaseq_aarhus_data=NULL; pig_rnaseq_wur_data=NULL
datasets<-list(exon=exon_data,gnf=gnf_data,mouse_gnf=mouse_gnf_data,mouse_gnfv3=mouse_gnfv3_data, rat_array=rat_array_data, pig_array=pig_array_data,rna=rna_data,hparna=hparna_data,mouse_rnaseq_encode=mouse_rnaseq_encode_data,mouse_rnaseq_mit=mouse_rnaseq_mit_data, rat_rnaseq_bodymap=rat_rnaseq_bodymap_data, rat_rnaseq_mit=rat_rnaseq_mit_data, pig_rnaseq_aarhus=pig_rnaseq_aarhus_data, pig_rnaseq_wur=pig_rnaseq_wur_data)
dataset_names<-list("exon","gnf","mouse_gnf","mouse_gnfv3","rat_array","pig_array","rna","hpa_rna","mouse_rnaseq_encode","mouse_rnaseq_mit","rat_rnaseq_bodymap","rat_rnaseq_mit", "pig_rnaseq_aarhus", "pig_rnaseq_wur")
names<-c("Human Exon Array","Human GNF","Mouse GNF","Mouse GNF V3", "Rat Array","Pig Array","Human RNA-seq atlas","Human HPA RNA-seq", "Mouse RNA-seq Encode", "Mouse RNA-seq MIT", "Rat RNA-seq BodyMap", "Rat RNA-seq MIT", "Pig RNA-seq Aarhus", "Pig RNA-seq WUR") 

human<-list("exon","gnf","rna","hpa_rna")
mouse<-list("mouse_gnf","mouse_gnfv3","mouse_rnaseq_encode","mouse_rnaseq_mit")
rat<-list("rat_array","rat_rnaseq_bodymap","rat_rnaseq_mit")
pig<-list("pig_array", "pig_rnaseq_aarhus", "pig_rnaseq_wur")

###############################################################################
#         Correlation analysis taking proteins common to all datasets         #
###############################################################################

protein_tissues_common_btw_all<-protein_tissues[protein_tissues$human!="na" & protein_tissues$mouse!="na" & protein_tissues$rat!="na" & protein_tissues$pig!="na", ]
for (i in 1:length(datasets)){
  datasets[[i]]<-filter_dataset(protein_tissues_common_btw_all,dataset_names[[i]])
}

#without rat array
datasets_to_remove_9<-c("rat_array")
names_to_remove_9<-c("Rat Array")
datasets_9<-datasets[!(names(datasets) %in% datasets_to_remove_9)]
names_9<-names[!(names %in% names_to_remove_9)]
#without rat array and pig wur
datasets_to_remove_8<-c("rat_array", "pig_rnaseq_wur")
names_to_remove_8<-c("Rat Array", "Pig RNA-seq WUR")
datasets_8<-datasets[!(names(datasets) %in% datasets_to_remove_8)]
names_8<-names[!(names %in% names_to_remove_8)]

#find the set of common proteins to all datasets
common_proteins_human<-intersect(intersect(x=intersect(x=datasets[[1]]$proteins,y=datasets[[2]]$proteins),y=datasets[[7]]$proteins),y=datasets[[8]]$proteins)
common_proteins_mouse<-intersect(intersect(x=intersect(x=datasets[[3]]$human,y=datasets[[4]]$human),y=datasets[[9]]$human),y=datasets[[10]]$human)
common_proteins_rat<-intersect(x=datasets[[11]]$human,y=datasets[[12]]$human)
common_proteins_pig<-intersect(x=intersect(x=datasets[[6]]$human,y=datasets[[13]]$human),y=datasets[[14]]$human)
common_proteins_all<-intersect(intersect(x=intersect(x=common_proteins_human,y=common_proteins_mouse),y=common_proteins_rat),y=common_proteins_pig)

p1<-calc_correlation_common_btw_all(datasets_9, names_9, common_proteins_all, c("liver"), "Liver", "liver",1000,1000,1.5,10)
p2<-calc_correlation_common_btw_all(datasets_9, names_9, common_proteins_all, c("nervous system"), "Nervous system", "nerv_sys",1000,1000,1.5,10)
p3<-calc_correlation_common_btw_all(datasets_8, names_8, common_proteins_all, c("heart"), "Heart", "heart",1000,1000,1.5,10)
p4<-calc_correlation_common_btw_all(datasets_8, names_8, common_proteins_all, c("kidney"), "Kidney", "kidney",1000,1000,1.5,10)

pearson_multi<-paste("../figures/pearson_indiv_tissues.png")
png(pearson_multi, height=1000, width=1000)   
multiplot(p3, p1, p4, p2, cols=2)
dev.off() 

pearson_multi<-paste("../figures/pearson_indiv_tissues.pdf")
pdf(pearson_multi, height=13, width=13)   
#pdf(pearson_multi)   
multiplot(p3, p1, p4, p2, cols=2)
dev.off() 

###############################################################################
#         Correlation analysis taking proteins common between pairs of datasets         #
###############################################################################

for (i in 1:length(datasets)){
  #datasets[[i]]<-filter_dataset_tissue(protein_tissues,dataset_names[[i]],c(),not_in=TRUE,c(2,3,4,5,6,7,8,9))
  datasets[[i]]<-filter_dataset(protein_tissues,dataset_names[[i]])
}
p<-calc_correlation_common_btw_pairs(datasets, dataset_names, names,"","",500,500,3,10)

###############################################################################
#         Correlation across different tissues within and between datasets         #
###############################################################################



#tissues<-c( "heart","kidney","liver","nervous system")
#tissue_names<-c("Heart","Kidney","Liver","Nervous system")

tissues<-c( "heart","kidney","liver","spleen","lung", "nervous system")
tissue_names<-c("Heart","Kidney","Liver","Spleen","Lung","Nervous system")

names_concat<-apply(expand.grid(names,tissue_names), 1, paste, collapse=" - ")
to_remove<-c("Rat Array - Liver", "Pig RNA-seq WUR - Heart", "Pig RNA-seq WUR - Kidney", "Human Exon Array - Spleen", "Human Exon Array - Lung","Rat Array - Lung")
names_concat<-names_concat[!(names_concat %in% to_remove)]

p<-calc_correlation_tissues(datasets, dataset_names, names, names_concat, "","",1500,1500,2,20)
