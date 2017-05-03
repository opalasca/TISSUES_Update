#Define the colors to be used for each dataset
get_colors<-function(is_label = F){
  #Names
  ds<-c("exon","gnf","rna","hpa_rna","tm", "mouse_gnf","mouse_gnfv3","mouse_rnaseq_encode","mouse_rnaseq_mit", "rat_array", "rat_rnaseq_mit", "rat_rnaseq_bodymap", "pig_array", "pig_rnaseq_aarhus", "pig_rnaseq_wur")
  if(is_label){
    ds<-c("Human Exon Array","Human GNF","Human RNA-seq atlas","Human HPA RNA-seq","Text-mining", "Mouse GNF", "Mouse GNF V3","Mouse RNA-seq Encode", "Mouse RNA-seq MIT", "Rat Array", "Rat RNA-seq MIT","Rat RNA-seq BodyMap", "Pig Array", "Pig RNA-seq Aarhus", "Pig RNA-seq WUR","Human","Mouse","Rat","Pig")
    }
  #assign colors to each dataset
  col<-c("#3cb371","#9acd32", "#2e8b57","#006400","#556b2f", "#4169e1","#4682b4","#1e90ff","#191970","#cd853f","#daa520","#a0522d","#ffc0cb","#ffa07a","#db7093","white","white","white","white" )
#            exon      gnf      rna      hparna                "m_gnf", "m_gnfv3","mencode", "m_mit", "r_array", r_mit   rat_bmap    p_array    p_aarhus    p_wur
            
   names(col) <- ds  
  return(col)
}

splitvalues<-function(row){
  x<-as.character(row)
  substr(x,1,regexpr("-",x)-1)
}

starts_with <- function(vars, match, ignore.case = TRUE) {
  if (ignore.case) match <- tolower(match)
  n <- nchar(match)
  if (ignore.case) vars <- tolower(vars)
  substr(vars, 1, n) == match
}

get_data<-function(consistency_file){
  in_file<-paste(consistency_file,".tsv",sep="")
  protein_tissues<-read.table(in_file,sep="\t",header=FALSE)
  names(protein_tissues)<-c("dataset","score","stars","tissue","proteins","human","mouse","rat","pig")
  #names(protein_tissues)<-c("dataset","stars","score","tissue","proteins","human","mouse","rat","pig")
  return(protein_tissues)
}

filter_dataset<-function(data,dataset_name){
    filtered_data<-data[data$dataset == dataset_name, ]
  return(filtered_data)
}

calc_pairwise_corr <- function(A, B){
  p = stats::cor(A,B, use="all.obs", method="pearson")   
  return(p) 
}

# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[lower.tri(cormat)] <- NA
  diag(cormat)<-NA
  return(cormat)
}

# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[upper.tri(cormat)]<- NA
  diag(cormat)<-NA
  return(cormat)
}

#adapted from http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization
draw_heatmap<-function(mat_data,title,output_file,min_val,w,h,bw,bh){
  #h=1000;w=1000
  if (min_val<0){
    col=c("blue4","whitesmoke","brown3")
  }
  else{
    col=c("whitesmoke","brown3")
  } 
  upper_tri <- get_upper_tri(mat_data)
  melted_mat_data<-melt(upper_tri,na.rm=TRUE)
  low<-min(melted_mat_data[,3])
  print(low)
  #if (min_val!=0){
    low=min_val
  #} 
  mid<-low+(1-low)/2
  plot <- ggplot(melted_mat_data, aes(Var1,ordered(Var2,levels=rev(sort(unique(Var2)))))) +   
    geom_tile(aes(fill = value))+
    labs(caption=title) + 
    theme(plot.caption = element_text(hjust=0.5, face="bold", size=rel(1.2)))+
    scale_fill_gradientn(colours=col, 
                         limits = c(low,1),
                         #space = "Lab", 
                         name="Pearson's\nCorrelation") +
    theme(axis.text.x = element_text(angle = 50, size=rel(1.2), vjust = 0.5, hjust = 0),
          axis.text.y = element_text(size=rel(1.2)))+
    scale_x_discrete(position = "top")+
    scale_y_discrete(position = "right")+
    coord_fixed()
    
  plot <- plot + 
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.justification = c(0, 1),
      legend.position = c(0, 0.7),
      legend.direction = "vertical")+
      guides(fill = guide_colorbar(barwidth = rel(bw), barheight = rel(bh),
                                 title.position = "top", title.hjust = 0))

  png(paste(output_file,"png",sep="."), height=h, width=w)      
  print(plot)
  dev.off()             
  
  pdf(paste(output_file,"pdf",sep="."), height=h/75, width=w/75)      
  print(plot)
  dev.off()             
  
  return(plot)
}  


# from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
# Multiple plot function
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#http://stackoverflow.com/questions/41127054/rotate-a-matrix-45-degrees-and-visualize-it-using-ggplot
rotate <- function(df, degree) {
  dfr <- df
  degree <- pi * degree / 180
  l <- sqrt(df$start1^2 + df$start2^2)
  teta <- atan(df$start2 / df$start1)
  dfr$start1 <- round(l * cos(teta - degree))
  dfr$start2 <- round(l * sin(teta - degree))
  return(dfr)
}
