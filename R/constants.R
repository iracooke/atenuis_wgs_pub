library(RColorBrewer)

gff_names <- function(){
  c("seqid","source","type","start","end","score","strand","phase","attributes")  
}

site_labels <- function(){
  c("MI"="Magnetic Island","PR"="Pandora Reef","DI"="Dunk Island","PI"="Pelorus Island","FI"="Flinders Island")
}


site_order <- function(){
  c('MI'=1,'PI'=3,'DI'=4,'PR'=2,'FI'=5)
}

site_colors <- function(){
  clrs <- c(brewer.pal(5,"YlOrBr")[3:5],brewer.pal(3,"Blues")[2:3])
  names(clrs) <- c("MI","PR","DI","PI","FI")
  clrs
}




