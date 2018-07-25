library(tidyverse)



load_contig_offsets <- function(folder = "_angsd"){
  contigs <- read_tsv(paste(c("raw_data/SF2",folder,"/contig_lens.txt"),collapse = ""),col_names = c("name","length"),col_types =cols()) %>% 
    arrange(desc(length))
  
  contig_offsets <- cumsum(contigs$length)-contigs$length[1]
  names(contig_offsets) <- contigs$name
  contig_offsets
}

load_sf2 <- function(folder = "_angsd"){

  read_contig_data <- function(path,offsets){
    scaff <- path %>% basename() %>% stringr::str_extract("[^\\.]*")
    raw_c <- read_tsv(path,col_types =cols()) %>% add_column(scaffold=scaff)

    if ( is.numeric(raw_c$location) ){
      raw_c$location <- raw_c$location + offsets[scaff]
    }
    raw_c
  }

  contigs <- read_tsv(paste(c("raw_data/SF2",folder,"/contig_lens.txt"),collapse = ""),col_names = c("name","length"),col_types =cols()) %>% 
    arrange(desc(length)) %>% 
    mutate(offset = cumsum(length) - first(length)) %>% 
    mutate(contig_index = row_number())

  contig_offsets <- load_contig_offsets()
  
  cachefile = paste(c("cache/sf2",folder,"_all.rds"),collapse = "")
  
  base_folder = paste(c("raw_data/SF2",folder),collapse = "")
  
  if ( ! file.exists(cachefile) ){

    sf2_fi <- do.call(rbind,lapply(list.files(file.path(base_folder,"FI"),pattern="*[0-9].sf2",full.names = TRUE),read_contig_data,contig_offsets))
    
    sf2_pi <- do.call(rbind,lapply(list.files(file.path(base_folder,"PI"),pattern="*[0-9].sf2",full.names = TRUE),read_contig_data,contig_offsets))
    
    sf2_di <- do.call(rbind,lapply(list.files(file.path(base_folder,"DI"),pattern="*[0-9].sf2",full.names = TRUE),read_contig_data,contig_offsets))
    
    sf2_pr <- do.call(rbind,lapply(list.files(file.path(base_folder,"PR"),pattern="*[0-9].sf2",full.names = TRUE),read_contig_data,contig_offsets))
    
    sf2_mi <- do.call(rbind,lapply(list.files(file.path(base_folder,"MI"),pattern="*[0-9].sf2",full.names = TRUE),read_contig_data,contig_offsets))
    
    sf2_all <- sf2_fi %>% add_column(pop="fi") %>% 
      rbind(sf2_pi %>% add_column(pop="pi")) %>% 
      rbind(sf2_di %>% add_column(pop="di")) %>%
      rbind(sf2_pr %>% add_column(pop="pr")) %>%
      rbind(sf2_mi %>% add_column(pop="mi")) 
    
    sf2_all$condition <- sf2_all$pop
    sf2_all$condition <- sub("fi","marine",sf2_all$condition)
    sf2_all$condition <- sub("pi","marine",sf2_all$condition)
    sf2_all$condition <- sub("di","plume",sf2_all$condition)
    sf2_all$condition <- sub("pr","plume",sf2_all$condition)
    
    sf2_all$catchment <- sf2_all$pop
    sf2_all$catchment <- sub("fi","north",sf2_all$catchment)
    sf2_all$catchment <- sub("pi","south",sf2_all$catchment)
    sf2_all$catchment <- sub("di","north",sf2_all$catchment)
    sf2_all$catchment <- sub("pr","south",sf2_all$catchment)
    
    sf2_all <- sf2_all %>% left_join(contigs, by = c("scaffold" = "name"))
    
    write_rds(sf2_all,cachefile)
    
  } else {
    sf2_all <- read_rds(cachefile)
  }
  sf2_all
}

load_sweep_genes <- function(folder = "_angsd"){
  base_folder = paste(c("raw_data/SF2",folder),collapse = "")
  
  cache_file = paste(c("cache/","sweep_genes",folder,".rds"),collapse = "")
  
  if ( !file.exists(cache_file)){
    sg_files <- list.files(base_folder,"*sweep_genes.gff",full.names = TRUE)
    
    read_sg_file <- function(path){
      pop <- basename(path) %>% str_extract("[^_]+") %>% toupper()
      score <- basename(path) %>% str_extract("[0-9]+")
      data <- read_table2(path,col_names = c("seqid","source","type","start","end","score","strand","phase","attributes","num_sweeps"), col_types = cols()) %>% add_column(population=pop) %>% add_column(lr_threshold=as.integer(score))
      data
    }
    
    sg_data <- do.call(rbind,lapply(sg_files,read_sg_file))
    write_rds(sg_data,path = cache_file)
  } else {
    sg_data <- read_rds(cache_file)
  }
  sg_data
}
