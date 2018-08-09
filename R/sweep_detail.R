library(tidyverse)

sweep_detail_plot <- function(focal_scaffold,xl,rl,gff,sf2_data,offsets = NULL, ygene=0){

  if (is.null(offsets)){
    data <- sf2_data %>% 
      filter(scaffold == focal_scaffold) %>% 
      filter(condition %in% c("marine","plume")) %>% 
      select(location,LR,pop) %>%  
      filter(location > xl) %>% 
      filter(location < rl) %>% 
      na.omit()    
  } else {
    data <- sf2_data %>% 
      filter(scaffold == focal_scaffold) %>% 
      filter(condition %in% c("marine","plume")) %>% 
      mutate(location = location - offsets[focal_scaffold]) %>% 
      select(location,LR,pop) %>%  
      filter(location > xl) %>% 
      filter(location < rl) %>% 
      na.omit()
  }
  
  
  
  anno_data <- gff %>% 
    filter(seqid==focal_scaffold) %>% 
    mutate(geneid=str_extract(attributes,"m[^\\;]+")) %>% 
    select(type,start,end,geneid) %>% 
    filter(( (start > xl) & (start < rl) ) |  ( (end > xl) & (end < rl)) ) %>% 
    na.omit()
  
  anno_data$type <- factor(anno_data$type)
  
#  browser()

  p <- ggplot(data,aes(x=location/1000000,y=LR)) + 
    geom_point(aes(color=pop)) +
    geom_line(aes(color=pop)) +
    guides(fill=FALSE) +
    theme_minimal() + 
    xlab("Position Mb") + 
    theme(legend.title = element_blank()) + 
    theme(text= element_text(size=16))
    
  if ( nrow(anno_data)>0){
    p <- p + geom_segment(data=anno_data %>% filter(type=="gene"),
                          size=2,
                          color="black",
                          aes(x=max(xl/1e6,start/1e6),xend=min(rl/1e6,end/1e6),y=-(10+ygene),yend=-(10+ygene), alpha = geneid)) +
      geom_text(data=anno_data %>% filter(type=="gene"), aes(x=(start/1e6+end/1e6)/2,y=-(5+ygene),label=geneid)) 
  } 
  p
}

