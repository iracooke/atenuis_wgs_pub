library(tidyverse)

sweep_detail_plot <- function(focal_scaffold,xl,rl,gff,sf2_data,offsets = NULL, ygene=20){

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
    filter(start > xl) %>% 
    filter(end < rl) %>% 
    na.omit()
  
  anno_data$type <- factor(anno_data$type)
  
#  browser()
  
  ggplot(data,aes(x=location/1000000,y=LR)) + 
    geom_point(aes(color=pop)) +
    xlim(xl/1000000,rl/1000000)  +
    geom_segment(data=anno_data %>% filter(type=="gene"),
                 size=1,
                 aes(x=start/1e6,xend=end/1e6,y=-(10*as.numeric(type)+ygene),yend=-(10*as.numeric(type)+ygene), alpha=geneid)) + 
    guides(fill=FALSE) +
    theme_minimal() + 
    xlab("Position Mb") + 
    theme(legend.title = element_blank()) + 
    theme(text= element_text(size=16))
}

