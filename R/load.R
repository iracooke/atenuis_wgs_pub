library(tidyverse)

load_sample_table <- function(){
  # location_type
  locations <- readxl::read_excel("raw_data/LocationDetails.xlsx")
  
  read_table("raw_data/accelerate.bamlist",col_names = c("path"),col_types = cols()) %>% 
    mutate(bamfile = str_remove(path,"../gatk3/")) %>% 
    add_column(dataset = "accelerate") %>% 
    mutate(sample_id = str_match(bamfile,"([^_]+)")[,2]) %>% 
    mutate(location_id = str_match(bamfile,"([A-z]+)")[,2]) %>% 
    left_join(locations)
}
