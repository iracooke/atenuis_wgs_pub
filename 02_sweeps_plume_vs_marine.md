Selective Sweeps
================

Calculating sweep intersections

We are looking for sweeps where we see a strong signal in both plume populations and a weak or no signal in the marine. (And vice versa) This is done in two steps.

1.  First define sweeps using a fairly relaxed minimal threshold of 15. This allows us to select all possible sweep regions.

2.  Find regions where there is a sweep in common (overlap &gt; 50%) between plume sites and require that both plume sites have some minimal score value (combination of LR and length)

3.  Now for those sweeps identified in 2 calculate the score (area under the sweep interval) for correspoding marine sites, and calculate the score differential between plume and marine for each of these sweeps.

### Steps 1 and 2

``` bash
TR=15 # LR Threshold
cd raw_data/SF2/
bedtools intersect -wo -f 0.5 -r  -a di_15_sweeps.gff -b   pr_15_sweeps.gff | awk '($5-$4)>2 && ($14-$13)>2 && $6>100 && $15>100 {print $0}' > common_plume.tsv
bedtools intersect -wo -f 0.5 -r  -a fi_15_sweeps.gff -b   pi_15_sweeps.gff | awk '($5-$4)>2 && ($14-$13)>2 && $6>100 && $15>100 {print $0}' > common_marine.tsv
```

### Step 3

Load sweeps from step 1

``` r
plume_sweeps_raw <- read_tsv("raw_data/SF2/common_plume.tsv", col_names = c(paste(gff_names(),"_a",sep=""),paste(gff_names(),"_b",sep=""),"overlap"), col_types = cols())
marine_sweeps_raw <- read_tsv("raw_data/SF2/common_marine.tsv", col_names = c(paste(gff_names(),"_a",sep=""),paste(gff_names(),"_b",sep=""),"overlap"), col_types = cols())
```

Convert plume sweeps into a more manageable form

``` r
get_sweep_intervals <- function(raw){
  raw %>% 
    mutate(start = pmin(start_a,start_b)) %>% 
    mutate(end = pmin(end_a,end_b)) %>%
    mutate(contig_offset = as.numeric(contig_offsets[seqid_a])) %>% 
    select(scaffold = seqid_a,contig_offset,score_a, score_b, overlap, start, end)
}

contig_offsets <- load_contig_offsets()

plume_sweeps <- plume_sweeps_raw %>% get_sweep_intervals()
marine_sweeps <- marine_sweeps_raw %>% get_sweep_intervals()
```

Load the raw sweepfinder scores

``` r
folder=""
sf2_all <- load_sf2(folder) %>% 
  arrange(desc(length)) %>% 
  mutate(contig_band = as.factor(contig_index %% 2))
```

For each sweep in plume sweeps we wish to calculate the sum of scores of each population in that interval

``` r
cachefile <- "cache/marine_plume_scores.rds"

score_sweeps <- function(sweeps,cachefile){
  if (!file.exists(cachefile)){
    
    score_interval <- function(start,end,pops=c("fi","pi")){
      score <- sf2_all %>%
        filter(pop %in% pops) %>% 
        filter((location >= start) & (location <= end)) %>% 
        summarise(score=sum(LR))
      as.numeric(score)
    }
    
    sweep_scores <- sweeps %>% 
      mutate(start_location = start+contig_offset) %>% 
      mutate(end_location = end+contig_offset) %>% 
      mutate(fi_score = 
               list(start_location,end_location) %>% pmap_dbl(score_interval,c("fi"))
      ) %>% 
      mutate(di_score = 
               list(start_location,end_location) %>% pmap_dbl(score_interval,c("di"))
      ) %>% 
      mutate(pi_score = 
               list(start_location,end_location) %>% pmap_dbl(score_interval,c("pi"))
      ) %>% 
      mutate(pr_score = 
               list(start_location,end_location) %>% pmap_dbl(score_interval,c("pr"))
      )  
    
    write_rds(sweep_scores, path = cachefile)
    
  } else {
    
    sweep_scores <- read_rds(cachefile)
    
  }
}

both_sweeps <- rbind(marine_sweeps,plume_sweeps)

marine_plume_scores <- both_sweeps %>% score_sweeps(cachefile)
```

OK now we can filter on the difference between marine and plume scores and vice versa

``` r
tt <- function(m1,m2,p1,p2){
  t.test(c(m1,m2),c(p1,p2))$p.value
}

test_marine_vs_plume <- function(scores){
  scores %>% 
    mutate(p_marine_vs_plume = 
             list(fi_score,pi_score,di_score,pr_score) %>%  pmap_dbl(tt)) %>% 
    mutate(fc_marine_vs_plume = 
             (di_score+pr_score)/(fi_score+pi_score)) 
}

sweeps_stronger_in_plume <- marine_plume_scores %>% 
  test_marine_vs_plume() %>% 
  filter((p_marine_vs_plume < 0.1) & (fc_marine_vs_plume > 2))

sweeps_stronger_in_marine <- marine_plume_scores %>% 
  test_marine_vs_plume() %>% 
  filter((p_marine_vs_plume < 0.1) & (fc_marine_vs_plume < 0.5))
```

Now save results to BED 3 column format

``` r
write_tsv(sweeps_stronger_in_plume %>% select(scaffold,start,end) %>% arrange(scaffold, start), path = "raw_data/SF2/sweeps_stronger_in_plume.bed",col_names = FALSE)
write_tsv(sweeps_stronger_in_marine %>% select(scaffold,start,end) %>% arrange(scaffold, start), path = "raw_data/SF2/sweeps_stronger_in_marine.bed",col_names = FALSE)
```

In order to generate plots with sweep details we use bedtools merge to define regions of interest that may include multiple sweeps

``` bash
cd raw_data/SF2
bedtools merge -d 300000 -i sweeps_stronger_in_plume.bed > sweeps_stronger_in_plume_merged.bed
bedtools merge -d 300000 -i sweeps_stronger_in_marine.bed > sweeps_stronger_in_marine_merged.bed
```

``` r
sweeps_stronger_in_plume_merged <- read_tsv("raw_data/SF2/sweeps_stronger_in_plume_merged.bed",col_names = c("scaffold","start","end"), col_types = cols()) %>% 
  rownames_to_column("row") %>% 
  unite("region",scaffold,row,sep = "_",remove = FALSE) %>% select(-row)

sweeps_stronger_in_marine_merged <- read_tsv("raw_data/SF2/sweeps_stronger_in_marine_merged.bed",col_names = c("scaffold","start","end"), col_types = cols()) %>% 
  rownames_to_column("row") %>% 
  unite("region",scaffold,row,sep = "_",remove = FALSE) %>% select(-row)
```

Read in gene annotation information

``` r
aten_gff <- read_tsv("raw_data/genome/aten_0.11.maker_post_001.genes.gff",
                     col_names = c("seqid","source","type","start","end","score","strand","phase","attributes"), col_types = cols()) %>% 
  mutate(gene_id = str_match(attributes,pattern = "ID=([^;]+)")[,2])

annotations <- read_tsv("raw_data/genome/annotation_table.tsv", col_types = cols()) %>% 
  mutate(gene_id = str_extract(aten_id,pattern = "aten_0.1.m[0-9]+.[0-9]+"))
```

Use fuzzyjoin to join together SF2, sweep region and annotation data

``` r
library(fuzzyjoin)

join_sweep_data_to_regions <- function(sf2_data,sweeps_data,gff_data,annotations){
  sf2_data %>% 
    mutate(start=location - offset) %>% 
    mutate(end=location - offset) %>% 
    genome_right_join(sweeps_data, by = c("scaffold","start","end"), maxgap = 20000) %>% 
    mutate(position = start.x) %>% select(-end.x,-start.x,-scaffold.x) %>% 
    rename_at(vars(matches('.y$')), funs(sub(".y","",.))) %>% 
    genome_left_join(gff_data, by = c("scaffold"="seqid","start","end")) %>%  filter(type %in% c("gene",NA)) %>% 
    rename_at(vars(matches('.y$')), funs(sub(".y","_gene",.))) %>%
    rename_at(vars(matches('.x$')), funs(sub(".x","",.))) %>% 
    left_join(annotations,by="gene_id")
}


best_name <- function(gene_name,acc,id){
  if ( !is.na(gene_name)){ return (gene_name)}
  if ( !is.na(acc)){return(acc)}
  return(id)
}

sweeps_stronger_in_plume_region_data <- sf2_all %>% join_sweep_data_to_regions(sweeps_stronger_in_plume_merged,aten_gff,annotations)
sweeps_stronger_in_marine_region_data <- sf2_all %>% join_sweep_data_to_regions(sweeps_stronger_in_marine_merged,aten_gff,annotations)

sweeps_stronger_in_plume_region_gene_data <- sweeps_stronger_in_plume_region_data %>% 
  select(region,start_gene,end_gene,`Protein names`,`Gene names`,saccver, gene_id) %>% 
  distinct() %>% 
  mutate(label = list(`Gene names`,saccver, gene_id) %>% pmap_chr(best_name))

sweeps_stronger_in_marine_region_gene_data <- sweeps_stronger_in_marine_region_data %>% 
  select(region,start_gene,end_gene,`Protein names`,`Gene names`,saccver, gene_id) %>% 
  distinct() %>% 
  mutate(label = list(`Gene names`,saccver, gene_id) %>% pmap_chr(best_name))


write_tsv(sweeps_stronger_in_plume_region_gene_data,"results/plume_only_sweep_genes.tsv")
write_tsv(sweeps_stronger_in_marine_region_gene_data,"results/marine_only_sweep_genes.tsv")
```

Stronger in Plume

``` r
library(ggrepel)
library(ggpubr)
```

    ## Loading required package: magrittr

    ## 
    ## Attaching package: 'magrittr'

    ## The following object is masked from 'package:purrr':
    ## 
    ##     set_names

    ## The following object is masked from 'package:tidyr':
    ## 
    ##     extract

``` r
plot_sweeps <- function(sweep_region_data,sweep_region_gene_data){
  gene_track_y <- 50
  
  ggplot(sweep_region_data,aes(y=LR,x=position)) + 
    geom_point(aes(color=pop)) + 
    scale_color_manual(values = c("fi"="#0044FC","pi"="#77D7FE","di"="#93510C","pr"="#FFD27D","mi"="#148D13")) +
    geom_segment(data=sweep_region_gene_data,aes(x=start_gene,xend=end_gene,y=-gene_track_y,yend=-gene_track_y), size=2) +
    geom_text_repel(data=sweep_region_gene_data,aes(label=label,x=(start_gene+end_gene)/2,y=-gene_track_y), size = 3) +
    geom_segment(data=sweep_region_gene_data,aes(x=(start_gene+end_gene-100)/2,xend=(start_gene+end_gene+100)/2,y=-gene_track_y,yend=-gene_track_y), size = 3, color="red") +
    facet_wrap(~region, scales = "free_x", ncol = 3) +
    theme_pubclean()
}

po_plot <- plot_sweeps(sweeps_stronger_in_plume_region_data,sweeps_stronger_in_plume_region_gene_data)
ggsave(po_plot,filename = "figures/plume_only.png", width=18,height = 24)
```

    ## Warning: Removed 6 rows containing missing values (geom_segment).

    ## Warning: Removed 6 rows containing missing values (geom_text_repel).

    ## Warning: Removed 6 rows containing missing values (geom_segment).

``` r
mo_plot <- plot_sweeps(sweeps_stronger_in_marine_region_data,sweeps_stronger_in_marine_region_gene_data)
ggsave(mo_plot,filename = "figures/marine_only.png", width=18,height = 24)
```

    ## Warning: Removed 1 rows containing missing values (geom_segment).

    ## Warning: Removed 1 rows containing missing values (geom_text_repel).

    ## Warning: Removed 1 rows containing missing values (geom_segment).