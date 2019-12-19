# Acropora tenuis whole genome resequencing

Commands and R code to accompany the paper

> Signatures of selection in the Acropora tenuis holobiont reveal complex adaptations to inshore environments driven by Holocene climate change

### Raw Data Processing

- [Read Alignment and Preprocessing](hpc/gatk3/README.md)
- [Variant calling](hpc/freebayes/README.md)
- [Variant filtering](hpc/freebayes_qc/README.md)
- [Host mitochondrial genome assembly](hpc/mitogenome/README.md)
- [Host mitochondrial haplotype network](hpc/mito_mapping/README.md)
- [Mitochondrial haplotype networks for symbionts](hpc/symbiodinium/README.md)


### Downstream Analysis and Plots

All of the following sections are provided as processed markdown files.  Clicking the link should display a web readable page with text, a few select commands and plots and tables. The code used to generate these pages is provided in the corresponding `.Rmd` file. If you would like to run the code in these files yourself you will need to;

1. Checkout this repository 
2. Download additional data not hosted on github due to size constraints
```bash
wget something
tar -zxvf something
```
3. Open the project file `atenuis_wgs_pub.Rproj` in RStudio and open the desired file (eg 01_population_structure.Rmd).  After installing any required R packages the code should run and produce plots and tables identical to those shown on the corresponding web link below.

- [Population Structure](01_population_structure.md) : Rmd file [01_population_structure.Rmd](01_population_structure.Rmd)
- [Estimating Mutation Rate](02_mutation_rates.md) : Rmd file [02_mutation_rates.Rmd](02_mutation_rates.Rmd)
- [Demographic History with MSMC](03_msmc.md) : Rmd file [03_msmc.Rmd](03_msmc.Rmd)
- [Demographic History with dadi](04_dadi.md) : Rmd file [04_dadi.Rmd](04_dadi.Rmd)
- [Selective Sweep Detection with SweepFinder 2](05_sf2_thresholds.md) : Rmd file [05_sf2_thresholds.Rmd](05_sf2_thresholds.Rmd)
- [Interpretation of SweepFinder results](06_sf2.md) : Rmd file [06_sf2.Rmd](06_sf2.Rmd)
- [Neutrality and Diversity statistics with ANGSD](07_popgen_stats.md) : Rmd file [07_popgen_stats.Rmd](07_popgen_stats.Rmd)
- [Symbiont Profiles with kraken](08_metagenome.md) : Rmd file [(08_metagenome.Rmd)](08_metagenome.Rmd)
- [Analysis of interspersed repeats](09_repeats.md) : Rmd file [(09_repeats.Rmd)](09_repeats.Rmd)

