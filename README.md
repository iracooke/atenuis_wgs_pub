# Acropora tenuis whole genome resequencing

Commands and R code to accompany the paper

> Signatures of selection in the Acropora tenuis holobiont reveal complex adaptations to inshore environments driven by Holocene climate change

### Raw Data Processing

### Downstream Analysis and Plots

All of the following sections are provided as processed markdown files.  Clicking the link should display a web readable page with text, a few select commands and plots and tables. The code used to generate these pages is provided in the corresponding `.Rmd` file. If you would like to run the code in these files yourself you will need to;

1. Checkout this repository 
2. Download additional data not hosted on github due to size constraints
```bash
wget something
tar -zxvf something
```
3. Open the project file `atenuis_wgs_pub.Rproj` in RStudio and open the desired file (eg 01_population_structure.Rmd).  After installing any require R packages the code should run and produce plots and tables identical to those shown on the corresponding web link below.

- [Population Structure](01_population_structure.md) : Rmd file [01_population_structure.md](01_population_structure.md)
- [Estimating Mutation Rate](02_mutation_rates.md) : Rmd file [02_mutation_rates.md](02_mutation_rates.md)
- [Demographic History with MSMC](03_msmc.md)
- [Demographic History with dadi](04_dadi.md)
- [Selective Sweep Detection with SweepFinder 2](05_sf2_thresholds.md)
- [Interpretation of SweepFinder results](06_sf2.md)
- [Neutrality and Diversity statistics with ANGSD](07_popgen_stats.md)
- [Symbiont Profiles with kraken](08_metagenome.md)
- [Analysis of interspersed repeats](09_repeats.md)

