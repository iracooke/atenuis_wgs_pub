Haplotype Network Statistics
================

Visual inspection of the haplotype network for symbionts shows a pattern
where plume sites tend to fall into a single haplotype whereas marine
sites tend to be more diverse. In addition there seems to be a general
pattern of differentiation these groups (plume and marine) since they
share few haplotypes (with the notable exception of pandora reef which
spans between both).

![symb\_hapnet](figures/AllSymbC1MitoConsensus_GoodCoverage_PopArt.png)

We use AMOVA (Excoffier, Smouse, and Quattro 1992) implemented in
[pegas](https://cran.r-project.org/web/packages/pegas/index.html)
(Paradis 2010) version 0.12 to test whether genetic differentiation is
significant between reefs and between plume and marine groupings.

The data used are aligned mitochondrial sequences for all 108 colonies
that had sufficient coverage (at least 5x over the *Cladocopium* (C1)
mitogenome). The alignment is trimmed to exclude all ambiguities and is
therefore 6170 bp in length.

![](10_haplonetworks_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

We define strata according to reef and also to `wq` which represents the
plume/marine grouping for reefs.

``` r
mhaps.gid <- multidna2genind(mhaps, mlst = TRUE)

sample_data <- tibble(sample=mhaps@labels) %>% 
  mutate(reef = str_extract(sample,pattern = "[A-Z]+")) %>% 
  mutate(wq = case_when(reef %in% c("DI","PR","MI") ~ "plume", reef %in% c("FI","PI") ~ "marine")) 

strata(mhaps.gid) <- sample_data %>% select(wq, reef)
setPop(mhaps.gid) <- ~reef
```

AMOVA is then performed based with reef nested within wq. This indicates
highly significant differentiation by reef but this is not significant
for the water quality (`wq`) grouping, which reflects low power due to
the small number of reefs sampled.

``` r
mhaps_dist <- dist.multidna(mhaps, pool = TRUE)
pegas::amova(mhaps_dist ~ wq/reef, data = strata(mhaps.gid), nperm = 100)
```

    ## 
    ##  Analysis of Molecular Variance
    ## 
    ## Call: pegas::amova(formula = mhaps_dist ~ wq/reef, data = strata(mhaps.gid), 
    ##     nperm = 100)
    ## 
    ##                SSD          MSD  df
    ## wq    0.0011405235 1.140524e-03   1
    ## reef  0.0001270251 4.234169e-05   3
    ## Error 0.0002958351 2.872185e-06 103
    ## Total 0.0015633836 1.461106e-05 107
    ## 
    ## Variance components:
    ##           sigma2 P.value
    ## wq    2.5715e-05  0.1089
    ## reef  2.0001e-06  0.0000
    ## Error 2.8722e-06        
    ## 
    ## Phi-statistics:
    ##   wq.in.GLOBAL (Phi_CT) reef.in.GLOBAL (Phi_ST)     reef.in.wq (Phi_SC) 
    ##               0.8407070               0.9060983               0.4105094 
    ## 
    ## Variance coefficients:
    ##        a        b        c 
    ## 19.73345 23.33669 42.42593

<div id="refs" class="references">

<div id="ref-Excoffier1992-pe">

Excoffier, L, P E Smouse, and J M Quattro. 1992. “Analysis of Molecular
Variance Inferred from Metric Distances Among DNA Haplotypes:
Application to Human Mitochondrial DNA Restriction Data.” *Genetics* 131
(2): 479–91.

</div>

<div id="ref-Paradis2010-ds">

Paradis, Emmanuel. 2010. “Pegas: An R Package for Population Genetics
with an Integrated-Modular Approach.” *Bioinformatics* 26 (3): 419–20.

</div>

</div>
