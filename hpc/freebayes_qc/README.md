# Variant Filtering and Quality Control

For practical purposes the raw vcf file produced by freebayes was first split by contig using vcftools.  This provided a simple basis for parallelisation (see [01_extract_contigs_down.sh](01_extract_contigs_down.sh))

We then used the [mdust](https://github.com/lh3/mdust) utility to create a map of low complexity regions as follows;
```bash
mdust aten_final_0.1.fasta -c > aten_final_0.1.mdust.txt
awk '{print $1,$3,$4}' aten_final_0.1.mdust.txt > aten_final_0.1.mdust.bed
```

We then applied the following filtering steps;

1. Only biallelic variants with a minimum quality score of 30 were retained
2. Sites with very low (<0.9) or very high (>5.4) depth were excluded
3. Genotypes with quality less than 20 were set to missing
4. Only variants with at least 50% called genotypes were retained
5. Variants with bad allele balance p<0.01 were removed
6. Variants that overlapped with low complexity regions (see mdust above) were removed
7. SNPs within a window of +/-10bp of an indel were removed
8. Small indels (<50bp) and SNPs were combined to produce the final filtered set.

For dadi analyses the same process was followed except that steps 3 and 4 were not applied as these are taken care of by the vcf2dadi script. In addition variants were filtered to ensure that they had a minor allele count of at least 2. (See [04_extra_filters.sh](04_extra_filters.sh))

Although small indels were retained at this step most downstream steps (eg dadi modelling, SweepFinder analyses) used only the SNP portion of the dataset.

The script [03_filter.sh](03_filter.sh) shows the commands used to implement the steps above.

