## Create exclusion filters
#
# First create merged regions representing areas of the genome where a sweep in either of the plume or marine sites exists
# These files are designed to capture the entire region in the vicinity of a sweep in marine or plume so that this 
# can be used as an exclusion filter (see below)
#
bedtools sort -i <(cat fi_10_sweeps.gff pi_10_sweeps.gff) > marine_10_sweeps.gff

bedtools merge -d 20000 -i marine_10_sweeps.gff -c 6 -o max > marine_merged_10_sweeps.gff

bedtools sort -i <(cat pr_10_sweeps.gff di_10_sweeps.gff) > plume_10_sweeps.gff

bedtools merge -d 20000 -i plume_10_sweeps.gff -c 6 -o max > plume_merged_10_sweeps.gff


### Find sweeps present in both marine

bedtools intersect -wo -a fi_10_sweeps.gff -b pi_10_sweeps.gff | \
	awk -f merge_sweeps.awk > marine_both.gff

### Or in both plume

bedtools intersect -wo -a di_10_sweeps.gff -b pr_10_sweeps.gff | \
	awk -f merge_sweeps.awk > plume_both.gff

### For sweeps present in both marine subtract any that are in plume

bedtools  subtract -A -a marine_both.gff -b plume_merged_10_sweeps.gff \
	> marine_only.gff

### And vice versa

bedtools  subtract -A -a plume_both.gff -b marine_merged_10_sweeps.gff \
	> plume_only.gff


### Fine genes-in-sweeps for these

bedtools window -a aten_0.11.maker_post_001.genes.gff \
	-b marine_only.gff | grep 'gene' > marine_only_sweep_genes.gff

bedtools window -a aten_0.11.maker_post_001.genes.gff \
	-b plume_only.gff | grep 'gene' > plume_only_sweep_genes.gff
