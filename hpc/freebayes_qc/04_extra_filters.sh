module load vcftools
module load parallel

# Rather than use a maf filter we use mac instead and require a count of 2 or better. 
# This based on observations in this paper https://www.biorxiv.org/content/biorxiv/early/2017/09/14/188623.full.pdf
# where we saw that inclusion of singleton observations was problematic but we should otherwise try to go to as low a MAF as possible
# because rare alleles form an important component of the SFS
# See also https://groups.google.com/forum/#!topic/dadi-user/sve_nNhVPjc
#
filter_mac2(){
	input=$1
	output=${input%.vcf}_mac2.vcf

	if [ ! -f $output ]; then
		vcftools --vcf $input --stdout --recode --recode-INFO-all \
		--mac 2 \
		--remove-indels \
		> $output
	fi
}

export -f filter_mac2

parallel -j 40 filter_mac2 ::: $(ls vcf_by_contig_dadi/*_filtered.vcf | tr '\n' ' ')
