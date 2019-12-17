module load vcftools

thin_pop(){
	pop=$1

	vcftools --gzvcf ${pop}_filtered_mac2.vcf.gz --out ${pop}_filtered_mac2_thin \
	--thin 1000 \
	--recode --recode-INFO-all
}

thin_pop 'dadi'
