# First import vcf files

module load vcftools
module load parallel
module load htslib


# vcf-concat ../freebayes_qc/vcf_by_contig_lowcov_down/*_filtered.vcf > vcf_by_contig_lowcov_down_filtered.vcf 
# vcftools --vcf vcf_by_contig_lowcov_down_filtered.vcf --stdout --recode --recode-INFO-all --keep mi.indiv.txt > mi.vcf
# vcftools --vcf vcf_by_contig_lowcov_down_filtered.vcf --stdout --recode --recode-INFO-all --keep nomi.indiv.txt > nomi.vcf

ziptab(){
	f=$1
	bgzip $f
	tabix ${f}.gz
}

export -f ziptab

parallel -j 2 ziptab ::: mi.vcf nomi.vcf
