module load vcftools
module load parallel
module load htslib

import_pop(){
	pop=$1
	vcf-concat ../freebayes_qc/vcf_by_contig_${pop}/*_filtered_mac2.vcf > ${pop}_filtered_mac2.vcf 
}

bgzip_pop(){
	pop=$1
	bgzip ${pop}_filtered_mac2.vcf
	tabix ${pop}_filtered_mac2.vcf.gz
}

export -f import_pop
export -f bgzip_pop

import_pop 'dadi'
bgzip_pop 'dadi'
