module load parallel
module load vcftools

RAW_VCF=lowcov_down.vcf

set -e

split_vcf_by_contig(){
	invcf=$1
	contig=$2
	mkdir -p vcf_by_contig_${invcf%.vcf}

	if [[ ! -f vcf_by_contig_${invcf%.vcf}/${contig}.vcf ]]; then
		vcftools --vcf $invcf  --chr $contig --recode --recode-INFO-all --stdout > vcf_by_contig_${invcf%.vcf}/${contig}.vcf
	fi
}

if [[ ! -f  contig_list.txt ]]; then
	cat aten_final_0.1.dict | awk '{print $2}' | sed s/SN:// | tail -n +2 > contig_list.txt
fi

export -f split_vcf_by_contig

parallel -j 40 split_vcf_by_contig $RAW_VCF ::: $(cat contig_list.txt | tr '\n' ' ')
