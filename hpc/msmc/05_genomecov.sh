module load bedtools2
module load samtools
module load parallel

getcov(){
	folder=$1
	prefix=$2
	bedtools genomecov -ibam ${folder}/${prefix}_merged_marked.bam | grep 'genome' > ${prefix}.genomecov_summary
}

export -f getcov

for ds in 'FI-1-3' 'MI-1-4';do
	bedtools genomecov -ibam ../freebayes/${ds}_merged_marked.bam | grep 'genome' > ${ds}.genomecov_summary
done
