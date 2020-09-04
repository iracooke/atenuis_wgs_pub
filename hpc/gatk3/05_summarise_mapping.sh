#!/bin/bash
module load samtools
module load parallel

summarise_mapping(){
	# um=`samtools view -f4 $1 | wc -l`
	# m=`samtools view -F4 $1 | wc -l`
	# echo $1 $um $m
	fs=$(samtools flagstat $1)
	echo $1 $fs
}
#echo "File" "Unmapped" "Mapped"
#for f in *_mapped.bam; do summarise_mapping $f;done

export -f summarise_mapping

parallel -j 40 summarise_mapping ::: *_merged_marked.bam > 05_summarise_mapping.log

cat 05_summarise_mapping.log | awk 'OFS="\t"{print $1,$28}' | sed 's/(//' | sed 's/%//' > mapping_rates.tsv