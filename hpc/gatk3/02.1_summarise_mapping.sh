#!/bin/bash
module load samtools

summarise_mapping(){
	um=`samtools view -f4 $1 | wc -l`
	m=`samtools view -F4 $1 | wc -l`
	echo $1 $um $m
}
echo "File" "Unmapped" "Mapped"
for f in *_mapped.bam; do summarise_mapping $f;done
