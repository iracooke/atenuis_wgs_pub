module load parallel
module load samtools
module load bcftools

make_mask(){
	c=$1
	ds=$2
	echo "Making mask for $c and $ds"
	samtools mpileup -q 20 -Q 20 -C 50 -u -r $c -f aten_final_0.1.fasta "../freebayes/${ds}_merged_marked.bam" | bcftools call -c -V indels |
./bamCaller.py 20.0 ${ds}_${c}_mask.bed.gz | gzip -c > ${ds}_${c}.vcf.gz	
}

export -f make_mask

mkdir -p workfiles

parallel --tempdir /fast/tmp -j 40 make_mask ::: $(cat contig_list_1M.txt | tr '\n' ' ') ::: 'FI-1-3' 'MI-1-4'