module load bwa
module load samtools
module load parallel/20170522 

map_reads(){

	mitogenome=$1

	bamfile=$2

	sample=${bamfile%_nonhost_sorted.bam}
	clade=${mitogenome%.MITO_seqs.fa}

	echo "${mitogenome} ${sample} ${clade}"

	if [ ! -f ${sample}_${clade}.bam ];then
		samtools fastq -F 1024 ${bamfile} | \
		bwa mem -t 16 ${mitogenome} - | \
		samtools view -b -F 4 - > ${sample}_${clade}.bam
	fi

}

export -f map_reads

parallel -j 40 map_reads SymbC1.MITO_seqs.fa ::: $(ls *_merged_marked.bam | tr '\n' ' ')
