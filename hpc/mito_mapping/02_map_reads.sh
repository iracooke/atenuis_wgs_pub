module load bwa
module load samtools
module load parallel

map_reads(){

	mitogenome=$1

	bamfile=$2

	if [ ! -f ${sample}_mitoreads.bam ];then

		sample=${bamfile%_merged_marked.bam}

		samtools sort -n -@ 8 ${bamfile} | \
		samtools fastq -F 1024 - | \
		bwa mem -t 16 ${mitogenome} - | \
		samtools view -b -F 4 - > ${sample}_mitoreads.bam
	fi

}

export -f map_reads

parallel -j 40 map_reads 'AF338425.fasta' ::: $(ls *_merged_marked.bam | tr '\n' ' ')