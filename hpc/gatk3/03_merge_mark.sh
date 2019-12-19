#!/bin/bash

module load picard
module load bwa
module load samtools
module load java

GENOME=aten_final_0.1.fasta

set -e


mark_duplicates(){
        input=$1
        output=$2

    java -Xmx32G -jar $PICARD_HOME/picard.jar MarkDuplicates \
    INPUT=$input.bam \
    OUTPUT=$output.bam \
    METRICS_FILE=${input}_markduplicates_txt \
    CREATE_INDEX=true \
    TMP_DIR=/tmp/sci-irc \
    MAX_RECORDS_IN_RAM=5000000 \
    MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=800
}


samples=$(ls *.fastq.gz | awk -F '_S' '{print $1}' | sort -u | tr '\n' \ '')

for sample in $samples;do
    if [ ! -f ${sample}_merged_marked.bam ];then
        echo "Merging $sample"
        samtools merge ${sample}_merged.bam ${sample}_*mapped_marked.bam
        echo "Marking duplicates for $sample"
        mark_duplicates ${sample}_merged ${sample}_merged_marked
    fi
done
