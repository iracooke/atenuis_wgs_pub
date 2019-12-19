#!/bin/bash

module load picard
module load bwa
module load samtools
module load java

GENOME=aten_final_0.1.fasta

set -e

# See here 
# http://gatkforums.broadinstitute.org/gatk/discussion/6472/read-groups#latest
# For a discussion of assigning read groups etc
#
fastq2ubam(){
	input1=$1
	input2=$2
	output=$3
	sample=$4
	barcode=$5
	lane=$6

	java -Xmx8G -jar $PICARD_HOME/picard.jar FastqToSam \
    	FASTQ=$input1 \
    	FASTQ2=$input2 \
    	OUTPUT=$output.bam \
    	READ_GROUP_NAME=$sample.$barcode.$lane \
    	SAMPLE_NAME=$sample \
    	LIBRARY_NAME=$sample \
    	PLATFORM_UNIT=$barcode.$lane.$sample \
    	PLATFORM=ILLUMINA
}

markadapters(){
	input=$1
	output=$2

    java -Xmx4G -jar $PICARD_HOME/picard.jar MarkIlluminaAdapters \
    I=$input.bam \
    M=${input}_txt \
    O=$output.bam 
}

map_reads(){
	input=$1
	fasta=$2
	ubam=$3
	output=$4

	java -Xmx4G -jar $PICARD_HOME/picard.jar SamToFastq \
    I=$input.bam \
    FASTQ=/dev/stdout \
    CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true \
    |  \
    bwa mem -M -t 16 -p $fasta /dev/stdin \
    | \
    java -Xmx8G -jar $PICARD_HOME/picard.jar MergeBamAlignment \
    ALIGNED_BAM=/dev/stdin \
    UNMAPPED_BAM=$ubam \
    OUTPUT=$output.bam \
    R=$fasta CREATE_INDEX=true ADD_MATE_CIGAR=true \
    CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true \
    INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 \
    PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS 
}

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

# First index the genome
bwa index -a bwtsw $GENOME

# Make the dict
java -jar $PICARD_HOME/picard.jar CreateSequenceDictionary REFERENCE=$GENOME OUTPUT=${GENOME%.fasta}.dict


for f in *_R1_00[0-9]_[12].fastq.gz; do

	flowcell=`gunzip -c $f | head -n 1 | awk -F ':' '{print $3}'`

	splitf=(${f//_/ })

	input1=$f
	input2=${f//R1/R2}

	input=${f//_R1_/_}
	outbase=${input%.fastq.gz}

	sample=${splitf[0]}
	barcode=${splitf[1]}
	lane=${splitf[2]}

	if [ ! -f $outbase.bam ]; then
		fastq2ubam $input1 $input2 $outbase $sample $flowcell $lane
	fi

	if [ ! -f ${outbase}_markadapters.bam ]; then
		markadapters $outbase ${outbase}_markadapters
	fi

	if [ ! -f ${outbase}_mapped.bam ]; then
		map_reads ${outbase}_markadapters $GENOME $outbase.bam ${outbase}_mapped 
	fi

	if [ ! -f ${outbase}_mapped_marked.bam ]; then
		mark_duplicates ${outbase}_mapped ${outbase}_mapped_marked
	fi
	
done




