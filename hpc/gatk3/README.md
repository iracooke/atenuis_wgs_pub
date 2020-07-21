# Read Alignment and Preprocessing

Raw sequence data was obtained as demultiplexed fastq files from the Ramaciotti Centre for Genomics.  A complete listing of all these files is provided as [raw_fastq.list](hpc/gatk3/raw_fastq.list)

Raw sequence data was processed according to GATK best practices.  Although GATK was not used for variant calling this pipeline is also useful for other variant callers because it ensures that appropriate metadata is embedded within bam files. It also incorporates adapter marking and duplicate marking, both of which are also used by Freebayes when calling as described [in the freebayes readme](https://github.com/ekg/freebayes)

The script [02_map_read.sh](02_map_read.sh) provides details of adaptor marking, read mapping and duplicate marking.  This script uses picard version 2.18.13, bwa mem version 0.7.17 and samtools 1.7

Since samples were multiplexed across multiple lanes it was necessary to merge these into a single bam file for input to Freebayes. The script [03_merge_mark.sh](hpc/gatk3/03_merge_mark.sh) performs this merging and then reruns the duplicate marking on the merged file.

Although most samples were sequenced at low (~3x) depth there were two high depth samples. For some downstream analyses we created down-sampled versions of these files using samtools to randomly select reads to achieve a desired depth of 3x. see [04_sampledeep.sh](hpc/gatk3/04_sampledeep.sh)
