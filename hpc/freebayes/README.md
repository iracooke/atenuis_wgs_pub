# Variant calling with Freebayes

Variants were called for all 146 low coverage samples and 2 high coverage samples downsampled to 3x. All sample names are listed in the file [lowcov_down_bams.txt](hpc/freebayes/lowcov_down_bams.txt). See [here](../gatk3/README.md) for details of preprocessing and read alignment for these samples.

To facilitate rapid variant calling we divided the *A. tenuis* genome into 100kb regions `fasta_generate_regions.py` script provided with freebayes version 1.1.1
```bash
fasta_generate_regions.py aten_final_0.1.fasta.fai 100000 > aten_regions.txt
```

Bam files for all samples were indexed with samtools (eg;)
```bash
samtools index DI-1-10_merged_marked.bam
```

Freebayes was then run as follows to call variants
```bash
freebayes-parallel aten_regions.txt 46 -f aten_final_0.1.fasta \
	-L lowcov_down_bams.txt \
	--genotype-qualities \
	-E -1 \
	-m 30 -q 20 \
	-K --strict-vcf > lowcov_down.vcf
```
This sets options  `-E -1` to avoid using the haplotype caller because this produces output that is difficult to parse for downstream SNP based scripts. It sets input base call and mapping qualities to Phred 20 and 30 respectively to ensure that only reliable reads are included in the read counts against each allele (used for allele frequency calculations). The `-K` option is also important for read-based allele frequency counting (for dadi inputs) because it ensures that all alleles passing input filters are output regardless of genotyping outcome.
