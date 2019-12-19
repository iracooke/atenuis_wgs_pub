# Host mitochondrial genotyping

Since our denovo assembly of the A. tenuis mitochondrial genome was extremely similar to the existing [sequence on genbank](https://www.ncbi.nlm.nih.gov/nuccore/AF338425) we used the reference assembly for mapping and genotype calling of all samples.

We obtained mitochondrial sequences for individual low coverage samples as follows;

- All host reads (except duplicates) for each sample were mapped to the mito reference using bwa mem (version 0.7.17).  This produced a small bam file containing only mito reads for each sample. See the script [02_map_reads.sh](02_map_reads.sh) for the exact commands used.  Most samples had sufficient reads to produce around 100x coverage of the mitochondrial genome with the lowest at 23x and highest at 900x. 
- Mitochondrial reads were then used to call variants against the reference sequence using samtools mpileup (version 1.7), followed by bcftools (version 1.9) to call a consensus sequence for the sample. (See script [03_call_consensus.sh](03_call_consensus.sh))
- Since consensus sequences were all the same length no alignment was necessary and they were simply concatenated to produce an alignment. 
- This set of aligned mitochondrial sequences was loaded into geneious and trimmed to remove four regions where one or more sequences had missing bases. The resulting trimmed alignment was 18163bp in length.
- The trimmed alignment was used as a set of haplotypes for network visualisation with [PopArt](http://popart.otago.ac.nz/index.shtml)

The resulting haplotype network is shown below

![host_mito_network](../../figures/AllSamplesTenuisMitoConsensus_PopArt.png)
