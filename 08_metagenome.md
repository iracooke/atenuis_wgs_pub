Symbiont Profiles
================

The relative abundance of major clades (genera) of Symbiodiniaceae was
profiled using [kraken](https://ccb.jhu.edu/software/kraken/) version
1.1.1 (Wood and Salzberg 2014) to classify raw reads from each sample.
To ensure that this made full use of available reads and also that it
was not affected by biased database composition we restricted our
analysis to taxa for which a complete reference genome was available.
This included the following data;

  - Clade A: *Symbiodinium microadriadicum* [from
    genbank](https://www.ncbi.nlm.nih.gov/assembly/GCA_001939145.1)
  - Clade B: *Breviolum sp.* [from
    OIST](https://marinegenomics.oist.jp/symb/download/symbB.v1.0.genome.fa.gz)
  - Clade C: (C1) *Cladocopium sp.* [from
    reefgenomics](http://symbs.reefgenomics.org/download/SymbC1.Genome.Scaffolds.fasta.gz)
  - Clade D: *Durusdinium sp.* provided courtesy of Assoc.
    Prof. Mauricio Rodriguez-Lanetty, Department of Biological Sciences
    Florida International University
  - Clade F: *Fugacium sp.* [from
    reefgenomics](http://symbs.reefgenomics.org/download/SymbF.Genome.Scaffolds.fasta.gz)

These genomes were combined with the host *A. tenuis* genome as well the
standard kraken bacterial sequences to build a kraken database (see
[07\_build\_kraken.sh](hpc/symbiodinium_profiles/07_build_kraken.sh))
using default values for kmer length (31) and minimiser length (15).

kraken was then used to classify all read pairs for all samples and the
raw outputs processed with `kraken-mpa-report`. This produces a report
in a format similar to MetaPhlAn’s output (see
[09\_run\_kraken\_genome.sh](hpc/symbiodinium_profiles/09_run_kraken_genome.sh)).

Irrespective of whether absolute read counts or proportion of reads is
used the dominant symbiont clade for all locations and all but one
sample was *Cladocopium*. A single sample from Dunk Island was dominated
by *Durusdinium*.

![](08_metagenome_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

<div id="refs" class="references">

<div id="ref-Wood2014-qp">

Wood, Derrick E, and Steven L Salzberg. 2014. “Kraken: Ultrafast
Metagenomic Sequence Classification Using Exact Alignments.” *Genome
Biol.* 15 (3): R46.

</div>

</div>
