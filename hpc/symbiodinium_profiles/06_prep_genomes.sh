# Taxonomy identifiers for various clades

# CladeA=311979
# CladeB=304322
# CladeC=293275
# CladeD=298137
# CladeF=298136
# Atenuis=70783

# Kraken FASTA headers need to look like this
#>sequence16|kraken:taxid|32630  Adapter sequence

cat aten_final_0.1.fasta | bioawk -c fastx -v taxid=70783 \
	'{printf(">sequence%s|kraken:taxid|%s %s\n%s\n",NR,taxid,$name,$seq)}' \
	> genome_krakenfa/aten.fa

cat CladeA/GCA_001939145.1_ASM193914v1_genomic.fna | bioawk -c fastx -v taxid=311979 \
	'{printf(">sequence%s|kraken:taxid|%s %s\n%s\n",NR,taxid,$name,$seq)}' \
	> genome_krakenfa/cladea.fa

cat CladeB/symbB.v1.0.genome.fa | bioawk -c fastx -v taxid=304322 \
	'{printf(">sequence%s|kraken:taxid|%s %s\n%s\n",NR,taxid,$name,$seq)}' \
	> genome_krakenfa/cladeb.fa

cat CladeC/SymbC1.Genome.Scaffolds.fasta | bioawk -c fastx -v taxid=293275 \
	'{printf(">sequence%s|kraken:taxid|%s %s\n%s\n",NR,taxid,$name,$seq)}' \
	> genome_krakenfa/cladec.fa

cat CladeD/Dtrenchii_assembly_v2.fasta | bioawk -c fastx -v taxid=298137 \
	'{printf(">sequence%s|kraken:taxid|%s %s\n%s\n",NR,taxid,$name,$seq)}' \
	> genome_krakenfa/claded.fa

cat CladeF/SymbF.Genome.Scaffolds.fasta | bioawk -c fastx -v taxid=298136 \
	'{printf(">sequence%s|kraken:taxid|%s %s\n%s\n",NR,taxid,$name,$seq)}' \
	> genome_krakenfa/cladef.fa

