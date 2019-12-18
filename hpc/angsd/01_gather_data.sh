# Make bamlists

# Exclude high coverage samples to avoid any potential biases
#
for pop in DI FI PI PR MI;do
	ls ../gatk3/${pop}*merged_marked.bam | grep -v 'MI-1-4_' | grep -v 'FI-1-3_' > ${pop}_bamlist.txt
done

# All northern pops
#
cp DI_bamlist.txt north_bamlist.txt

for pop in FI PI PR;do
	cat ${pop}_bamlist.txt >> north_bamlist.txt
done

ln -s ../freebayes_qc/aten_final_0.1.fasta .

# This makes a list of sites suitable for analysis
# These are sites that passed our filters
#
for f in ../freebayes_qc/vcf_by_contig_lowcov_down/*_filtered.vcf; do cat $f | awk '!/^#/{print $1,$2}';done > good_sites.txt

