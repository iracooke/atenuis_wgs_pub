
echo "BEGIN" > all.poplist.txt
grep 'CHROM' ../freebayes_qc/vcf_by_contig_lowcov_down/Sc0000000_filtered.vcf | tr '\t' '\n' | tail -n +10 | grep -v '^M' | awk '{printf("%s\tnomi\n",$1)}' >> all.poplist.txt
grep 'CHROM' ../freebayes_qc/vcf_by_contig_lowcov_down/Sc0000000_filtered.vcf | tr '\t' '\n' | tail -n +10 | grep '^M' | awk '{printf("%s\tmi\n",$1)}' >> all.poplist.txt
echo "END" >> all.poplist.txt

cat all.poplist.txt | grep -v 'MI-2-9' | grep -v 'MI-1-16' | grep -v 'MI-1-1[[:space:]]' | grep -v 'PI-1-16' | grep -v 'DI-2-4'  > clean.poplist.txt
