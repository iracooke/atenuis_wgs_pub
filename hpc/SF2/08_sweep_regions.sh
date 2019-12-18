
#cat aten_final_0.1.fasta | bioawk -c fastx 'BEGIN{OFS="\t"}{print $name,length($seq)}' > aten_final_0.1.genome

# Big 200kb slop around each sweep region.  We tend to view them on a large scale


# bedtools merge \
# 	-i <(bedtools slop -b 100000 -i nomi_30_sweeps.gff -g aten_final_0.1.genome) |\
# 	awk '{printf("%s:%s-%s\n",$1,$2,$3)}' > nomi_30_sweep_regions.txt

for t in 30 50 100;do
	bedtools merge \
		-i <(bedtools slop -b 100000 -i nomi_${t}_sweeps.gff -g aten_final_0.1.genome) >\
		nomi_${t}_sweep_regions.txt


	bedtools merge \
		-i <(bedtools slop -b 100000 -i mi_${t}_sweeps.gff -g aten_final_0.1.genome) >\
		mi_${t}_sweep_regions.txt		
done