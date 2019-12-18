
#cat aten_final_0.1.fasta | bioawk -c fastx '{print $name,length($seq)}' > aten_chromsizes.txt


for report in "LR" "alpha";do
	./sf2_to_wig.sh mi ${report} > mi_${report}.wig 
	./sf2_to_wig.sh nomi ${report} > nomi_${report}.wig 
	wigToBigWig nomi_${report}.wig aten_chromsizes.txt nomi_${report}.bw
	wigToBigWig mi_${report}.wig aten_chromsizes.txt mi_${report}.bw
done