# We use bedtools to find values of stats overlapping sweeps
# First we import the sweeps file 

#ln ../SF2/nomi_10_sweeps.gff .


# In order to use bedtools we need to convert the relevant stats to bed
#
#
# Conversion for thetas
#
for thetaf in *.thetasWindow.gz.pestPG;do
	pop=${thetaf%_theta.thetasWindow.gz.pestPG}

	cat $thetaf | awk '{OFS="\t";print $2,$3-1,$3,"TajimaD",$9}' | grep -v 'Chr' > ${pop}_TajimaD.bed

	bedtools window -w 10000 -a nomi_10_sweeps.gff -b ${pop}_TajimaD.bed | \
		awk '{OFS="\t";print $1,$4,$6,$12,$14}'  > ${pop}_Tajima_sweeps.tsv

done

#
# Conversion for Fst
#
for fst in *.fst_windows.txt;do
	pp=${fst%.fst_windows.txt}
	cat $fst | awk '{OFS="\t";print $2,$3,$3,"Fst",$5}' | grep -v 'chr' > ${pp}_Fst.bed
	bedtools window -w 10000 -a nomi_10_sweeps.gff -b ${pp}_Fst.bed |\
		awk '{OFS="\t";print $1,$4,$6,$12,$14}'  > ${pp}_Fst_sweeps.tsv
done
