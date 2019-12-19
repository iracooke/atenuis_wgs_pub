
if [[ ! -f aten_final_0.1.mdust.txt ]]; then
	module load mdust
	mdust aten_final_0.1.fasta -c > aten_final_0.1.mdust.txt
	awk '{print $1,$3,$4}' aten_final_0.1.mdust.txt > aten_final_0.1.mdust.bed
fi



