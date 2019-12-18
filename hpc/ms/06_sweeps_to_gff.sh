module load parallel
module load bedtools2

get_sweep_regions(){
	LR=$1
	model=$2
	pop=$3
	f=$4
	bn=$(basename $f); 
	../../bin/sf2gff.py -t ${LR} ${bn%.sf2} $f > ${f%.sf2}_sweeps.gff 
}

export -f get_sweep_regions

LR=10


for model in $(cat ../dadi/best_models.tsv | sort -k 3 -n -r | awk '{printf("ms_%s ",$1)}');do
	num_sf2=$(ls ${model}/*_a.sf2 | wc -l)
	if [[ ${num_sf2} -gt 400 ]]; then
		echo $model
		rm ${model}_mi_${LR}_sweeps.gff
		rm ${model}_nomi_${LR}_sweeps.gff	
		parallel -j 40 get_sweep_regions ${LR} ${model} 'mi' ::: $(ls ${model}/*_a.sf2 | tr '\n' ' ')
		parallel -j 40 get_sweep_regions ${LR} ${model} 'nomi' ::: $(ls ${model}/*_b.sf2 | tr '\n' ' ')
		cat ${model}/*_a_sweeps.gff > ${model}_mi_${LR}_sweeps.gff
		cat ${model}/*_b_sweeps.gff > ${model}_nomi_${LR}_sweeps.gff
	else
		echo "Skipping $model $num_sf2"
	fi
done


parallel -j 40 get_sweep_regions ${LR} ms_msmc_MI 'mi' ::: $(ls ${model}/*.sf2 | tr '\n' ' ')
parallel -j 40 get_sweep_regions ${LR} ms_msmc_FI 'nomi' ::: $(ls ${model}/*.sf2 | tr '\n' ' ')
cat ms_msmc_MI/*_sweeps.gff > ms_msmc_mi_${LR}_sweeps.gff
cat ms_msmc_FI/*_sweeps.gff > ms_msmc_nomi_${LR}_sweeps.gff
