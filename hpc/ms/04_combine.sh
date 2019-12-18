module load SF2
module load parallel

# Concatenate scaffold-wise spectra into a combined file

# For msmc (single populations)
for pop in 'ms_msmc_FI' 'ms_msmc_MI';do
	if [[ ! -f ${pop}.combinedfreqfile ]]; then
		echo "Combining freqs for ${pop}"
		echo "position	x	n	folded" > ${pop}.combinedfreqfile
		for f in ${pop}/*.af;do
			tail -n +2 $f >> ${pop}.combinedfreqfile
		done
	fi
done

# And for two pop simulations
for model in $(cat ../dadi/best_model_mscommands.tsv | awk '{printf("ms_%s ",$1)}');do
	if [[ -d ${model} ]]; then
		num_afs=$(ls ${model}/*_a.af | wc -l)
		if [[ ! -f ${model}_a.combinedfreqfile ]] && [[ ${num_afs} -eq 500 ]]; then
			echo "Combining freqs for ${model}"
			echo "position	x	n	folded" > ${model}_a.combinedfreqfile
			for f in ${model}/*_a.af;do
				tail -n +2 $f >> ${model}_a.combinedfreqfile
			done

			echo "position	x	n	folded" > ${model}_b.combinedfreqfile
			for f in ${model}/*_b.af;do
				tail -n +2 $f >> ${model}_b.combinedfreqfile
			done
		fi
	fi
done


generate_empirical_spectrum(){
	model=$1
	if [[ ! -f ${model}.spectfile ]]; then
		SweepFinder2 -f ${model}.combinedfreqfile ${model}.spectfile
	fi
}

generate_empirical_spectrum_twopop(){
	model=$1
	if [[ -d ${model} ]]; then
		if [[ ! -f ${model}_a.spectfile ]]; then
			SweepFinder2 -f ${model}_a.combinedfreqfile ${model}_a.spectfile
		fi
		if [[ ! -f ${model}_b.spectfile ]]; then
			SweepFinder2 -f ${model}_b.combinedfreqfile ${model}_b.spectfile
		fi
	fi
}


export -f generate_empirical_spectrum generate_empirical_spectrum_twopop


parallel -j 2 generate_empirical_spectrum ::: 'ms_msmc_FI' 'ms_msmc_MI'

parallel -j 6 generate_empirical_spectrum_twopop ::: $(cat ../dadi/best_model_mscommands.tsv | awk '{printf("ms_%s ",$1)}')
