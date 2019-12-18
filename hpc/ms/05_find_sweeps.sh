module load SF2/1.0
module load parallel

find_sweeps(){
	pop=$1
	contig=$2

	mkdir -p gridfiles

	if [ ! -e  gridfiles/${contig}.grid.txt ]; then
		make_gridfile ${contig} > gridfiles/${contig}.grid.txt
	fi

	SweepFinder2 -lu gridfiles/${contig}.grid.txt ${pop}/${contig}.af ${pop}.spectfile ${pop}/${contig}.sf2

}


find_sweeps_twopop(){
	model=$1
	contig=$2

	mkdir -p gridfiles

	if [ ! -e  gridfiles/${contig}.grid.txt ]; then
		make_gridfile ${contig} > gridfiles/${contig}.grid.txt
	fi

	if [[ ! -e ${model}/${contig}_a.sf2 ]]; then
		SweepFinder2 -lu gridfiles/${contig}.grid.txt ${model}/${contig}_a.af ${model}_a.spectfile ${model}/${contig}_a.sf2
		SweepFinder2 -lu gridfiles/${contig}.grid.txt ${model}/${contig}_b.af ${model}_b.spectfile ${model}/${contig}_b.sf2
	fi

}

make_gridfile(){
	contig=${1}
	min_pos=1000000000
	max_pos=0
	for pop in 'ms_msmc_FI';do
		mp=$(head -n 2 ${pop}/${contig}.af | grep -v 'position' | awk '{print $1}')
		if ((mp < min_pos));then
			min_pos=$mp
		fi

		mp=$(tail -n 1 ${pop}/${contig}.af | awk '{print $1}')
		if ((mp > max_pos));then
			max_pos=$mp
		fi
	done

	grid_pos=$((min_pos+1000))
	while ((grid_pos < max_pos)); do
		echo ${grid_pos}
		grid_pos=$((grid_pos+1000))
	done

}


export -f find_sweeps find_sweeps_twopop make_gridfile

for pop in 'ms_msmc_FI' 'ms_msmc_MI';do
	parallel -j 40 find_sweeps ${pop} ::: $(for s in $(seq 1 500);do echo "sample_$s";done | tr '\n' ' ')
done

for model in $(cat ../dadi/best_models.tsv | sort -k 3 -n -r | awk '{printf("ms_%s ",$1)}');do
	if [[ -e ${model}_a.spectfile ]]; then
		parallel -j 40 find_sweeps_twopop ${model} ::: $(for s in $(seq 1 500);do echo "sample_$s";done | tr '\n' ' ')
	fi
done
