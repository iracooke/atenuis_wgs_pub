module load ms
module load parallel

simulate_msmc(){
	pop=$1
	sample_num=$2

	command=$(./msmc2ms.py msmc_av_$pop.txt --form snps)

	mkdir -p ms_msmc_${pop}

	echo "ms 60 1 $command"

	ms 60 1 $command | \
	awk -f ms2sf2.awk | \
	sort -n -k 1 > ms_msmc_${pop}/sample_${sample_num}.af
}

export -f simulate_msmc

parallel -j 40 simulate_msmc 'MI' ::: $(seq 1 500)

parallel -j 40 simulate_msmc 'FI' ::: $(seq 1 500)

