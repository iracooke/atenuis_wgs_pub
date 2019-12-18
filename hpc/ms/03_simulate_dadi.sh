module load ms
module load parallel

simulate_model(){
	model=$1
	sample_num=$2

	m=$(cat ../dadi/best_model_mscommands.tsv | awk -F '\t' -v model=$model '$1==model{print $0}')

	model=$(printf "$m\n" | awk -F '\t' '{print $1}')
	command=$(printf "$m\n" | awk -F '\t' '{print $2}')



	if [[ ! -f ms_${model}/sample_${sample_num}.ms ]]; then
		printf "Running $model : $command : ${sample_num}\n"
		mkdir -p ms_${model}
		ms $command -seeds $RANDOM $RANDOM $RANDOM > ms_${model}/sample_${sample_num}.ms
		cat ms_${model}/sample_${sample_num}.ms | awk -f ms2sf2_twopop.awk -v outfile=ms_${model}/sample_${sample_num}
		cat ms_${model}/sample_${sample_num}_a.txt | sort -n -k 1 > ms_${model}/sample_${sample_num}_a.af
		cat ms_${model}/sample_${sample_num}_b.txt | sort -n -k 1 > ms_${model}/sample_${sample_num}_b.af
	else
		printf "Skipping $model : ${sample_num}\n"
	fi

}

export -f simulate_model

models=$(cat ../dadi/best_models.tsv | sort -k 3 -n -r | awk '{printf("%s ",$1)}')

parallel -j 40 simulate_model \
 	::: $models \
 	::: $(seq 1 500)
