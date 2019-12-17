conda activate dadi
module load parallel

run_model(){
	model=$2
	round=$1
	echo $model $round
	python dd_${model}.py $round
}

export -f run_model

parallel --lb -j 12 run_model\
 ::: $(seq 1 10) \
 ::: 'no_mig' 'sym_mig' 'asym_mig' 'priorsize_asym_mig' 'asym_mig_size' 'isolation_asym_mig'

	  
