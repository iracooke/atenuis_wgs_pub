generate_bootstrap_data(){
	../sw/msmc-tools/multihetsep_bootstrap.py -n 100 -s 500000 --chunks_per_chromosome 20 --nr_chromosomes 20 FI_bootstrap workfiles/*FI*.multihetsep.txt
	../sw/msmc-tools/multihetsep_bootstrap.py -n 100 -s 500000 --chunks_per_chromosome 20 --nr_chromosomes 20 MI_bootstrap workfiles/*MI*.multihetsep.txt	
}


run_replicate(){
	bsn=$1

	msmc2 -t 10 -I 0,1 -o ${bsn}_withinFI_msmc FI_bootstrap_${bsn}/*.txt
	msmc2 -t 10 -I 0,1 -o ${bsn}_withinMI_msmc MI_bootstrap_${bsn}/*.txt

}

module load parallel
module load msmc

export -f run_replicate generate_bootstrap_data

generate_bootstrap_data

parallel -j 20 run_replicate ::: $(seq -s ' ' 100)

mkdir -p bootstrap_results

mv *_msmc.final.txt bootstrap_results/
mv *_msmc.txt bootstrap_results/
