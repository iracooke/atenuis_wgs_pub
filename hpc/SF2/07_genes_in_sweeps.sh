

# These bedtools commands should report 

get_sweep_genes(){
	POP=$1
	SCORE=$2
	bedtools window -b <(cat ${POP}_${SCORE}_sweeps.gff| sed s/_mixed//) \
					-a aten_0.11.maker_post_001.genes.gff | grep 'gene' \
					> ${POP}_${SCORE}_sweep_genes.gff

#					-a aten_0.11.maker_post_001.genes.gff -c | awk '$10>0{print}' | grep 'gene' \
}



for LR in 10 20 30 50 100;do
	for pop in nomi mi;do
		echo ${LR} ${pop}
		get_sweep_genes $pop $LR
	done
done

