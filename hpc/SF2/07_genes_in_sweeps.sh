get_sweep_genes(){
	POP=$1
	SCORE=$2
	bedtools window -b <(cat ${POP}_${SCORE}_sweeps.gff| sed s/_mixed//) \
					-a aten_0.11.maker_post_001.genes.gff | grep 'gene' \
					> ${POP}_${SCORE}_sweep_genes.gff
}

LR=10

for pop in nomi;do
	echo ${LR} ${pop}
	get_sweep_genes $pop $LR
done

