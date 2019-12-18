module load parallel
module load bedtools2

get_sweep_regions(){
	LR=$1
	pop=$2
	f=$3
	bn=$(basename $f); 
	../../bin/sf2gff.py -t ${LR} ${bn%.sf2} $f >> ${pop}_${LR}_sweeps.gff
}

export -f get_sweep_regions

LR=10
for pop in 'nomi' 'pi' 'pr' 'mi' 'fi' 'di';do	
	echo $pop
	rm ${pop}_${LR}_sweeps.gff
	parallel -j 40 get_sweep_regions ${LR} ${pop} ::: $(ls ${pop}/*.sf2 | tr '\n' ' ')
done
