source activate popgen
module load parallel

generate_sf2_inputs_for_pop(){
	pop=$1
	f=$2
	mkdir -p $pop
	../../bin/vcf2dadi.py --sweepfinder -p clean.poplist.${pop}.txt $f | awk '$2>0' > ${pop}/$(basename $f).af
}


export -f generate_sf2_inputs_for_pop

for pop in 'fi' 'pi' 'di' 'pr' 'mi' 'nomi';do
 	parallel -j 40 generate_sf2_inputs_for_pop $pop ::: $(ls -1 ../freebayes_qc/vcf_by_contig_lowcov_down/*_filtered.vcf | tr '\n' ' ')
done
