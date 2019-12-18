
# This script implements section 4 of the SweepFinder2 manual

module load SF2
module load parallel

# Concatenate scaffold-wise spectra into a combined file

for pop in 'nomi' 'pi' 'pr' 'mi' 'fi' 'di';do
	if [[ ! -f ${pop}.combinedfreqfile ]]; then
		echo "Combining freqs for ${pop}"
		echo "position	x	n	folded" > ${pop}.combinedfreqfile
		for f in ${pop}/*.af;do
			tail -n +2 $f >> ${pop}.combinedfreqfile
		done
	fi
done

generate_empirical_spectrum(){
	pop=$1
	SweepFinder2 -f ${pop}.combinedfreqfile ${pop}.spectfile
}

# Run SweepFinder to generate empirical frequency spectrum across the genome

export -f generate_empirical_spectrum

parallel -j 6 generate_empirical_spectrum ::: 'pi' 'pr' 'fi' 'di' 'nomi' 'mi'

