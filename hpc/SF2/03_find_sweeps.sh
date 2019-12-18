module load SF2/1.0
module load parallel

find_sweeps(){
	pop=$1
	contig=$2

	mkdir -p gridfiles

	if [ ! -e  gridfiles/${contig}.grid.txt ]; then
		make_gridfile ${contig} > gridfiles/${contig}.grid.txt
	fi

	SweepFinder2 -lu gridfiles/${contig}.grid.txt ${pop}/${contig}_filtered.vcf.af ${pop}.spectfile ${pop}/${contig}.sf2

}

make_gridfile(){
	contig=${1}
	min_pos=1000000000
	max_pos=0
	for pop in 'nomi' 'pi' 'pr' 'mi' 'fi' 'di';do
		mp=$(head -n 2 ${pop}/${contig}_filtered.vcf.af | grep -v 'position' | awk '{print $1}')
		if ((mp < min_pos));then
			min_pos=$mp
		fi

		mp=$(tail -n 1 ${pop}/${contig}_filtered.vcf.af | awk '{print $1}')
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


export -f find_sweeps make_gridfile

for pop in 'nomi' 'pi' 'pr' 'mi' 'fi' 'di';do
	parallel -j 40 find_sweeps ${pop} ::: $(cat contig_list.txt | tr '\n' ' ')
done
