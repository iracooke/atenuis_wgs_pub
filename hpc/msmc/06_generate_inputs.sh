module load parallel
module load samtools

generate_multihetsep(){
	c=$1
	ds=$2

	../sw/msmc-tools/generate_multihetsep.py  --mask=workfiles/${ds}_${c}_mask.bed.gz \
                          --mask=workfiles/chr${c}.mask.bed \
                          workfiles/${ds}_${c}.vcf.gz > workfiles/${c}_${ds}.multihetsep.txt
}

export -f generate_multihetsep

parallel -j 40 generate_multihetsep ::: $(cat contig_list_1M.txt | tr '\n' ' ') ::: 'FI-1-3' 'MI-1-4'
