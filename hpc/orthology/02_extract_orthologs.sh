module load mafft
module load samtools

# Extract single copy orthos in adi / aten pairs
#
cat acropora.proteinortho | awk '($1 == 3) && ($2==3) && ($3==1){print $4,$6}' > single_copy_orthos.txt

extract_pair(){
	p1=$2
	p2=$3
	dest=$1
	cat <(samtools faidx adi.faa $p1) <(samtools faidx aten.faa $p2) | sed s/X// > $dest/${p1}_${p2}.fasta
	cat <(samtools faidx adi_cds.fna $p1) <(samtools faidx aten_cds.fna $p2) > $dest/${p1}_${p2}_cds.fasta
	echo "($p1,$p2);" > $dest/${p1}_${p2}.tree
	mafft --auto $dest/${p1}_${p2}.fasta > $dest/${p1}_${p2}.aln
}


export -f extract_pair

while read pair;do extract_pair "codeml" $pair;done < single_copy_orthos.txt

