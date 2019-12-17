module load pal2nal
module load paml

prepare_codeml_files(){
	p1=$1
	p2=$2
	cat test.cnt | sed -E "s/test./${p1}_${p2}./" > codeml/${p1}_${p2}.cnt
	pal2nal.pl codeml/${p1}_${p2}.aln codeml/${p1}_${p2}_cds.fasta -output paml -nogap > codeml/${p1}_${p2}.codon
}


while read pair;do prepare_codeml_files $pair;done < single_copy_orthos.txt

cd codeml

for f in *.cnt;do
	codeml $f
done
