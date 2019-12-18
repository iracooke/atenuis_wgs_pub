module load angsd/0.928
module load parallel

# Since we don't know the ancestral state we use the ref as the ancestral for Fst calculations
# See here https://github.com/ANGSD/angsd/issues/65
#
do_pop(){
	POP=$1

	REF=aten_final_0.1.fasta

    angsd -nThreads 32  -b ${POP}_bamlist.txt -anc $REF -ref $REF -out ${POP}_af \
     -minMapQ 5 -minQ 20 -GL 1 -doSaf 1 \
     -sites good_sites.txt

}

export -f do_pop

if [[ ! -f good_sites.txt.idx ]]; then
	angsd sites index good_sites.txt
fi

parallel -j 5 do_pop ::: DI FI PI PR MI

do_pop north


