module load angsd/0.928
module load parallel

# Since we don't know the ancestral state we use the ref as the ancestral for Fst calculations
# See here https://github.com/ANGSD/angsd/issues/65
#

#The .pestPG file is a 14 column file (tab separated). The first column contains information about the region. The second and third column is the reference name and the center of the window.

#We then have 5 different estimators of theta, these are: Watterson, pairwise, FuLi, fayH, L. And we have 5 different neutrality test statistics: Tajima's D, Fu&Li F's, Fu&Li's D, Fay's H, Zeng's E. The final column is the effetive number of sites with data in the window.

do_pop(){
	POP=$1

	REF=aten_final_0.1.fasta

    angsd -nThreads 32  -b ${POP}_bamlist.txt -anc $REF -out ${POP}_folded_af \
     -minMapQ 5 -minQ 20 -GL 1 -doSaf 1 -fold 1 -sites good_sites.txt

	realSFS ${POP}_folded_af.saf.idx -P 24 > ${POP}_folded_af.sfs

	angsd -nThreads 32  -b ${POP}_bamlist.txt -out ${POP}_folded_af -doThetas 1 -doSaf 1 -pest ${POP}_folded_af.sfs \
		-minMapQ 5 -minQ 20 -anc $REF -GL 1 -fold 1 -sites good_sites.txt

	thetaStat do_stat ${POP}_folded_af.thetas.idx

	thetaStat do_stat ${POP}_folded_af.thetas.idx -win 50000 -step 10000  -outnames ${POP}_theta.thetasWindow.gz

}

export -f do_pop

parallel -j 5 do_pop ::: DI FI PI PR MI

do_pop north

