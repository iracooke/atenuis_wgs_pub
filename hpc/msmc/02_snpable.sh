module load bioawk
module load seqbility
module load parallel
module load bwa

genome=aten_final_0.1.fasta

splitfa ${genome} 100 | split -l 20000000

do_bwa(){
	in=$1
	bwa aln -t 8 -R 1000000 -O 3 -E 3 aten_final_0.1.fasta ${in} > ${in}.sai 	
}

do_sampe(){
	input=${1%.sai}
	bwa samse aten_final_0.1.fasta ${input}.sai ${input} | gzip -c > ${input}.sam.gz
}

export -f do_bwa do_sampe


parallel -j 6 do_bwa ::: $(ls x* | tr '\n' ' ')

parallel -j 40 do_sampe ::: $(ls *.sai | tr '\n' ' ')

gzip -dc x??.sam.gz | gen_raw_mask.pl > rawMask_100.fa

gen_mask -l 100 -r 0.5 rawMask_100.fa > mask_100_50.fa

# Create a cut-down version with only 1M or longer contigs

cat mask_100_50.fa | bioawk -c fastx 'length($seq)>1000000{printf(">%s\n%s\n",$name,$seq)}' > mask_100_50_1M.fa

grep '>' mask_100_50_1M.fa | sed 's/>//' > contig_list_1M.txt
