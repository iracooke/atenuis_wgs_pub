module load vcftools
module load vcflib
module load parallel
module load bedtools2

#
# Min QUAL Score 30
# Biallelic
# Min depth (av/3) and Max depth (av*2)
# ---- Not using -----
# At least 3 individuals required for the alternate allele
#		--mac 2 \
#
do_filtering_1(){
	input=$1
	output=${input%.vcf}_Q30_biallelic_MD09_54_mingq20_mmiss0.5.vcf

	if [ ! -f $output ]; then
		vcftools --vcf $input --stdout --recode --recode-INFO-all \
		--minQ 30 \
		--min-alleles 2 --max-alleles 2 \
		--min-meanDP 0.9 --max-meanDP 5.4 \
		--minGQ 20 \
		--max-missing 0.5 \
		> $output
	fi

	echo $output
}

# Same basic filtering as filtering_1 but without removing low GQ 
# variants or sites with excessive numbers of missing variants
# This is because we can make use of the allele counts even for 
# very low coverage variants
#
# We can also incorporate a coverage bias filter between populations
# directly in our vcf2dadi script
#
#
do_filtering_1_dadi(){
	input=$1
	output=${input%.vcf}_Q30_biallelic_MD09_54.vcf

	if [ ! -f $output ]; then
		vcftools --vcf $input --stdout --recode --recode-INFO-all \
		--minQ 30 \
		--min-alleles 2 --max-alleles 2 \
		--min-meanDP 0.9 --max-meanDP 5.4 \
		> $output
	fi

	echo $output
}



# This will filter variants with bad AB values (ie not 0.5 for hets)
# We choose a threshold of 20 because it means our chance of falsely 
# removing a good variant is only 1%
do_filtering_ab(){
	input=$1
	output=${input%.vcf}_AB.vcf

	if [ ! -f $output ]; then
		vcffilter -s -f "ABP < 20" $input > $output
	fi

	echo $output
}

do_filtering_mdust(){
	input=$1
	output=${input%.vcf}_mdust.vcf

	if [ ! -f $output ]; then
		head -n 10000 $input | awk '/^#/{print $0}' > $output
		bedtools subtract -a $input -b aten_final_0.1.mdust.bed -A >> $output
	fi

	echo $output
}

do_filtering_snps(){
	input=$1
	output=${input%.vcf}_snps.vcf

	if [ ! -f $output ]; then
		vcftools --vcf $input --stdout --recode --recode-INFO-all \
		--remove-indels > $output
	fi

	echo $output
}

do_filtering_indels(){
	input=$1
	output=${input%.vcf}_indels.vcf

	if [ ! -f $output ]; then
		vcftools --vcf $input --stdout --recode --recode-INFO-all \
		--keep-only-indels > $output
	fi

	echo $output
}

do_filtering_snpgap(){
	input=$1
	indels=$2
	indelwindows=${indels%.vcf}_windows.bed
	output=${input%.vcf}_snpgap.vcf

	if [ ! -f $output ]; then
		# Also make a window file around indels
		cat $indels | awk -F '\t' '!/^#/{printf("%s\t%s\t%s\n",$1,$2-10,$2+10)}' > $indelwindows


		head -n 10000 $input | awk '/^#/{print $0}' > $output
		bedtools subtract -a $input -b $indelwindows -A >> $output
	fi

	echo $output
}

do_filtering_smlindels(){
	input=$1
	output=${input%.vcf}_smlindels.vcf

	if [ ! -f $output ]; then
		vcffilter -s -f "LEN < 50" $input > $output
	fi

	echo $output
}

do_combine_snps_indels(){
	output=$1
	unsorted_out=${output%_filtered.vcf}_unsorted_filtered.vcf

	if [ ! -f $output ]; then
		vcf-concat $2 $3 > ${unsorted_out}
		vcf-sort -t /fast/tmp ${unsorted_out} > $output
		rm ${unsorted_out}
	fi

	echo $output
}


do_filtering_lowcov_down(){
	contig=$1

	out1=$(do_filtering_1 vcf_by_contig_lowcov_down/${contig}.vcf)

	out2=$(do_filtering_ab $out1)

	out3=$(do_filtering_mdust $out2)

	snps=$(do_filtering_snps $out3)

	indels=$(do_filtering_indels $out3)

	snpgap=$(do_filtering_snpgap $snps $indels)

	smlindels=$(do_filtering_smlindels $indels)

	final_filtered=$(do_combine_snps_indels vcf_by_contig_lowcov_down/${contig}_filtered.vcf $snpgap $smlindels)
}

do_filtering_dadi(){
	contig=$1

	out1=$(do_filtering_1_dadi vcf_by_contig_dadi/${contig}.vcf)

	out2=$(do_filtering_ab $out1)

	out3=$(do_filtering_mdust $out2)

	snps=$(do_filtering_snps $out3)

	indels=$(do_filtering_indels $out3)

	snpgap=$(do_filtering_snpgap $snps $indels)

	smlindels=$(do_filtering_smlindels $indels)

	final_filtered=$(do_combine_snps_indels vcf_by_contig_dadi/${contig}_filtered.vcf $snpgap $smlindels)

}


# functions for filtering steps
#

export -f do_filtering_1 do_filtering_ab do_filtering_mdust do_filtering_snps do_filtering_indels do_filtering_snpgap do_filtering_smlindels do_combine_snps_indels

export -f do_filtering_dadi do_filtering_1_dadi

parallel -j 40 do_filtering_dadi ::: $(cat contig_list.txt | tr '\n' ' ')

export -f do_filtering_lowcov_down

parallel -j 40 do_filtering_lowcov_down ::: $(cat contig_list.txt | tr '\n' ' ')


