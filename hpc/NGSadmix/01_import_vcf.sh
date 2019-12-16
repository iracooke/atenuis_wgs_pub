module load parallel
module load angsd

import_vcf(){
	prefix=$1
	f=$2
	bn=$(basename $f)
	contig=${bn%_filtered.vcf}

  # Here aten_final_0.1.fasta.fai is the genome index. Outputs are in Beagle format for use with PCAngst and NGSAdmix
  #
	angsd -vcf-gl $f -out ${prefix}_freebayes_bgl/${contig} -fai aten_final_0.1.fasta.fai -nind 148 -doMaf 1 -doGlf 2 -doMajorMinor 1 -SNP_pval 1e-6	
}


export -f import_vcf

mkdir -p accelerate_freebayes_bgl

parallel -j 40 import_vcf 'accelerate' ::: $(ls ../freebayes_qc/vcf_by_contig_lowcov_down/*_filtered.vcf | tr '\n' ' ')


# Rather convoluted way of concatenating the beagle files
#
for prefix in 'accelerate';do
	cat ${prefix}_freebayes_bgl/*.beagle.gz > ${prefix}_freebayes_bgl.beagle.gz
	gunzip -c ${prefix}_freebayes_bgl.beagle.gz | head -n 1 >> ${prefix}_freebayes_bgl_clean.beagle
	gunzip -c ${prefix}_freebayes_bgl.beagle.gz | grep -v '^marker' >> ${prefix}_freebayes_bgl_clean.beagle
	gzip ${prefix}_freebayes_bgl_clean.beagle
	mv ${prefix}_freebayes_bgl_clean.beagle.gz ${prefix}_freebayes_bgl.beagle.gz

	# We should extract our list of individuals from the vcf file to be sure it is in the correct order
	#
	grep 'CHROM' ../freebayes_qc/vcf_by_contig_lowcov_down/Sc0000000_filtered.vcf | awk '{for(i=10;i<=NF;i++){printf("%s\n",$i)}}' > ${prefix}_individuals.txt
done


