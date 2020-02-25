
# Annotation table
echo raw_data/annotation_table.tsv > data.list

# NGSADdmix
echo hpc/NGSadmix/accelerate_fb_pcangst_2.K2.a625.0.qopt >> data.list
echo hpc/NGSadmix/accelerate_fb_pcangst_2.cov >> data.list

# msmc
echo hpc/msmc/bootstrap_results >> data.list

# dadi
echo hpc/dadi/optim >> data.list
echo hpc/dadi/dadi.thin1k.txt >> data.list

#OrthoFinder results
echo raw_data/Orthogroups.GeneCount.csv >> data.list

# ms
ls hpc/ms/*.gff >> data.list

# SF2
ls hpc/SF2/*.gff >> data.list

# angsd
ls hpc/angsd/*.txt >> data.list
ls hpc/angsd/*.pestPG >> data.list
ls hpc/angsd/*.tsv >> data.list

# symbiodinium_profiles
echo hpc/symbiodinium_profiles/genome_kraken_mpa >> data.list

# repeats
echo hpc/repeats/ >> data.list


# symbiont_mito
echo hpc/symbiodinium/AllSymbC1MitoConsensus_GoodCoverage.fasta   >> data.list   
echo hpc/symbiodinium/AllSymbC1MitoConsensus_GoodCoverageTrim.fasta >> data.list

tar -zcvf data.tgz -T data.list

rm data.list