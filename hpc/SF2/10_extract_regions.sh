module load htslib
module load bcftools


# bcftools view -R nomi_100_sweep_regions.txt nomi.vcf.gz -o nomi_100_sweep_regions.vcf
# bcftools sort nomi_100_sweep_regions.vcf -o nomi_100_sweep_regions.sorted.vcf
# bgzip nomi_100_sweep_regions.sorted.vcf
# tabix nomi_100_sweep_regions.sorted.vcf.gz



bcftools view -R nomi_100_sweep_regions.txt mi.vcf.gz -o mi_100_sweep_regions.vcf
bcftools sort mi_100_sweep_regions.vcf -o mi_100_sweep_regions.sorted.vcf
bgzip mi_100_sweep_regions.sorted.vcf
tabix mi_100_sweep_regions.sorted.vcf.gz
