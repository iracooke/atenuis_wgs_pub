# Prepare gff files for proteinortho

```
cat adi_v2/GCF_000222465.1_Adig_1.1_genomic.gff | grep 'CDS' | sed 's/protein_id/ID/' > adi.gff
cat amil/amil_1.0.maker_006.genes.gff | grep 'CDS' | sed 's/Parent/Name/' > amil.gff
cat aten/aten_0.11.maker_post_001.genes.gff | grep 'CDS' | sed 's/Parent/Name/' > aten.gff
```

# Prepare cds files for pal2nal

```
cat adi_v2/GCF_000222465.1_Adig_1.1_cds_from_genomic.fna | sed -E 's/.*_cds_([A-Z]+_[^_]+).*/>\1/' > adi_cds.fna
cat amil/amil_1.0.maker_006.cds.fasta > amil_cds.fna
```

# Extract synonymous substitution rates from codeml results 

```
grep 'dN/dS=' codeml/*.codeml | sed 's/:/ /' | awk 'BEGIN{OFS="\t";printf("genes\tt\tS\tN\tdNdS\tdN\tdS\n")}{print $1,$3,$5,$7,$9,$12,$15}'| sed 's/codeml\///' | sed s/.codeml// > adi_aten_codeml.txt 
```