for pop in DI FI MI PR PI;do

   samtools merge ${pop}_sweep_regions.bam $(find ./ -name "${pop}*.bam" | tr '\n' ' ')
	samtools index ${pop}_sweep_regions.bam

done
