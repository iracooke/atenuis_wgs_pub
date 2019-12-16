#!/usr/bin/bash

genome=aten_final_0.1.fasta

for set in accelerate;do

	python /usr/local/sci-irc/sw/pcangsd/pcangsd.py -beagle ${set}_freebayes_bgl.beagle.gz \
		-threads 40 \
		-admix \
		-admix_save \
		-admix_auto 10000 \
		-o ${set}_fb_pcangst_2

done
