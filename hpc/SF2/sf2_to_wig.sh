#!/bin/bash

pop=$1
report=$2

for f in ${pop}/*.sf2;do
	scaffold=$(basename $f)
	scaffold=${scaffold%.sf2}
	echo "variableStep chrom=${scaffold}"
	if [[ $report == "LR" ]]; then
		cat $f | grep -v 'location' | awk '{printf("%i\t%s\n",$1,$2)}'
	elif [[ $report == "alpha" ]]; then
		cat $f | grep -v 'location' | awk '{printf("%i\t%s\n",$1,$3)}'
	fi
done
