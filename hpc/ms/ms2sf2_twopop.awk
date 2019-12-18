function min(a,b){
	if ( a > b){ return b}
	return a
}

BEGIN { 
	start=0
	if (outfile==""){
		printf "You must provide a filename for output\n"
		exit
	}
}

$1~/^ms/ {
	nsam_a=$8
	nsam_b=$9
	nsites=$12
	next
}

$1~/^positions/ {
	for (i=2;i<=NF;i++){
		sites[i-1]=$i
		AFS_a[i-1]=0
		AFS_b[i-1]=0
	}
	start=NR
	FS=""
	next
}

{
	if ( start > 0){

		if (NR <= (start + nsam_a))
		{
			for (i=1;i<=NF;i++)
			{
				AFS_a[i]+=$i
			} 
		} else {
			for (i=1;i<=NF;i++)
			{
				AFS_b[i]+=$i
			}			
		}
	}
}
END {
	if (outfile!=""){
		outfile_a=sprintf("%s_a.txt",outfile)
		outfile_b=sprintf("%s_b.txt",outfile)

		printf("Sending sfs for population a to %s\n",outfile_a)

		printf("position\tx\tn\tfolded\n") > outfile_a
		for (i in sites){
			maf = min(AFS_a[i],nsam_a-AFS_a[i])
			printf("%d\t%d\t%d\t1\n",sites[i]*nsites, maf, nsam_a) > outfile_a
		}

		printf("Sending sfs for population b to %s\n",outfile_b)

		printf("position\tx\tn\tfolded\n") > outfile_b
		for (i in sites){
			maf = min(AFS_b[i],nsam_b-AFS_b[i])
			printf("%d\t%d\t%d\t1\n",sites[i]*nsites, maf, nsam_b) > outfile_b
		}

	}
}