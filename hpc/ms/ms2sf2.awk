function min(a,b){
	if ( a > b){ return b}
	return a
}

BEGIN { 
	start=0 
}

$1~/^ms/ {
	nsam=$2
	nsites=$8
	next
}

$1~/^positions/ {
	for (i=2;i<=NF;i++){
		sites[i]=$i
		AFS[i]=0
	}
	start=1
	next
}

{
	if ( start == 1){
		FS=""
		for (i=2;i<=NF;i++){
			AFS[i]+=$i
		}		
	}
}

END {
	printf("position\tx\tn\tfolded\n")
	for (i in sites){
		maf = min(AFS[i],nsam-AFS[i])
		printf("%d\t%d\t%d\t1\n",sites[i]*nsites, maf, nsam)
	}
}