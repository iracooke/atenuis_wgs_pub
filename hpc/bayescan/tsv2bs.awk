BEGIN{
	pop=""
	snp_count=-1
	pop_count=0
	printf("[loci]=27109\n\n[populations]=4\n")
}

$1>0{
	if( pop != $3){
		pop=$3
		snp_count=1
		pop_count+=1
		printf("\n[pop]=%s\n",pop_count)
	} else {
		snp_count+=1
	}
	print(snp_count,$1,"2",$2,$1-$2)
	if ( pop_count==1 ){
		print(snp_count,$4) > "snp2id.txt"
	}
}
