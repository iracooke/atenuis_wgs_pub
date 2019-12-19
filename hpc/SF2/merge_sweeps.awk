function min(a,b){
	if ( a > b){ return b}
	return a
}

function max(a,b){
	if ( a < b){ return b}
	return a
}

 BEGIN{OFS="\t"}

{
	print $1,min($4,$13),max($5,$14),max($6,$16)
}