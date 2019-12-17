BEGIN {
	max_likelihood = -1000000000
}

{
	if ( $3 > max_likelihood){
		max_likelihood=$3
		best_model=$0
	}
}

END {
	print best_model
}