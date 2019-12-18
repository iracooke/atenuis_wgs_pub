# 
# Our goal here is to simulate large enough genomic regions to capture typical sweeps. We use 1Mb as our region size
#
# We need to incorporate recombination and assume a rate of 1e-8 per base
#
# Mutation rate is set to 2.9e-8 which was obtained from Mao et al
#
# We generate samples of 30 individuals which corresponds to our population sizes
#
# Ne is set to 50000 in accordance with estimates from MSMC
#

mkdir -p no_demography

# Theta = 4NµL = 4 * 50000 * 2.9e-8 * 1e6 = 5800

# Rho = 4NrL = 4 * 50000 * 1e-8 * 1e6 = 2000

simulate_no_demography(){
	sample_num=$1
	ms 60 1 -t 5800 -r 2000 1000000 -p 9 | \
	awk -f ms2sf2.awk | \
	sort -n -k 1 > no_demography/sample_${sample_num}.af
}

export -f simulate_no_demography

parallel -j 40 simulate_no_demography ::: $(seq 1 500)


