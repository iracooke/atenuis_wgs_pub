#!/usr/bin/env python3

# takes a PSMC' (msmc) output file as input and outputs an input file for ms

import argparse, math

parser = argparse.ArgumentParser()
parser.add_argument("msmc", help="path/prefix for MSMC files to take as input")
parser.add_argument("--form",help="desired format of ms output", type=str,choices=['trees','snps','both'], default='both') 
parser.add_argument("--chromL", help="Length of chromosome to simulate",type=int,default=1000000)
args = parser.parse_args()


MSNeT = []

with open(args.msmc,'r') as infile:
	# skip the header
	for line in infile:
		if line[0] != 't': #skip the header line
			MSMCparams = [float(x) for x in line.split()]
			# MSMCparams should be [index, mu t_left, mu t_right, 1/(2Ne(t) mu)]
			# ms needs [t_left/(4N0), Ne/N0]
			if MSMCparams[0]==0: #the first line of data; only use it to normalize the rest
				norm = MSMCparams[-1] # = 1/(2N0 mu)
			else:
				if MSNeT == [] or norm / MSMCparams[-1] != MSNeT[-1][-1]: #there has been a pop size change from the previous time step
					MSNeT.append([MSMCparams[1] * norm / 2.0, norm / MSMCparams[-1]])

# find the estimated mutation and recombination rates so that we can take the ratio:			
# We use mutation rate obtained as the average of all bootstraps
# grep 'mutationRate' *.log | awk '{av+= $2}END{print av/NR}'
#
# And similarly for recombination rate
# tail -q -n 1 *.loop.txt | awk '{av+= $1}END{print av/NR}'

m0 = 0.00292
r = 0.000881625

#print("Assuming mutation rate "+str(m0)+" and recombination rate "+str(r))
			
if args.form in ('trees','both'):
	print('-T', end=' ')			
if args.form in ('snps','both'):
	print('-t '+str(2.0* args.chromL/norm), end=' ')
print('-r '+str(2.0*args.chromL/norm * r/m0)+ ' ' + str(args.chromL) + ' -p ' + str( math.ceil( math.log10(args.chromL) ) ) + ' -eN ' + ' -eN '.join(' '.join(str(x) for x in timestep) for timestep in MSNeT))


