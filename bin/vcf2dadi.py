#!/usr/bin/env python3

import argparse
import sys
import vcf
import pdb
import math, statistics

parser = argparse.ArgumentParser(
	formatter_class=argparse.RawDescriptionHelpFormatter,
	description='Convert vcf to dadi format',
	epilog='''
Population list file assigns each sample to a population.
It should be whitespace separated, like this:

sample1    population1
sample2    population1
sample3    population2

The vcf-field argument should be used to select the method
for calculating allele frequencies. Use genotypes (GT) only for 
variants called with deep sequencing to avoid accumulation 
of heterozygote bias. The RO/AO option should be used for shallow 
sequenced samples as it will estimate allele frequencies directly
from numbers of observations in the AO and RO fields (ie reads).

''')

parser.add_argument('-p popfile',
	dest='popfile',
	metavar='FILE',
	required=True,
	type=argparse.FileType('r'),
	help='Population list (see below)')

parser.add_argument('-r ref',
	dest='ref_genome',
	metavar='FILE',
	type=argparse.FileType('r'),
	help='Reference genome (ingroup) in FASTA format')

parser.add_argument('-a ancestral',
	dest='anc_genome',
	metavar='FILE',
	type=argparse.FileType('r'),
	help='Reference genome (ancestral/outgroup) in FASTA format')

parser.add_argument('--ref-name name',
	dest='ref_name',
	default='ref',
	help='Name for reference species')

parser.add_argument('--anc-name name',
	dest='anc_name',
	default='anc',
	help='Name for ancestral or outgroup species')

parser.add_argument('--vcf-field',
	dest='vcf_field',
	choices=['ROAO','GT'],
	default='ROAO',
	help='The VCF field to use for calculating allele frequencies. (Default: ROAO)')

parser.add_argument('--max-dev value',
	dest='max_dev',
	type=float,
	default='0.5',
	help='Maximum value by which call ratios from populations may differ. Use this setting to filter out sites with biased call ratios perhaps indicating differences in mapping rate between populations')

parser.add_argument('--min-call value',
	dest='min_call',
	type=float,
	default='0.5',
	help='Minimum proportion of samples with genotype calls at each site')

parser.add_argument('--strict-pops',
	action='store_true',	
	dest='strict_pops',
	help='Exit if there are any samples in the vcf that are not in the popfile')

parser.add_argument('--sweepfinder',
	action='store_true',
	dest='sweepfinder',
	help='Write outputs in Sweepfinder2 format rather than dadi. Also allow non-snp alleles')



parser.add_argument('invcf',nargs='?',
	type=argparse.FileType('r'),
	help='VCF File',default=sys.stdin)

args  = parser.parse_args()

# Parse the popfile
#
rows = ( line.strip().split() for line in args.popfile )
sample2pop = { row[0]:row[1] for row in rows if len(row)>1}


vcf_reader = vcf.Reader(args.invcf)

populations = [ p for p in set(sample2pop.values()) ]

def popsize(popname,sample2pop):
	return(sum(1 if p==popname else 0 for p in sample2pop.values()))


population_sample_sizes = { pname:popsize(pname,sample2pop) for pname in populations }
# Write the header
# Example of the dadi format
# Human	Chimp Allele1 	YRI	CEU Allele2 YRI CEU Gene Position 
# ACG 	ATG 	C 		29 	24 	T 		1 	0 	abcb 	1289
# CCT 	CCT 	C 		29 	23 	G 		3 	2 	abcb 	1345

if args.sweepfinder:
	if len(population_sample_sizes)!=1:
		exit("Sweepfinder2 mode requires a single population")
	sys.stdout.write("position\tx\tn\tfolded\n")
else:
	sys.stdout.write(args.ref_name+'\t'+args.anc_name+'\t'+'Allele1'+'\t')
	for p in populations:
		sys.stdout.write(p+'\t')

	sys.stdout.write('Allele2'+'\t')

	for p in populations:
		sys.stdout.write(p+'\t')
	sys.stdout.write('Gene\tPosition\n')

num_rejected_max_dev=0
num_rejected_min_call=0

# This script was designed to work from freebayes outputs which have this info in the FORMAT string
#
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality, the Phred-scaled marginal (or unconditional) probability of the called genotype">
##FORMAT=<ID=GL,Number=G,Type=Float,Description="Genotype Likelihood, log10-scaled likelihoods of the data given the called genotype for each possible genotype generated from the reference and alternate alleles given the sample ploidy">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=DPR,Number=A,Type=Integer,Description="Number of observation for each allele">
##FORMAT=<ID=RO,Number=1,Type=Integer,Description="Reference allele observation count">
##FORMAT=<ID=QR,Number=1,Type=Integer,Description="Sum of quality of the reference observations">
##FORMAT=<ID=AO,Number=A,Type=Integer,Description="Alternate allele observation count">
##FORMAT=<ID=QA,Number=A,Type=Integer,Description="Sum of quality of the alternate observations">
##FORMAT=<ID=MIN,Number=1,Type=Integer,Description="Minimum depth in gVCF output block.">


for record in vcf_reader:

	# Our method for determining allele frequencies is based on the ROAO field
	# RO is reference observations
	# AO is alternate observations
	# In general RO and AO should add up to the read depth 
	# If they don't it means we have reads that don't match either allele
	# Which could be because the site is multi-allelic or because there are errors
	# We exclude these sites
	# Could this be a cause of spectrum distortions?
	#
	if args.vcf_field=='ROAO':

		# Ref,Alt observation counts for each population
		pop_obs = { p:[0,0] for p in populations }

		# Total number of genotypes involved in calculating pop_obs
		num_pop_calls = { p:0 for p in populations }
		num_unaccounted_depth = 0

		for call in record.samples:
			call_pop=sample2pop.get(call.sample,None)

			if call_pop==None:
				if args.strict_pops:
					exit("Sample, "+call.sample+" found in vcf does not exist in popfile")
			else:
				if call.data.AO!=None and call.data.RO!=None:
					if (call.data.RO + call.data.AO != call.data.DP):
						num_unaccounted_depth+=1
					else:
						pop_obs[call_pop][0]+=call.data.RO
						pop_obs[call_pop][1]+=call.data.AO
						num_pop_calls[call_pop]+=1

		# Check for invalid input data
		#
		if len(record.ALT)>1:
			exit('SNPs for should be biallelic but found multiple alternate alleles, '+str(record.ALT))

		if not args.sweepfinder:	
			if len(record.REF)>1:
				exit('Only SNPs allowed but found non-SNP ref allele, '+record.REF)

			if len(record.ALT[0])>1:
				exit('Only SNPs allowed but found non-SNP alt allele, '+str(record.ALT[0]))

		# Reject samples that do not meet QC filters
		#
		call_ratios = { p:num_pop_calls[p]/population_sample_sizes[p] for p in populations}

		if min(call_ratios.values()) < args.min_call:
			num_rejected_min_call+=1
			continue

		if max([ abs(r - statistics.mean(call_ratios.values())) for r in call_ratios.values()]) > args.max_dev:
			num_rejected_max_dev+=1
			continue

		pop_freqs = { p:[0,0] for p in populations}

		for p in populations:
			pop_depth = sum(pop_obs[p])
			pop_freqs[p][0] = int(round(2*num_pop_calls[p]*pop_obs[p][0]/pop_depth))		
			pop_freqs[p][1] = int(round(2*num_pop_calls[p]*pop_obs[p][1]/pop_depth))
			# Check arithmetic
			if sum(pop_freqs[p]) != num_pop_calls[p]*2:
				exit("Arithmetic error. "+str(pop_freqs[p])+" does not add up to "+str(num_pop_calls[p]))
		#
		# Write output
		#
		if args.sweepfinder:
			p = populations[0]
			sys.stdout.write(str(record.POS)+"\t"+str(min(pop_freqs[p]))+"\t"+str(sum(pop_freqs[p]))+"\t1\n")

			
		else:
			# Ingroup TODO: use ref genome
			sys.stdout.write("-"+record.REF+"-"+'\t')
			# Outgroup TODO: use ancestral genome
			sys.stdout.write("-"+record.REF+"-"+'\t')
			sys.stdout.write(record.REF+'\t')
			for p in populations:
				sys.stdout.write(str(pop_freqs[p][0])+'\t')


			sys.stdout.write(str(record.ALT[0])+'\t')
			
			for p in populations:
				sys.stdout.write(str(pop_freqs[p][1])+'\t')
			sys.stdout.write(record.CHROM+'\t')
			sys.stdout.write(str(record.POS)+'\n')

	elif args.vcf_field=='GT':
			exit("Using the GT field to calculate allele frequencies is not implemented yet")
	else:
		exit('Unrecognised vcf_field option')

if ( num_rejected_min_call>0):
	sys.stderr.write("Skipped "+str(num_rejected_min_call)+" sites with too few called samples\n")

if ( num_rejected_max_dev>0):
	sys.stderr.write("Skipped "+str(num_rejected_max_dev)+" sites with biased call ratios\n")


