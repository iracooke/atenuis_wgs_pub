import sys
import os
import numpy
import dadi
import pylab
import dadi.Numerics
import dd_models
import csv

import argparse

parser = argparse.ArgumentParser(
	formatter_class=argparse.RawDescriptionHelpFormatter,
	description='Print ms commands and plottable demographies for dadi models')

parser.add_argument('--mode', choices=['plot', 'ms','params'], default='ms')

args  = parser.parse_args()

best_models={}
with open('best_models.tsv') as tsvfile:
	reader = csv.reader(tsvfile, delimiter='\t')
	best_models = { row[0]:[row[5],row[6]] for row in reader}

mu = 1.86e-8
rho = 1e-8
L = 1e6 # Length of region we want to simulate
P = 6

# Effective L for dadi using 15k of 14million sites is
Leff = 450e6*(15e3/14e6)

# We want to simulate blocks of length L = 1e5 so we need to rescale theta

if args.mode == 'plot':
	sys.stdout.write('model'+'\t'+'T'+'\t'+'mi'+'\t'+'nomi'+'\n')

if args.mode == 'params':
	sys.stdout.write('model'+'\t'+'divtime'+'\t'+'m12'+'\t'+'m21'+'\n')

for model,params in best_models.iteritems():


	dadi_theta,core_params = params
	mscore_func_name = model + '_mscore'
	mscore_func=None
	msplot_func_name = model + '_msplot'
	msplot_func=None
	msparams_func_name = model+'_keyparams'
	msparams_func=None

	try:
		mscore_func = getattr(dd_models,mscore_func_name)
	except Exception as e:
		print "Warning. No mscore command for "+mscore_func_name
		pass

	try:
		msplot_func = getattr(dd_models,msplot_func_name)
	except Exception as e:
		print "Warning. No plot command for "+msplot_func_name
		pass

	try:
		msparams_func = getattr(dd_models,msparams_func_name)
	except Exception as e:
		print "Warning. No msparams command for "+msparams_func_name
		pass


	Nref = float(dadi_theta)/(Leff*mu*4)


	if args.mode == 'ms' and (mscore_func!=None):

		core = mscore_func([float(f) for f in core_params.split(',')])



		sub_dict = {'theta':4*Nref*mu*L,
					'recomb':4*Nref*rho*L,
					'L':L,
					'P':P,
					'core':core}


		ms_command = "300 1 -t %(theta)f -I 2 60 240 -r %(recomb)f %(L)i -p %(P)i %(core)s" % sub_dict

		print model+"\t"+ms_command

	if args.mode == 'plot' and (msplot_func!=None):

		sys.stdout.write(msplot_func(Nref,5,[float(f) for f in core_params.split(',')]))

	if args.mode == 'params' and (msparams_func!=None):
		# import pdb;pdb.set_trace()
		sys.stdout.write(msparams_func(Nref,5,[float(f) for f in core_params.split(',')]))


















