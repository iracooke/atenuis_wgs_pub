import sys
import os
import numpy
import dadi
import pylab
import dadi.Numerics
import dd_models
import csv

best_models={}
with open('best_models.tsv') as tsvfile:
	reader = csv.reader(tsvfile, delimiter='\t')
	best_models = { row[0]:row[6] for row in reader}




snps = "dadi.thin1k.txt"
dd = dadi.Misc.make_data_dict(snps)
pop_ids=["mi", "nomi"]
proj = [45,180]
fs = dadi.Spectrum.from_data_dict(dd, pop_ids=pop_ids, projections = proj, polarized = False)

fs.mask[1,0]=True
fs.mask[0,1]=True

pts = [120,160,200]


for model,params in best_models.iteritems():
	func = getattr(dd_models,model)

	print "Plotting residuals for " + func.func_name

	func_ex = dadi.Numerics.make_extrap_log_func(func)

	popt = [float(f) for f in params.split(',')]

	print "With optimal params " + str(popt)

	model = func_ex(popt, proj, pts)

	ll_model = dadi.Inference.ll_multinom(model, fs)
	print 'Maximum log composite likelihood: {0}'.format(ll_model)
	# The optimal value of theta given the model.
	theta = dadi.Inference.optimal_sfs_scaling(model, fs)
	print 'Optimal value of theta: {0}'.format(theta)

	fig = pylab.figure(1)
	fig.clear()

	modelp = dadi.Inference.optimally_scaled_sfs(model, fs)

	dadi.Plotting.plot_2d_comp_Poisson(modelp, fs,vmin=1,resid_range=10, pop_ids = ('MI',"North"), show=False)

	#dadi.Plotting.plot_2d_comp_multinom(model, fs, vmin=1, resid_range=10,pop_ids =('MI','North'), show=False)

	fig.savefig('residuals/'+func.func_name+'.png')                               

	