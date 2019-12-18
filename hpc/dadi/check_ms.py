import dadi
import pylab
import os, sys
import dd_models
import csv

def get_func(func_name):
	try:
		func = getattr(dd_models,func_name)
	except Exception as e:
		print "Unable to find function named " + func_name
		return None

	return func



def test_mscore_for_model(model_name,params,fs):
	func = get_func(model_name)
	mscore_func_name = model_name+'_mscore'

	mscore_func = get_func(mscore_func_name)

	if( (mscore_func==None) or (func == None)):
		return

	print "Testing "+func.func_name+" and "+mscore_func.func_name

#	import pdb;pdb.set_trace()
	popt = [float(f) for f in params.split(',')]

	func_ex = dadi.Numerics.make_extrap_log_func(func)
	model = func_ex(popt,proj,pts)


	ll_model = dadi.Inference.ll_multinom(model, fs)
	print('Maximum log composite likelihood: {0}'.format(ll_model))
	# The optimal value of theta given the model.
	theta = dadi.Inference.optimal_sfs_scaling(model, fs)
	print('Optimal value of theta: {0}'.format(theta))

	# Let's generate some data using ms, if you have it installed.
	mscore = mscore_func(popt)
	# I find that it's most efficient to simulate with theta=1, average over many
	# iterations, and then scale up.
	mscommand = dadi.Misc.ms_command(1., proj, mscore, int(1e5))
	# If you have ms installed, uncomment these lines to see the results.

	return_code = os.system('{0} > test.msout'.format(mscommand))
	# We check the return code, so the script doesn't crash if you don't have ms
	# installed
	msdata = dadi.Spectrum.from_ms_file('test.msout')

	fig = pylab.figure(1)
	fig.clear()


	modelp = dadi.Inference.optimally_scaled_sfs(model, fs)

	dadi.Plotting.plot_2d_comp_Poisson(modelp, theta*msdata.fold(),vmin=1,
		pop_ids = pop_ids, show=False)

	fig.savefig('mscore_residuals/'+func.func_name+'.png') 



snps = "dadi.thin1k.txt"
pop_ids=["mi", "nomi"]
proj = [45,180]
pts = [120,160,200]
dd = dadi.Misc.make_data_dict(snps)
fs = dadi.Spectrum.from_data_dict(dd, pop_ids=pop_ids, projections = proj, polarized = False)
fs.mask[0,1] = True
fs.mask[1,0] = True
	

best_models={}
with open('best_models.tsv') as tsvfile:
	reader = csv.reader(tsvfile, delimiter='\t')
	best_models = { row[0]:row[6] for row in reader}

if len(sys.argv) > 1:
	model=sys.argv[1]
	params = best_models[model]
	test_mscore_for_model(model,params,fs)
else:
	for model,params in best_models.iteritems():
		print "------------------------"
		test_mscore_for_model(model,params,fs)
		print "----------- Done ----------"

