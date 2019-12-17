import sys
import os
import numpy
import dadi
from dadi import Numerics, PhiManip, Integration
from dadi.Spectrum_mod import Spectrum
import Optimize_Functions
import dd_models

func=dd_models.isolation_sym_mig

snps = "dadi.thin1k.txt"

dd = dadi.Misc.make_data_dict(snps)

pop_ids=["mi", "nomi"]

proj = [45,180]

fs = dadi.Spectrum.from_data_dict(dd, pop_ids=pop_ids, projections = proj, polarized = False)

fs.mask[1,0]=True
fs.mask[0,1]=True

# These are the grid point settings will use for extrapolation.
pts = [120,160,200]

prefix = "optim_d/isolation_sym_mig"

p_labels = "nu1a, nu2a, T1, nu1b, nu2b, T2, m"

upper =     [100,    100,    10,    100,    100,    10,     100]
lower =     [1e-1,      1e-1,      0,    1e-2,   1e-2,   0,      0]
p0 =        [90.8216,1.3095,0.8504,0.0109, 0.1973, 0.0084, 20]

reps = [10,10,10,10]
maxiters = [5,5,10,20]
folds = [2,1,1,1]


i=sys.argv[1]
Optimize_Functions.Optimize_Routine(fs, pts, prefix+"_{}".format(i), func.func_name, func, len(reps), len(p0), 
   fs_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, 
   reps = reps, maxiters = maxiters, folds = folds)






