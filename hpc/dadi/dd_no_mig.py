import sys
import os
import numpy
import dadi
from dadi import Numerics, PhiManip, Integration
from dadi.Spectrum_mod import Spectrum
import Optimize_Functions
import dd_models

func=dd_models.no_mig

snps = "dadi.thin1k.txt"

dd = dadi.Misc.make_data_dict(snps)

pop_ids=["mi", "nomi"]

proj = [45,180]

fs = dadi.Spectrum.from_data_dict(dd, pop_ids=pop_ids, projections = proj, polarized = False)

# This mask reflects our MAC>1 requirement
fs.mask[0,1] = True
fs.mask[1,0] = True

# These are the grid point settings will use for extrapolation.
pts = [120,160,200]

prefix = "optim/no_mig"

# Parameters are: (nu1, nu2, T)
p_labels = "nu_mi, nu_nomi, T_split"
upper =     [100,   100,        10]
lower =     [1e-3,  1e-3,       0]
p0 =        [1,     1,          0.1]


reps = [10,10,10,10]
maxiters = [5,5,10,20]
folds = [5,3,2,1]


i=sys.argv[1]
Optimize_Functions.Optimize_Routine(fs, pts, prefix+"_{}".format(i), func.func_name, func, len(reps), len(p0), 
   fs_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, 
   reps = reps, maxiters = maxiters, folds = folds)


