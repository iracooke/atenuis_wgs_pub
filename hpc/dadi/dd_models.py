import dadi
import math
from dadi import Numerics, PhiManip, Integration
from dadi.Spectrum_mod import Spectrum
import Optimize_Functions

def keyparams_string(model,Nref,gen_time,T,m12,m21):
    keyparams = "%(model)s\t%(dt)f\t%(m12)f\t%(m21)f\n"

    sub_dict = {
    'model':model,
    'dt':T*Nref*gen_time,
    'm12':m12/(2*Nref),
    'm21':m21/(2*Nref)
    }

    return keyparams % sub_dict


###############
#
# no_mig
#
###############


def no_mig(params, ns, pts):
    """
    Split into two populations, no migration.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T: Time in the past of split (in units of 2*Na generations) 
    """
    nu1, nu2, T = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T, nu1, nu2, m12=0, m21=0)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def no_mig_mscore(params):
    """
    ms core command corresponding to asym_mig
    """
    nu1, nu2, T = params

    command = "-n 1 %(nu1)f -n 2 %(nu2)f "\
            "-ej %(t1)f 2 1 "\
            "-en %(t1)f 1 1"

    # There are several factors of 2 necessary to convert units between dadi
    # and ms.
    sub_dict = {
    'nu1':nu1, 'nu2':nu2, 
    't1':T/2}

    return command % sub_dict

def no_mig_msplot(Nref,gen_time,params):
    nu1, nu2, T = params

    plot_coords = "%(model)s\t0\t%(nu1)f\t%(nu2)f\n"\
                "%(model)s\t%(t1)f\t%(nu0)f\t%(nu0)f\n" \
                "%(model)s\t%(t0)f\t%(nu0)f\t%(nu0)f\n"

    sub_dict = {
    'nu1':nu1*Nref, 'nu2':nu2*Nref, 
    't1':T*Nref*gen_time, 't0':(T+0.5)*Nref*gen_time,
    'nu0':Nref, 'model':'no_mig'}

    return plot_coords % sub_dict  

def no_mig_keyparams(Nref,gen_time,params):
    nu1, nu2, T = params
    return keyparams_string('no_mig',Nref,gen_time,T,0,0)

###############
#
# sym_mig
#
###############


def sym_mig(params, ns, pts):
    """
    Split into two populations, with symmetric migration.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T: Time in the past of split (in units of 2*Na generations) 
    m: Migration rate between populations (2*Na*m)
    """
    nu1, nu2, T, m = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T, nu1, nu2, m12=m, m21=m)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def sym_mig_mscore(params):
    """
    ms core command corresponding to asym_mig
    """
    nu1, nu2, T, m = params

    command = "-n 1 %(nu1)f -n 2 %(nu2)f "\
            "-ma x %(m12)f %(m21)f x "\
            "-ej %(t1)f 2 1 "\
            "-en %(t1)f 1 1"

    # There are several factors of 2 necessary to convert units between dadi
    # and ms.
    sub_dict = {
    'nu1':nu1, 'nu2':nu2, 
    'm12':2*m,'m21':2*m,
    't1':T/2}

    return command % sub_dict

def sym_mig_msplot(Nref,gen_time,params):
    nu1, nu2, T, m = params

    plot_coords = "%(model)s\t0\t%(nu1)f\t%(nu2)f\n"\
                "%(model)s\t%(t1)f\t%(nu0)f\t%(nu0)f\n" \
                "%(model)s\t%(t0)f\t%(nu0)f\t%(nu0)f\n"

    sub_dict = {
    'nu1':nu1*Nref, 'nu2':nu2*Nref, 
    't1':T*Nref*gen_time, 't0':(T+0.5)*Nref*gen_time,
    'nu0':Nref, 'model':'sym_mig'}

    return plot_coords % sub_dict    


def sym_mig_keyparams(Nref,gen_time,params):
    nu1, nu2, T, m = params
    return keyparams_string('sym_mig',Nref,gen_time,T,m,m)

###############
#
# asym_mig
#
###############

def asym_mig(params, ns, pts):
    """
    Split into two populations, with different migration rates.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T: Time in the past of split (in units of 2*Na generations) 
    m12: Migration from pop 2 to pop 1 (2*Na*m12)
    m21: Migration from pop 1 to pop 2
    """
    nu1, nu2, T, m12, m21 = params
    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    phi = Integration.two_pops(phi, xx, T, nu1, nu2, m12=m12, m21=m21)
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    
    return fs    

def asym_mig_mscore(params):
    """
    ms core command corresponding to asym_mig
    """
    nu1, nu2, T, m12, m21 = params

    command = "-n 1 %(nu1)f -n 2 %(nu2)f "\
            "-ma x %(m12)f %(m21)f x "\
            "-ej %(t1)f 2 1 "\
            "-en %(t1)f 1 1"

    # There are several factors of 2 necessary to convert units between dadi
    # and ms.
    sub_dict = {
    'nu1':nu1, 'nu2':nu2, 
    'm12':2*m12,'m21':2*m21,
    't1':T/2}

    return command % sub_dict

def asym_mig_msplot(Nref,gen_time,params):
    nu1, nu2, T, m12, m21 = params

    plot_coords = "%(model)s\t0\t%(nu1)f\t%(nu2)f\n"\
                "%(model)s\t%(t1)f\t%(nu0)f\t%(nu0)f\n" \
                "%(model)s\t%(t0)f\t%(nu0)f\t%(nu0)f\n"

    sub_dict = {
    'nu1':nu1*Nref, 'nu2':nu2*Nref, 
    't1':T*Nref*gen_time, 't0':(T+0.5)*Nref*gen_time,
    'nu0':Nref, 'model':'asym_mig'}

    return plot_coords % sub_dict    


def asym_mig_keyparams(Nref,gen_time,params):
    nu1, nu2, T, m12, m21 = params
    return keyparams_string('asym_mig',Nref,gen_time,T,m12,m21)

###############
#
# priorsize_asym_mig
#
###############

def priorsize_asym_mig(params, ns, pts):
    """
    Size change followed by split with asymmetric migration

    nua: First Size of population before split.
    T1: Duration of first time before split
    nu1b: Size of population 1 after split.
    nu2b: Size of population 2 after split.
    T2: Time in the past of split (units of 2*Na generations)
    m12: Migration from pop 2 to pop 1 (2*Na*m12)
    m21: Migration from pop 1 to pop 2
    """
    nua, T1,nu1b, nu2b, T2, m12, m21 = params
    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx)

    phi = Integration.one_pop(phi, xx, T1, nua)

    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    phi = Integration.two_pops(phi, xx, T2, nu1b, nu2b, m12=m12, m21=m21)
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    
    return fs

def priorsize_asym_mig_mscore(params):
    """
    ms core command corresponding to priorsize_asym_mig
    """
    nua, T1,nu1b, nu2b, T2, m12, m21 = params

    command = "-n 1 %(nu1b)f -n 2 %(nu2b)f "\
            "-ma x %(m12)f %(m21)f x "\
            "-ej %(t2)f 2 1 "\
            "-eN %(t2)f %(nua)f "\
            "-eN %(t1)f 1"

    # There are several factors of 2 necessary to convert units between dadi
    # and ms.
    sub_dict = {
    'nu1b':nu1b, 'nu2b':nu2b, 
    'nua':nua,
    'm12':2*m12,'m21':2*m21,
    't2':(T2)/2,
    't1':(T2+T1)/2}

    return command % sub_dict


def priorsize_asym_mig_msplot(Nref,gen_time,params):
    nua, T1,nu1b, nu2b, T2, m12, m21 = params

    plot_coords = "%(model)s\t0\t%(nu1b)f\t%(nu2b)f\n"\
                "%(model)s\t%(t2)f\t%(nua)f\t%(nua)f\n"\
                "%(model)s\t%(t21)f\t%(nu0)f\t%(nu0)f\n"\
                "%(model)s\t%(t0)f\t%(nu0)f\t%(nu0)f\n"

    sub_dict = {
    'nu1b':nu1b*Nref, 'nu2b':nu2b*Nref, 
    'nua':nua*Nref, 
    'nu0':Nref,       
    't2':T2*Nref*gen_time,
    't21':(T2+T1)*Nref*gen_time,
    't0':(T2+T1+0.5)*Nref*gen_time,
    'model':'priorsize_asym_mig'}

    return plot_coords % sub_dict    

def priorsize_asym_mig_keyparams(Nref,gen_time,params):
    nua, T1,nu1b, nu2b, T2, m12, m21 = params
    return keyparams_string('priorsize_asym_mig',Nref,gen_time,(T2),m12,m21)


###############
#
# asym_mig_size
#
###############

def asym_mig_size(params, ns, pts):
    """
    Split followed by two size changes with continuous asymmetric migration

    nu1a: Size of population 1 after split in first time interval.
    nu2a: Size of population 2 after split.
    T1:  First time interval after split (units of 2*Na generations)
    nu1b: Size of population 1 after split in second time interval.
    nu2b: Size of population 2 after split.
    T2:  Second time interval after split (units of 2*Na generations)
    m12: Migration from pop 2 to pop 1 (2*Na*m12)
    m21: Migration from pop 1 to pop 2
    """
    nu1a, nu2a, T1,nu1b, nu2b, T2, m12, m21 = params
    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx)

    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    phi = Integration.two_pops(phi, xx, T1, nu1a, nu2a, m12=m12, m21=m21)

    phi = Integration.two_pops(phi, xx, T2, nu1b, nu2b, m12=m12, m21=m21)
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    
    return fs


def asym_mig_size_mscore(params):
    """
    ms core command corresponding to asym_mig_size
    """
    nu1a, nu2a, T1,nu1b, nu2b, T2, m12, m21 = params

    command = "-n 1 %(nu1b)f -n 2 %(nu2b)f "\
            "-ma x %(m12)f %(m21)f x "\
            "-en %(t2)f 1 %(nu1a)f -en %(t2)f 2 %(nu2a)f "\
            "-ej %(t1)f 2 1 "\
            "-eN %(t1)f 1"

    # There are several factors of 2 necessary to convert units between dadi
    # and ms.
    sub_dict = {
    'nu1b':nu1b, 'nu2b':nu2b, 
    'nu1a':nu1a, 'nu2a':nu2a,
    'm12':2*m12,'m21':2*m21,
    't2':T2/2,
    't1':(T1+T2)/2
    }

    return command % sub_dict


def asym_mig_size_msplot(Nref,gen_time,params):
    nu1a, nu2a, T1,nu1b, nu2b, T2, m12, m21 = params

    plot_coords = "%(model)s\t0\t%(nu1b)f\t%(nu2b)f\n"\
                "%(model)s\t%(t2)f\t%(nu1a)f\t%(nu2a)f\n"\
                "%(model)s\t%(t21)f\t%(nu0)f\t%(nu0)f\n" \
                "%(model)s\t%(t0)f\t%(nu0)f\t%(nu0)f\n"

    sub_dict = {
    'nu1b':nu1b*Nref, 'nu2b':nu2b*Nref, 
    'nu1a':nu1a*Nref, 'nu2a':nu2a*Nref, 
    'nu0':Nref,
    't2':T2*Nref*gen_time,
    't21':(T2+T1)*Nref*gen_time,
    't0':(T2+T1+0.5)*Nref*gen_time,
    'model':'asym_mig_size'}

    return plot_coords % sub_dict    

def asym_mig_size_keyparams(Nref,gen_time,params):
    nu1a, nu2a, T1,nu1b, nu2b, T2, m12, m21 = params
    return keyparams_string('asym_mig_size',Nref,gen_time,(T2+T1),m12,m21)

###############
#
# isolation_asym_mig
#
###############

def isolation_asym_mig_base(nu1a, nu2a, T1,nu1b, nu2b, T2, m12, m21, ns, pts):
    """
    Split followed by two size changes with continuous asymmetric migration

    nu1a: Size of population 1 after split in first time interval.
    nu2a: Size of population 2 after split.
    T1:  First time interval after split (units of 2*Na generations)
    nu1b: Size of population 1 after split in second time interval.
    nu2b: Size of population 2 after split.
    T2:  Second time interval after split (units of 2*Na generations)
    m12: Migration from pop 2 to pop 1 (2*Na*m12)
    m21: Migration from pop 1 to pop 2
    """
    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx)

    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    phi = Integration.two_pops(phi, xx, T1, nu1a, nu2a, m12=0, m21=0)

    phi = Integration.two_pops(phi, xx, T2, nu1b, nu2b, m12=m12, m21=m21)
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    
    return fs    


def isolation_asym_mig(params, ns, pts):
    """
    Split followed by two size changes with continuous asymmetric migration

    nu1a: Size of population 1 after split in first time interval.
    nu2a: Size of population 2 after split.
    T1:  First time interval after split (units of 2*Na generations)
    nu1b: Size of population 1 after split in second time interval.
    nu2b: Size of population 2 after split.
    T2:  Second time interval after split (units of 2*Na generations)
    m12: Migration from pop 2 to pop 1 (2*Na*m12)
    m21: Migration from pop 1 to pop 2
    """
    nu1a, nu2a, T1,nu1b, nu2b, T2, m12, m21 = params
    return isolation_asym_mig_base(nu1a, nu2a, T1,nu1b, nu2b, T2, m12, m21, ns, pts)


def isolation_asym_mig_mscore(params):
    """
    ms core command corresponding to asym_mig_size
    """
    nu1a, nu2a, T1,nu1b, nu2b, T2, m12, m21 = params

    command = "-n 1 %(nu1b)f -n 2 %(nu2b)f "\
            "-ma x %(m12)f %(m21)f x "\
            "-en %(t2)f 1 %(nu1a)f -en %(t2)f 2 %(nu2a)f "\
            "-ema %(t2)f 2 x 0 0 x "\
            "-ej %(t1)f 2 1 "\
            "-eN %(t1)f 1"

    # There are several factors of 2 necessary to convert units between dadi
    # and ms.
    sub_dict = {
    'nu1b':nu1b, 'nu2b':nu2b, 
    'nu1a':nu1a, 'nu2a':nu2a,
    'm12':2*m12,'m21':2*m21,
    't2':T2/2,
    't1':(T1+T2)/2
    }

    return command % sub_dict


def isolation_asym_mig_msplot(Nref,gen_time,params):
    nu1a, nu2a, T1,nu1b, nu2b, T2, m12, m21 = params

    plot_coords = "%(model)s\t0\t%(nu1b)f\t%(nu2b)f\n"\
                "%(model)s\t%(t2)f\t%(nu1a)f\t%(nu2a)f\n"\
                "%(model)s\t%(t21)f\t%(nu0)f\t%(nu0)f\n" \
                "%(model)s\t%(t0)f\t%(nu0)f\t%(nu0)f\n"

    sub_dict = {
    'nu1b':nu1b*Nref, 'nu2b':nu2b*Nref, 
    'nu1a':nu1a*Nref, 'nu2a':nu2a*Nref, 
    'nu0':Nref,
    't2':T2*Nref*gen_time,
    't21':(T2+T1)*Nref*gen_time,
    't0':(T2+T1+0.5)*Nref*gen_time,
    'model':'isolation_asym_mig'}

    return plot_coords % sub_dict    

def isolation_asym_mig_keyparams(Nref,gen_time,params):
    nu1a, nu2a, T1,nu1b, nu2b, T2, m12, m21 = params
    return keyparams_string('isolation_asym_mig',Nref,gen_time,(T2+T1),m12,m21)






