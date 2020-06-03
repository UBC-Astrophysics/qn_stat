# -*- coding: utf-8 -*-

__version__ = "1.0.1.dev0"
__author__ = "Jeremy Heyl (heyl@phas.ubc.ca)"
__contributors__ = [
    # Alphabetical by first name.
    "John F. Monaghan (NC State)",
    "Florian Plaza Oñate @fplaza",
]


import os, ctypes
try:
    _sum = ctypes.CDLL(os.path.dirname(__file__)+os.sep+'robust_c.so')
except OSError as e:
    print('The qn_stat package uses a C shared library.  Please go to the qn_stat directory and execute "make".')
    raise e
else:    
    _sum.qn_calc.argtypes = (ctypes.POINTER(ctypes.c_double), ctypes.c_int)
    _sum.qn_calc.restype = ctypes.c_double

    _sum.hlqest.argtypes = (ctypes.POINTER(ctypes.c_double), ctypes.c_int)
    _sum.hlqest.restype = ctypes.c_double


def qn_calc(x):
    global _sum
    lx=len(x)
    array_type=ctypes.c_double * lx
    result=_sum.qn_calc(array_type(*x),ctypes.c_int(lx))
    return result

def hlqest(x):
    global _sum
    lx=len(x)
    array_type=ctypes.c_double * lx
    result=_sum.hlqest(array_type(*x),ctypes.c_int(lx))
    return result

def robustchi2(param,f,xdata,ydata,sdata):
    rat=(ydata-f(*((xdata,)+tuple(param))))/sdata
    return len(xdata)*(hlqest(rat)**2+qn_calc(rat)**2)

def chi2(param,f,xdata,ydata,sdata):
    rat=(ydata-f(*((xdata,)+tuple(param))))/sdata
    return sum(rat**2)

from scipy.optimize import minimize
def robustcurve_fit(f,xdata,ydata,p0=None,sigma=1,
                    absolute_sigma=False,method='BFGS',
                    return_res=False,dorobust=True):
    # the parameters work the same as for scipy.optimize.curve_fit
    m=len(xdata)
    # if you don't give starting values, try to find number of function
    # parameters
    if p0 is None:
        try:
            from inspect import signature 
            sig=signature(f)
            n=len(sig.parameters)-1
        except ImportError:
            from inspect import getargspec
            n=len(getargspec(f))-1
        p0=[1]*n
    else:
        n=len(p0)

    # use the robust fitter or just chi2
    if dorobust:
        res=minimize(robustchi2,p0,args=(f,xdata,ydata,sigma),method=method)
    else:
        res=minimize(chi2,p0,args=(f,xdata,ydata,sigma),method=method)
    if method=='BFGS':
        if absolute_sigma:
            cov=res.hess_inv*2
        else:
            cov=res.hess_inv*res.fun/(m-n)*2
    else:
        cov=0
    if return_res:
        return res.x,cov,res
    else:
        return res.x,cov
    

