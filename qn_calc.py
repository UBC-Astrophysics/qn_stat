import numpy as np
import subprocess

def qn_calc_old(x):
    np.save("dummy",x)
    output=subprocess.check_output(["qn_stat-master/qn_calc","dummy.npy"])
    return float(output)

import ctypes
_sum = ctypes.CDLL('qn_stat-master/qn_stat.so')
_sum.qn_calc.argtypes = (ctypes.POINTER(ctypes.c_double), ctypes.c_int)
_sum.qn_calc.restype = ctypes.c_double

def qn_calc(x):
    global _sum
    lx=len(x)
    array_type=ctypes.c_double * lx
    result=_sum.qn_calc(array_type(*x),ctypes.c_int(lx))
    return result
