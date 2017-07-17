import numpy as np
import subprocess

import ctypes
_sum = ctypes.CDLL('qn_stat.so')
_sum.qn_calc.argtypes = (ctypes.POINTER(ctypes.c_double), ctypes.c_int)
_sum.qn_calc.restype = ctypes.c_double

def qn_calc(x):
    global _sum
    lx=len(x)
    array_type=ctypes.c_double * lx
    result=_sum.qn_calc(array_type(*x),ctypes.c_int(lx))
    return result
