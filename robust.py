import ctypes
_sum = ctypes.CDLL('robust_c.so')
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
