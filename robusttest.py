from qn_calc import qn_calc
from hlqest import hlqest

import numpy as np
ans_qn=[]
ans_std=[]
ans_med=[]
ans_mean=[]
ans_hl=[]
ans_mad=[]
for i in range(10000):
    rr=np.random.normal(size=100)
    ans_qn.append(qn_calc(rr))
    ans_std.append(np.std(rr))
    med=np.median(rr)
    ans_med.append(med)
    ans_mad.append(1.4826*np.median(np.abs(rr-med)))
    ans_mean.append(np.mean(rr))
    ans_hl.append(hlqest(rr))
print('Mean Qn value: ',np.mean(ans_qn))
print('Mean Std value: ',np.mean(ans_std))
print('Mean MADvalue: ',np.mean(ans_mad))
rq=np.std(ans_qn)
rs=np.std(ans_std)
rmad=np.std(ans_mad)
print('Std of Qn: ',rq)
print('Std of MAD: ',rmad)
print('Std of std: ',rs,(len(rr)*2.)**-0.5)
print('Efficiency of Qn:',rs/rq)
print('Efficiency of MAD:',rs/rmad)

print('Mean median value: ',np.mean(ans_med))
print('Mean HL value: ',np.mean(ans_hl))
print('Mean mean value: ',np.mean(ans_mean))

rq=np.std(ans_med)
rs=np.std(ans_mean)
rl=np.std(ans_hl)
print('Std of med: ',rq,(len(rr)*np.pi/2)**-0.5)
print('Std of hl: ',rl)
print('Std of mean: ',rs,(len(rr)*1.)**-0.5)
print('Efficiency of med:',rs/rq,(np.pi/2)**-0.5)
print('Efficiency of hl:',rs/rl)
