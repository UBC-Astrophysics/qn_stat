from robust import qn_calc, hlqest, robustcurve_fit
import numpy as np

if True:
    print('Perform Efficiency Calculations:')
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
    print('Mean Qn value:     %g' % np.mean(ans_qn))
    print('Mean Std value:    %g' % np.mean(ans_std))
    print('Mean MADvalue:     %g' % np.mean(ans_mad))
    rq=np.std(ans_qn)
    rs=np.std(ans_std)
    rmad=np.std(ans_mad)
    print('Std of Qn:         %g' % rq)
    print('Std of MAD:        %g' % rmad)
    print('Std of std:        %g (%g)' % (rs,(len(rr)*2.)**-0.5))
    print('Efficiency of Qn:  %g' % (rs/rq))
    print('Efficiency of MAD: %g' % (rs/rmad))

    print('Mean median value: %g' % np.mean(ans_med))
    print('Mean HL value:     %g' % np.mean(ans_hl))
    print('Mean mean value:   %g' % np.mean(ans_mean))

    rq=np.std(ans_med)
    rs=np.std(ans_mean)
    rl=np.std(ans_hl)
    print('Std of med:        %g (%g)' % (rq,(len(rr)*np.pi/2)**-0.5))
    print('Std of hl:         %g' % rl)
    print('Std of mean:       %g (%g)' % (rs,(len(rr)*1.)**-0.5))
    print('Efficiency of med: %g (%g)' % (rs/rq,(np.pi/2)**-0.5))
    print('Efficiency of hl:  %g' % (rs/rl))

if True:
    print('Perform fitting demonstration p=20,40,25:')

    def funk(x,p0,p1,p2):
        return np.cos(x)*p2+np.sin(x)*p1+p0

    from scipy.optimize import curve_fit
    xdata=np.linspace(0,20,201)
    ydata=funk(xdata,20.0,40.0,25.0)+np.random.normal(size=len(xdata))
    print('robust.robustcurve_fit with dorobust=False:')
    print(robustcurve_fit(funk,xdata,ydata,dorobust=False))
    print('robust.robustcurve_fit with dorobust=True:')
    print(robustcurve_fit(funk,xdata,ydata))
    print('scipy.optimize.curve_fit with starting value:')
    print(curve_fit(funk,xdata,ydata,p0=[4.0,4.0,4.0],sigma=0*xdata+1))
    print('robust.robustcurve_fit with starting value, dorobust=True, different method and return res:')
    print(robustcurve_fit(funk,xdata,ydata,p0=[4.0,4.0,4.0],sigma=1,return_res=True,method='Nelder-Mead'))
    print('Changing one of the ydata to 1e6')
    ydata[0]=1e6
    print('scipy.optimize.curve_fit:')
    print(curve_fit(funk,xdata,ydata))
    print('robust.robustcurve_fit with defaults')
    print(robustcurve_fit(funk,xdata,ydata))
    
