#!/usr/bin/python3
# vim: set et ts=4 sw=4 tw=100:
 
# imported file name format: cval_D1.900_T1.540002.dat
# having data in it as: x/L \ c from extrapolation \  errorbar"

from glob import glob
import numpy as np
import mpmath as mpm
from scipy.optimize import curve_fit


def DedekindEta(tau):
    q = mpm.exp(2 * np.pi * 1j * tau)
    dede = mpm.qp(q) * (q**(1.0/24.0))
    return dede
    
    
    
def Jfunc(x, param):
    res = mpm.fp.fabs(DedekindEta(x*1j) * DedekindEta((1-x)*1j)/ DedekindEta(0.5j)**2)
    return (np.log(res) * param * 0.5)

#def func(x, yintercept, slope):
    #return slope*x + yintercept


files = glob('cval*')
Temps=[]
D=[]
fdict = {}
for i in files:
    t = i.split('_')
    temp =t[2][1:]
    dval = float(t[1] [1:])
    #print(temp)
    #print(dvals)
    # T=float(temp[0:4])
    if temp not in Temps:
        Temps.append(temp)
    if dval not in D:
        D.append(dval)
    fdict[('{0:.4f}'.format(dval),temp)] = i

for temp in Temps:
        tt = temp.split('.')
        temp1 = tt[0][0:]
        temp2 = tt[1][0:]
        T = float(temp1 + "." + temp2)
        #print(T)
        resarr = []
        for dval in D:
            fdata = np.loadtxt(fdict[('{0:.4f}'.format(dval),temp)])
            val1  = fdata[:,1]
            l=len(val1)
            xarray = np.array(list( fdata[i,0] for i in range(0,l)))
            yarray = np.array(list(fdata[i,1] for i in range(0,l)))
            popt, pcov = curve_fit(Jfunc, xarray, yarray)
            perr = np.sqrt(np.diag(pcov)) # this is error estimate array for fitting params
            #print(r,popt[0], perr[0])

            resarr.append( [r,popt, perr] )

        res = np.reshape(resarr, (-1,3))
#popt[0] = yintercept     
#popt[1] = slope
        np.savetxt("fit_c_D{0:.3f}_T{1:.3f}.dat".format(dval,T), resarr,fmt=' % 15.4E % 15.4E % 15.4E')


        #for i in range(0,l-1):
            #xdata.append(1/fdata[i,0])
            #ydata.append(fdata[i,4]])



 
