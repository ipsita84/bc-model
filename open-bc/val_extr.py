#!/usr/bin/python3
# vim: set et ts=4 sw=4 tw=100:
 
# imported file name format: extra_x0.357140_D1.900000_T1.560014.dat
# having data in it as: L \ x \ x/L \ 1/L \ I2shifted \ errorbar")

from glob import glob
import numpy as np
from scipy.optimize import curve_fit


def func(x, yintercept, slope):
    return slope*x + yintercept


D=1.9
files = glob('extra*')
Temps=[]
ratios=[]
fdict = {}
for i in files:
    t = i.split('_')
    temp =t[3][1:]
    #print(temp)
    ratio = float(t[1] [1:])
#print(ratio)
# T=float(temp[0:4])
    if temp not in Temps:
        Temps.append(temp)
    if ratio not in ratios:
        ratios.append(ratio)
    fdict[(ratio,temp)] = i

for temp in Temps:
        resarr = []
        for r in ratios: 
            fdata = np.loadtxt(fdict[(ratio,temp)])
            val1  = fdata[:,2]
            l=len(val1)
            xarray = np.array(list(1/fdata[i,0] for i in range(0,l)))
            yarray = np.array(list(fdata[i,4] for i in range(0,l)))
            popt, pcov = curve_fit(func, xarray, yarray)
            perr = np.sqrt(np.diag(pcov)) # this is error estimate array for fitting params
            print(r,popt[0], perr[0])

            resarr.append( [r,popt[0], perr[0]] )

        res = np.reshape(resarr, (-1,3))
#popt[0] = yintercept     
#popt[1] = slope
        np.savetxt("cval_D%f_T%s"%( D,temp), resarr,fmt=' % 15.4E % 15.4E % 15.4E')


        #for i in range(0,l-1):
            #xdata.append(1/fdata[i,0])
            #ydata.append(fdata[i,4]])



 
