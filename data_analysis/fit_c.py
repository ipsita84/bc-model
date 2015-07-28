# imported file name format: cval_D1.900_T1.540002.dat
# having data in it as: x/L \ c from extrapolation \  errorbar"

# output format: D_c T_c chi_sq-value

import numpy as np
from glob import glob
import glob
import mpmath as mpm
from scipy.optimize import curve_fit

# imported file name format: cval_D1.900_T1.540002.dat
# having data in it as: x/L \ c from extrapolation \  errorbar"

import glob
#from glob import glob
import numpy as np
import mpmath as mpm
from scipy.optimize import curve_fit


def DedekindEta(tau):
    q = mpm.exp(2 * np.pi * 1j * tau)
    dede = mpm.qp(q) * (q**(1.0/24.0))
    return dede
    
    
    
def Jfunc(x, param):
    yarray = []
    l=len(x)
    for i in range(0,l):
        yin= mpm.fp.fabs(DedekindEta(x[i]*1j) * DedekindEta((1-x[i])*1j)/ DedekindEta(0.5j)**2)
        y= np.log(yin) * param * 0.5
        yarray.append(y)
    #res = mpm.fp.fabs(DedekindEta(x*1j) * DedekindEta((1-x)*1j)/ DedekindEta(0.5j)**2)
    return (yarray)

#def func(x, yintercept, slope):
    #return slope*x + yintercept


                
def main():
    fout = open('c.txt','w')
    names = glob.glob("cval_*")
    res = []
    Temps=[]
    D=[]
    for fname in names: 
        resarr = []
        t = fname.split('_')
        temp =t[2][1:]
        dval = float(t[1] [1:])
        tt = temp.split('.')
        temp1 = tt[0][0:]
        temp2 = tt[1][0:]
        T = float(temp1 + "." + temp2)
        fdata = np.loadtxt(fname)
        val1  = fdata[:,1]
        l=len(val1)
        xarray = np.array(list( fdata[i,0] for i in range(0,l)))
        yarray = np.array(list(fdata[i,1] for i in range(0,l)))
        popt, pcov = curve_fit(Jfunc, xarray, yarray)
        print(popt)
        perr = np.sqrt(np.diag(pcov)) # this is error estimate array for fitting params
        fout.write( "{0:.3f} {1:.3f} {2:.3f}\n".format(dval,T/2, popt[0]))
        resarr.append( [dval,T/2,popt] )
    fout.close()

        
if __name__ == "__main__":
    main()
