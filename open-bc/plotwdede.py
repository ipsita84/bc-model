# vim: set et ts=4 sw=4 tw=80:
"""
File format
MI2_L32_B0.82101

Requires mpmath module to compute the DedekindEta function at realtime
instead of reading from datafile.
"""


import numpy as np
import mpmath as mpm
#import glob

from matplotlib import use
use('agg')

import matplotlib.pyplot as plt
# plt.rc('text', usetex=True)

from glob import glob
#from matplotlib import pyplot as plt
#import numpy as np

def DedekindEta(tau):
    q = mpm.exp(2 * np.pi * 1j * tau)
    dede = mpm.qp(q) * (q**(1.0/24.0))
    return dede

def JFunc(x, param):
    res = mpm.fp.fabs(DedekindEta(x*1j) * DedekindEta((1-x)*1j)
                      / DedekindEta(0.5j)**2)
    return (np.log(res) * param * 0.5)

def main():
<<<<<<< HEAD
    D = 1.965
=======
    D = 1.9
>>>>>>> cef3a894f59b9f9660f9741f902ada03c3ea1483
    files = glob('MI*')
    sizes = []
    betas = []
    fdict = {}
    for i in files:
        t = i.split('_')
        s = int(t[1][1:])
        b = t[2][1:]
        if s not in sizes:
            sizes.append(s)
        if b not in betas:
            betas.append(b)
        fdict[(s,b)] = i


    for b in betas:
        plt.clf()
<<<<<<< HEAD
        plt.title('D= %f, T = %f' % (D, 1./float(b)))
=======
        plt.title('D = %f, T = %f' % (D, 1./float(b)))
>>>>>>> cef3a894f59b9f9660f9741f902ada03c3ea1483
        plt.xlabel('x/L')
        plt.ylabel('I_2(x/L) - I_2(1/2)')

#       DedekindEta data
        xarray = np.array([(i*1./100) for i in range(1,100)])
        ypts   = np.array([JFunc(x, param=0.7) for x in xarray])
        plt.plot(xarray, ypts, linestyle='solid',
          linewidth=2, label="CFT curve")
       
        for s in sizes:
            d = np.loadtxt(fdict[(s,b)])
            xarray = np.array([i*1./s for i in range(s+1)])
            #plt.errorbar(xarray, d[:,0], yerr=d[:,1],linestyle='none',
                           #marker='o', markersize=4, label='L = {:n}'.format(s))
            plt.errorbar(xarray, d[:,0], yerr=d[:,1],linestyle='dotted', 
                          marker='o',markersize=4, label='L = {:n}'.format(s))



        plt.xlim([0,0.5])
        plt.ylim([-0.5,0])
        plt.legend(loc='lower right', handlelength=3, numpoints=1)
        plt.show()
        plt.savefig('test1_%f.pdf' % (1./float(b)) )

if __name__ == "__main__":
    main()