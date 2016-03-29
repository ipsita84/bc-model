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

import matplotlib as mpl

import matplotlib.pyplot as plt
# plt.rc('text', usetex=True)

from glob import glob

def DedekindEta(tau):
    q = mpm.exp(2 * np.pi * 1j * tau)
    dede = mpm.qp(q) * (q**(1.0/24.0))
    return dede

def JFunc(x, param):
    res = mpm.fp.fabs(DedekindEta(x*1j) * DedekindEta((1-x)*1j)
                      / DedekindEta(0.5j)**2)
    return (np.log(res) * param * 0.5)

def main():
    plt.style.use('ggplot')
    D = 1.8754
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
        sizes.sort()

    for b in betas:
        colorlist = plt.cm.gist_rainbow(np.linspace(0, 1, len(sizes)))
        colidx = 0
        plt.clf()
        plt.axes([0.1, 0.1, 0.7, 0.85])
        plt.title('$D_c$$/$$J$ = %.3f, $k$ $T_c$$/$$J$ = %.2f' 
                  % (D, 0.5/float(b)), fontsize=18)
        plt.xlabel('$\ell/L$',fontsize=18)
        plt.ylabel('$I_2$($\ell,L$) - $I_2$($L/2,L$)', fontsize=18)

#       DedekindEta data
        xarray = np.array([(i*1./100) for i in range(1,100)])
        ypts   = np.array([JFunc(x, param=0.7) for x in xarray])
        plt.plot(xarray, ypts, linestyle='solid', color='black',
                 linewidth=2, label="CFT curve")
       
        for s in sizes:
            d = np.loadtxt(fdict[(s,b)])
            xarray = np.array([i*1./s for i in range(s+1)])
            #plt.errorbar(xarray, d[:,0], yerr=d[:,1],linestyle='none',
                           #marker='o', markersize=4, label='L = {:n}'.format(s))
            plt.errorbar(xarray, d[:,0], yerr=d[:,1],linestyle='--', 
                         marker='o',markersize=(6.0 + float(s-6)/3.0),
                         color = colorlist[colidx % len(colorlist)],
                         label='L = {:n}'.format(s))
            colidx += 1

        plt.xlim([0,0.5])
        plt.xticks([0,0.1,0.2,0.3,0.4,0.5], fontsize = 11)
        plt.ylim([-0.5,0])
        plt.yticks([-0.5,-0.4,-0.3,-0.2,-0.1,0], fontsize = 11)
        plt.legend(loc='center right', bbox_to_anchor=(1.28, 0.5),
                   handlelength=4, numpoints=1, fontsize=10,
                   labelspacing=1.5)
        plt.show()
        plt.savefig('test1_%f.pdf' % (1./float(b)) )

if __name__ == "__main__":
    main()

