"""
File format
MI2_L32_B0.82101
"""


import numpy as np
#import glob

from matplotlib import use
use('agg')

import matplotlib.pyplot as plt
plt.rc('text', usetex=True)

from glob import glob
#from matplotlib import pyplot as plt
#import numpy as np

def main():
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
        plt.title('$T = %f$' % (1./float(b)))
        plt.xlabel('$x/L$')
        plt.ylabel('$I_2(x/L) - I_2(1/2)$')
        #fun = lambda x: 0.7*np.log(np.sin(x*np.pi))
        fun = lambda x: 0.7*np.log(np.sin(x*np.pi))
        for s in sizes:
            d = np.loadtxt(fdict[(s,b)])
            xarray = np.array([i*1./s for i in range(s+1)])
            plt.errorbar(xarray, d[:,0], d[:,1], label='$L = {:n}$'.format(s))

#       DedekindEta data
        ddkdata = np.loadtxt("dedekindeta_07.dat");
        plt.plot(ddkdata[:,0], ddkdata[:,1], 'r--',label="Exact curve")


        plt.xlim([0,0.5])
        plt.ylim([-0.5,0])
        plt.legend(loc='upper left', handlelength=3, numpoints=1)
        plt.show()
        plt.savefig('test_%f.pdf' % (1./float(b)) )

if __name__ == "__main__":
    main()
