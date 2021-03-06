"""
File format
MI2_L32_B0.82101
"""


from glob import glob
from matplotlib import pyplot as plt
import numpy as np

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
        plt.xlim([0,0.5])
        plt.ylim([-1,0])
        #fun = lambda x: 0.7*np.log(np.sin(x*np.pi))
        for s in sizes:
            d = np.loadtxt(fdict[(s,b)])
            xarray = np.array([i*1./s for i in range(s+1)])
            plt.errorbar(xarray, d[:,0], d[:,1])
        xarray = np.array([i*1./100. for i in range(1,101)])
        # DedekindEta data
        ddkdata = np.loadtxt("dedekindeta_07.dat");
        plt.plot(ddkdata[:,0], ddkdata[:,1], 'r--', label="Exact curve")
        plt.show()

if __name__ == "__main__":
    main()
