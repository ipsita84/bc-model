#!/usr/bin/python3
# vim: set ts=4 sw=4 tw=80:

from glob import glob
import numpy as np
import sys

def main():
    files = glob('MI*')
    sizes = []
    betas = []
    fdict = {}
    resarr = [];
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
        T=1./float(b)
        for s in sizes:
            d = np.loadtxt(fdict[(s,b)])
            r=int(s/2)
            for i in range(1,r):
                resarr.append(s, i*1./s, d[:,0], d[:,1])

    res = np.reshape(resarr,(-1, 4))
# Write data to file using numpy.savetxt
    np.savetxt("tot.dat", res,
           fmt='%3i % 15.4E % 15.4E % 15.4E')

# header="L\x/L\I2shifted\errorbar")
if __name__ == "__main__":
    main()
