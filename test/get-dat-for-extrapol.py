#!/usr/bin/python3
# vim: set et ts=4 sw=4 tw=80:


import numpy as np

def main():
    fdata = np.loadtxt("tot.dat")
    val1  = fdata[:,1]
    l=len(val1)
    thres=1E-10

    for i in range(0,l-1):
        resarr = []
        comval=val1[i] 
        resarr.append([1/fdata[i,0], fdata[i,2], fdata[i,3]])
        for j in range(i+1,l-1):
            if (np.abs(val1[j]-comval) < thres):
                resarr.append([1/fdata[j,0], fdata[j,2], fdata[j,3]])
        len_resarr = len(resarr[:])
        if (len_resarr > 1):
             res = np.reshape(resarr,(-1, 3))
             np.savetxt("extra_{:0.3f}.dat".format(comval), res,
                        fmt='% 15.4E % 15.4E % 15.4E')

if __name__ == "__main__":
    main()
