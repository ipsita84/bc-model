#!/usr/bin/python3
# vim: set et ts=4 sw=4 tw=100:
# header="L \ x \ x/L \ 1/L \ I2shifted \ errorbar")

import numpy as np

# Store matched points in a global list
global matvals
matvals = []

def flcomp(val1, val2, thres=1E-10):
    if (np.abs(val1-val2) < thres):
        res = True
    else:
        res = False
    return res

def in_list(val, arr=matvals):
    if len(matvals) == 0:
        return False
    else:
        for x in matvals:
            if flcomp(val, x):
                return True

def main():
    fdata = np.loadtxt("tot.dat")
    val1  = fdata[:,2]
    l=len(val1)

    for i in range(0,l-1):
        resarr = []
        comval=val1[i]
        if in_list(comval):
            continue

        matvals.append(comval)
        resarr.append([fdata[i,0], fdata[i,1], comval,
                       1/fdata[i,0], fdata[i,3], fdata[i,4]])
         
        for j in range(i+1,l-1):
            if (flcomp(val1[j],comval)):
                if (j < len(fdata[:,0])):
                    resarr.append([fdata[j,0], fdata[j,1], comval,
                                   1/fdata[j,0], fdata[j,3], fdata[j,4]])

        res = np.reshape(resarr, (-1,6))
        len_res = len(res[:,0])
        if (len_res > 1):
                    np.savetxt("extra_{:0.3f}.dat".format(comval), resarr,
                               fmt='%3i %3i % 15.4E % 15.4E % 15.4E % 15.4E')

if __name__ == "__main__":
    main()
