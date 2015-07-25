#!/usr/bin/python3
# vim: set et ts=4 sw=4 tw=100:
# header="L \ x \ x/L \ 1/L \ I2shifted \ errorbar")

import numpy as np

def main():
    fdata = np.loadtxt("tot.dat")
    val1  = fdata[:,2]
    l=len(val1)
    print(l)
    thres=1E-10

    for i in range(0,l-1):
        resarr = []
        len_resarr = 1
        comval=val1[i] 
        resarr.append([fdata[i,0], fdata[i,1],comval,
         1/fdata[i,0], fdata[i,3], fdata[i,4]])
         
        for j in range(i+1,l-1):
            if (np.abs(val1[j]-comval) < thres):
                len_resarr = len_resarr + 1 
                resarr.append([fdata[j,0], fdata[j,1], comval, 
                1/fdata[j,0], fdata[j,3], fdata[j,4]])
                
        #len_resarr = len(resarr[:])
        #print(len_resarr)
        if (len_resarr > 1):
             res = np.reshape(resarr,(-1, 6))
             np.savetxt("extra_{:0.3f}.dat".format(comval), res,
                        fmt='%3i %3i % 15.4E % 15.4E % 15.4E % 15.4E')

if __name__ == "__main__":
    main()
