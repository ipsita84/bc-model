#!/usr/bin/python3
# vim: set ts=4 sw=4 tw=80:


import numpy as np

def main():
    val1 = []

    fdata = np.loadtxt("tot.dat")
    l=len(fdata[:,1])
    thres=1E-10


    for i in range(0,l-1):
        resarr = []
        comval=val1[i] 
        resarr.append([1/fdata[i,0], fdata[i,2], fdata[i,3]])
        for j in rangge(i+1,len-1):
            if  (math.abs(val1[j]-comval) < thres):
               resarr.append([1/fdata[j,0], fdata[j,2], fdata[j,3]])
        len_resarr = len(resarr[:])
        if (len_resarr > 1):
             res = np.reshape(resarr,(-1, 3))
             np.savetxt("extra_%f.dat"  % (comval), res,
               fmt='% 15.4E % 15.4E % 15.4E')





# header="L data errorbar")
if __name__ == "__main__":
    main()
