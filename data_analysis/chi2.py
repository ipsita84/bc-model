# imported file name format: cval_D1.900_T1.540002.dat
# having data in it as: x/L \ c from extrapolation \  errorbar"

# output format: D_c T_c chi_sq-value

import numpy as np
import glob
from scipy.optimize import leastsq

class MyException(Exception):
    pass

class linFit():
    def __init__(self):
        self.data = np.loadtxt("dedekindeta_07.dat")
        self.minx = self.data[0,0]
        self.maxx = self.data[-1,0]
    
    def val(self,x):
        if x<self.minx:
            raise MyException("Out of domain, %f < %f" % (x, self.minx))
        if x>self.maxx:
            raise MyException("Out of domain, %f > %f" % (x, self.maxx))
        for i in range(len(self.data)):
            if x < self.data[i,0]:
                x1 = self.data[i-1,0]
                x2 = self.data[i,0]
                y1 = self.data[i-1,1]
                y2 = self.data[i,1]
                return y1 + (y2-y1)*(x-x1)/(x2-x1)

def err(X,Y,dY):
    return (((X-Y)/dY)**2).sum()
                
def main():
    fout = open('chi2.txt','w')
    dede = linFit() 
    names = glob.glob("cval_*")
    res = []
    Temps=[]
    D=[]
    for fname in names: 
        t = fname.split('_')
        temp =t[2][1:]
        dval = float(t[1] [1:])
        tt = temp.split('.')
        temp1 = tt[0][0:]
        temp2 = tt[1][0:]
        T = float(temp1 + "." + temp2)
        d = np.loadtxt(fname)
        fitvals = [dede.val(i) for i in d[:,0]]
        fitvals = np.array(fitvals)
        chi2 = err(fitvals, d[:,1], d[:,2])
        #fout.write("%f %f %f\n" % (D, T ,chi2))
        fout.write( "{0:.3f} {1:.3f} {2:.3f}\n".format(dval,T/2, chi2))
        #T = float(fname.split('T')[1].split('.dat')[0])
        res.append((dval,T,chi2))
    fout.close()
    res = np.array(res)
    fitf = lambda C,x,y: C[0]*(C[1] - x)**2 + C[2] - y
    print (res[2:])
    x = leastsq(fitf,[1.,1.58,10.],(res[1:,0],res[1:,1]))
    print (x)
        
if __name__ == "__main__":
    main()
