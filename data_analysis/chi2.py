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
    names = glob.glob("*1.900_*")
    res = []
    for fname in names:
        d = np.loadtxt(fname)
        fitvals = [dede.val(i) for i in d[:,0]]
        fitvals = np.array(fitvals)
        chi2 = err(fitvals, d[:,1], d[:,2])
        fout.write("%s Total chi^2 = %f\n" % (fname,chi2))
        T = float(fname.split('T')[1].split('.dat')[0])
        res.append((T,chi2))
    fout.close()
    res = np.array(res)
    fitf = lambda C,x,y: C[0]*(C[1] - x)**2 + C[2] - y
    print (res[1:])
    x = leastsq(fitf,[1.,1.58,10.],(res[1:,0],res[1:,1]))
    print (x)
        
if __name__ == "__main__":
    main()
