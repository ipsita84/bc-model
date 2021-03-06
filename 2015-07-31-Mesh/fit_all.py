# imported file name format: cval_D1.900_T1.540002.dat
# having data in it as: x/L \ c from extrapolation \  errorbar"

# output format: D_c T_c chi_sq-value

import numpy as np
import glob
from scipy.optimize import leastsq
import operator
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl

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

# x here are the system sizes
def func(x,C0, C1): 
    return C0 + 1./x*C1

def main():
    dede = linFit()
    allfiles = glob.glob("D*/MI2*")
    #allfiles = allfiles[:10000]
    fitdata = {}
    Lset = [8,16,24,32]
    #Lset = [8,16,24]
    Dvalset = []
    betaset = []
    for f in allfiles:
        Dval = float(f.split('/')[0][1:])
        beta = float(f.split('/')[1].split('_')[2][1:])
        L = int(f.split('/')[1].split('_')[1][1:])
        if L not in Lset:
            continue
        if Dval not in Dvalset:
            Dvalset.append(Dval)
        if beta not in betaset and beta < 1.0:
            betaset.append(beta)
        d = np.loadtxt(f)
        N = len(d)/2
        d = d[N/4:N:N/4]
        if np.any(d[:,1] <= 0):
            continue
        fitdata[(Dval,beta,L)] = d[:]
    fitres = {}
    Lset = np.array(Lset)
    err = {}
    for D in Dvalset:
        for B in betaset:
            complete = True
            valid = True
            for L in Lset:
                if (D,B,L) not in fitdata:
                    complete = False
                    break
            if complete:
                allY = []
                for z in xrange(3):
                    fitY = []
                    fitdY = []
                    for L in Lset:
                        fitY.append(fitdata[(D,B,L)][z,0])
                        fitdY.append(fitdata[(D,B,L)][z,1])
                    fitY = np.array(fitY)
                    fitdY = np.array(fitdY)
                    popt, perr = curve_fit(func, Lset, fitY, [-1.,1.], fitdY, absolute_sigma=True)
                    err = (((func(Lset, *popt) - fitY)/fitdY)**2).sum()
                    if err > 100.:
                        valid = False
                        break
                    allY.append((popt[0], perr[0,0]**0.5, err))
                    #plt.clf()
                    #plt.plot(1./Lset, func(Lset, *popt))
                    #plt.plot(1./Lset, fitY)
                    #plt.title('D=%1.4f, B=%1.4f, z=%d' % (D,B,z))
                    #plt.show()
                if valid:
                    fitres[(D,B)] = np.array(allY)
    devals = np.array([dede.val(i) for i in [1./8, 2./8, 3./8]])
    chivals = []
    cvals = []
    cfit = lambda x,c: c*devals
    for i in fitres:
        chi2 = (((devals - fitres[i][:,0])/fitres[i][:,1])**2).sum()
        chivals.append((i[0],i[1],chi2))
        popt, perr = curve_fit(cfit, [], fitres[i][:,0], [1.], fitres[i][:,1], absolute_sigma=True)
        chi2 = (((cfit([],popt[0]) - fitres[i][:,0])/fitres[i][:,1])**2).sum()
        cvals.append((i[0], i[1], popt[0], perr[0,0]**0.5, chi2))
    chivals = sorted(chivals, key=operator.itemgetter(0,1))
    Dres = []
    poly2 = lambda x,a,b,c: a*x**2 + b*x + c
    for D in Dvalset:
        chimin = (-1,-1)
        for iv,v in enumerate(chivals):
            if v[0] == D:
                if chimin[0] == -1:
                    chimin = (iv,v[2])
                elif v[2] < chimin[1]:
                    chimin = (iv,v[2])
        if chimin[0] != -1:
            xvec = [chivals[i][1] for i in range(chimin[0]-1, chimin[0]+2)]
            yvec = [chivals[i][2] for i in range(chimin[0]-1, chimin[0]+2)]
            popt = np.polyfit(xvec,yvec,2)
            xmin =  -1*popt[1]/(2*popt[0])
            ymin = poly2(xmin,*popt)
            #print xvec
            #print yvec
            #print xmin
            #print ymin
            #print ''
            Dres.append((D,xmin,ymin))
    Dres = np.array(Dres)
    #plt.plot(Dres[:,0],Dres[:,1])
    #plt.xlabel('D')
    #plt.ylabel('Beta_c')
    #plt.show()
    #plt.clf()
    #plt.plot(Dres[:,0],Dres[:,2])
    #plt.xlabel('D')
    #plt.ylabel('Chi^2')
    #plt.show()
    #print cvals
    d = np.array(chivals)
    #d = np.array(cvals)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    #ax.scatter(d[:,0],d[:,1],d[:,4])
    #ax.scatter(d[:,0],d[:,1],d[:,4])
    crange = (d[:,1]-d[0,1])/(d[-1,1] - d[0,1])
    ax.scatter(d[:,0],d[:,1],d[:,2],c=crange,cmap=mpl.cm.spring)
    ax.scatter(d[:,0],d[:,1],d[:,2],c=crange,cmap=mpl.cm.spring)
    ax.set_zlim3d((0,2000))
    ax.set_xlabel('$D$')
    ax.set_ylabel('$\\beta$')
    ax.set_zlabel('Total Chi^2')
    #ax.set_ylim3d((0.35,0.42))
    plt.show()
    #for i in fitres:
    #    print fitres[i][:,2]
    #results = np.array(sorted(results, key=operator.itemgetter(0,1)))
    #np.savetxt('chi2.txt',results)

if __name__ == "__main__":
    main()
