"""
Analyze script for the Blume Capel model
"""

import numpy as np
import glob

def getmean(filename):
    fin = open(filename)
    x = 0
    x2 = 0
    c = 0
    for i in fin:
        t = float(i)
        x += t
        x2 += t*t
        c += 1
    try:
        x /= c
        x2 /= c
    except:
        print filename
        raise
    return x, ((x2-x**2)/c)**0.5

def main():
    ftemplate = 'L%03d/ratio_%d'
    sizes = [8,16,32,64]
    #sizes = [8,16,32]
    for S in sizes:
        s2 = [0 for i in range(S+1)]
        s2e = [0 for i in range(S+1)]
        mi2 = [0 for i in range(S+1)]
        mi2e = [0 for i in range(S+1)]
        for R in range(S):
            x, dx = getmean(ftemplate % (S,R))
            s2[R+1] = s2[R] - np.log(x)
            s2e[R+1] = (s2e[R]**2 + (dx/x)**2)**0.5
        for R in range(S-1):
            mi2[R+1] = s2[R+1] + s2[-1-R-1] - s2[-1]
            e1 = (s2e[R+1]**2 + (s2e[-1]**2 - s2e[-1-R-1]**2)**2)**0.5
            e2 = (s2e[-1-R-1]**2 + (s2e[-1]**2 - s2e[R+1]**2)**2)**0.5
            mi2e[R+1] = e1 if e1 < e2 else e2
        fout = open('S2_%03d.dat' % S,'w')
        for i in range(S+1):
            fout.write('%1.8f %1.8f\n' % (mi2[i], mi2e[i]))
        fout.close()
        fout = open('ent2_%03d.dat' % S,'w')
        for i in range(S+1):
            fout.write('%1.8f %1.8f\n' % (s2[i], s2e[i]))
        fout.close()

def findtemps(foldername):
    files = glob.glob(foldername + 'ratio_*')
    temps = []
    for i in files:
        temps.append(i.split('ratio_')[1])
    return temps

def main():
    ftemplate = 'L%03d/R%03d/ratio_%s'
    folder = 'L%03d/R%03d/'
    sizes = [8,16,32,64]
    betas = findtemps(folder % (8,0))
    data = {}
    for S in sizes:
        for B in betas:
            s2 = [0 for i in range(S+1)]
            s2e = [0 for i in range(S+1)]
            mi2 = [0 for i in range(S+1)]
            mi2e = [0 for i in range(S+1)]
            for R in range(S):
                x, dx = getmean(ftemplate % (S,R,B))
                s2[R+1] = s2[R] - np.log(x)
                s2e[R+1] = (s2e[R]**2 + (dx/x)**2)**0.5
            for R in range(S-1):
                mi2[R+1] = s2[R+1] + s2[-1-R-1] - s2[-1]
                e1 = (s2e[R+1]**2 + (s2e[-1]**2 - s2e[-1-R-1]**2)**2)**0.5
                e2 = (s2e[-1-R-1]**2 + (s2e[-1]**2 - s2e[R+1]**2)**2)**0.5
                mi2e[R+1] = e1 if e1 < e2 else e2
            fout = open('S2_%03d_%s.dat' % (S,B),'w')
            for i in range(S+1):
                fout.write('%1.8f %1.8f\n' % (mi2[i], mi2e[i]))
            fout.close()

if __name__ == "__main__":
    main()
