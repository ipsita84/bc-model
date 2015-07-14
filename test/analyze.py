"""
New analysis scripy for parallel tempering in
the BLume-Capel model

file format:
L008/R004/ratio_0.65681445000000005
"""

import numpy as np
from glob import iglob
from itertools import izip

class NoDataError(Exception):
    pass

def getAvg(filename):
    fin = open(filename)
    x = 0
    x2 = 0
    c = 0
    for i in fin:
        t = float(i)
        x += t
        x2 += t*t
        c += 1
    fin.close()
    if c==0:
        raise NoDataError("File was empty")
    x /= c
    x2 /= c
    return x, ((x2-x**2)/c)**0.5

def main():
    file_dict = {}
    #sizes = [8, 16, 32]
    sizes = [8, 16]
    betas = []
    for fnames in iglob('L*/R*/ratio*'):
        t = fnames.split('/')[0]
        size = int(t.split('L')[1])
        t = fnames.split('/')[1]
        ratio = int(t.split('R')[1])
        t = fnames.split('/')[2]
        beta = t.split('_')[1][:7]
        file_dict[(size,beta,ratio)] = fnames
        if beta not in betas:
            betas.append(beta)
    for S in sizes:
        for b in betas:
            S2  = [0 for i in range(S+1)]
            S2e = [0 for i in range(S+1)]
            noData = False
            for R in range(S):
                try:
                    x, dx = getAvg(file_dict[(S,b,R)])
                except NoDataError:
                    noData = True
                    break
                S2[R+1]  = S2[R] - np.log(x)
                S2e[R+1] = (S2e[R]**2 + (dx/x)**2)**0.5
            if not noData:
                fout = open('S2_L%d_B%s' % (S,b),'w')
                for x,dx in izip(S2, S2e):
                    fout.write('%1.8f %1.8f\n' % (x,dx))
                fout.close()
                MI  = [0 for i in S2]
                MIe = [0 for i in S2]
                MIMax = 0
                MIeMax = 0
                for R in range(S-1):
                    f = R+1
                    r = S-f
                    MI[f]  = S2[f] + S2[r] - S2[S]
                    v1 = S2e[S]**2 - S2e[f]**2 + S2e[r]**2
                    v2 = S2e[S]**2 - S2e[r]**2 + S2e[f]**2
                    MIe[f] = v1**0.5 if v1 < v2 else v2**0.5
                    if MI[f] > MIMax:
                        MIMax = MI[f]
                    if MIe[f] > MIeMax:
                        MIeMax = MIe[f]
                # Rectify the error, shift the middle point to zero
                for R in range(S+1):
                    MI[R] -= MIMax
                    MIe[R] = (MIeMax**2 - MIe[R]**2)**0.5
                fout = open('MI2_L%d_B%s' % (S,b),'w')
                for x,dx in izip(MI, MIe):
                    fout.write('%1.8f %1.8f\n' % (x,dx))
                fout.close()

if __name__ == "__main__":
    main()
