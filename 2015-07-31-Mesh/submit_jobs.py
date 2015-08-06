"""
Script for submitted a batch of jobs for the Blume-Capel model
Main job is run with the command

./transfer L regionSize D

L - system size
regionSize - number of strips in region A
Dval - coupling parameter D

for a complete simulations, regionSize should run from [0..L-1]
allowing one to calculate the mutual information at a fixed temperature 
"""

import subprocess as sp
import os

def main():
    folder_name = 'D%1.6f/L%03d/R%03d'
    sizes = [6,8,10,12,14,16,18,20,24,28,30,32]
    Dmin = 1.88
    Dmax = 1.91
    N_D = 30
    Dvals = [Dmin + (Dmax-Dmin)*i*1./N_D for i in xrange(N_D+1)]
    jnum = 0
    jnum2 = 0
    jmin = 1000
    jmax = 2000
    for D in Dvals:
        for S in sizes:
            for R in range(S):
                cur_path = folder_name % (D,S,R)
                if jnum >= jmin and jnum < jmax and S != 14 and S != 28:
                    try:
                        os.makedirs(cur_path)
                    except:
                        pass
                    sp.call('cp ../../../transfer .', cwd=cur_path, shell=True)
                    sp.call('cp ../../../betas.dat .', cwd=cur_path, shell=True)
                    sp.call('sqsub --mpp 2500M -r 168h -e run%d.err -o run%d.out ./transfer %d %d %f' % (jnum, jnum, S, R, D), cwd=cur_path, shell=True)
                jnum += 1
                if S != 14 and S != 28:
                    jnum2 += 1

if __name__=="__main__":
    main()
