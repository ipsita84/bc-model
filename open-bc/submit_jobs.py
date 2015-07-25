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
    folder_name = 'L%03d/R%03d'
    #sizes = [8,16,32,64]
    #sizes = [32]
    sizes = [8,12,16,20,24]
    Dval = 1.965
    jnum = 0
    jmin = 0
    jmax = 800
    for S in sizes:
        for R in range(S):
            cur_path = folder_name % (S,R)
            if jnum >= jmin and jnum < jmax:
                try:
                    os.makedirs(cur_path)
                except:
                    pass
                sp.call('cp ../../transfer .', cwd=cur_path, shell=True)
                sp.call('cp ../../betas.dat .', cwd=cur_path, shell=True)
                sp.call('sqsub --mpp 500M -r 96h -e run%d.err -o run%d.out ./transfer %d %d %f' % (jnum, jnum, S, R, Dval), cwd=cur_path, shell=True)
            jnum += 1

if __name__=="__main__":
    main()
