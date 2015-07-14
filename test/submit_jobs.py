"""
Script for submitted a batch of jobs for the Blume-Capel model
Main job is run with the command

./transfer L regionSize

L - system size
regionSize - number of strips in region A

for a complete simulations, regionSize should run from [0..L-1]
allowing one to calculate the mutual information at a fixed temperature 
"""

import subprocess as sp
import os

def main():
    folder_name = 'L%03d/R%03d'
    #sizes = [8,16,32,64]
    #sizes = [32]
    sizes = [8,16]
    jnum = 0
    jmin = 0
    jmax = 400
    for S in sizes:
        for R in range(S):
            cur_path = folder_name % (S,R)
            try:
                os.makedirs(cur_path)
            except:
                pass
            sp.call('cp ../../transfer .', cwd=cur_path, shell=True)
            sp.call('cp ../../betas.dat .', cwd=cur_path, shell=True)
            if jnum >= jmin and jnum < jmax:
                sp.call('sqsub --mpp 500M -r 96h -e run%d.err -o run%d.out ./transfer %d %d' % (jnum, jnum, S, R), cwd=cur_path, shell=True)
            jnum += 1

if __name__=="__main__":
    main()
