"""
Script for submitted a batch of jobs for the Blume-Capel model
Main job is run with the command

./transfer bmin bmax dbeta L regionSize

beta - inverse temperature
L - system size
regionSize - number of strips in region A

for a complete simulations, regionSize should run from [0..L-1]
allowing one to calculate the mutual information at a fixed temperature 

kT/J=0.609, D/J=1.965 at critical point
"""

import subprocess as sp
import os

def main():
    bc = 1.0/0.609
    bmin = bc/10.
    bmax = 2*bc+bc/20.
    dbeta = bc/10.
    folder_name = 'L%03d/R%03d'
    sizes = [8,16,32,64]
    #sizes = [64]
    jnum = 0
    jmin = 2
    jmax = 400
    for S in sizes:
        for R in range(S):
            cur_path = folder_name % (S,R)
            try:
                os.makedirs(cur_path)
            except:
                pass
            if jnum >= jmin and jnum < jmax:
                sp.call('cp ../../transfer .', cwd=cur_path, shell=True)
                sp.call('sqsub --mpp 500M -r 2d -e run%d.err -o run%d.out ./transfer %f %f %f %d %d' % (jnum, jnum, bmin, bmax, dbeta, S, R), cwd=cur_path, shell=True)
            jnum += 1

if __name__=="__main__":
    main()
