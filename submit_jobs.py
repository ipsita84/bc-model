"""
Script for submitted a batch of jobs for the Blume-Capel model
Main job is run with the command

./transfer beta L regionSize

beta - inverse temperature
L - system size
regionSize - number of strips in region A

for a complete simulations, regionSize should run from [0..L-1]
allowing one to calculate the mutual information at a fixed temperature 

kT/J=0.695, D/J=1.965 at critical point
"""

import subprocess as sp
import os

def main():
    beta = 1.0/0.695
    folder_name = 'L%03d'
    sizes = [8,16,32,64]
    #sizes = [64]
    jnum = 0
    jmin = 0
    jmax = 400
    for S in sizes:
        cur_path = folder_name % S
        try:
            os.makedirs(cur_path)
        except:
            pass
        sp.call('cp ../transfer .', cwd=cur_path, shell=True)
        for R in range(S):
            if jnum >= jmin and jnum < jmax:
                sp.call('sqsub --mpp 500M -r 6000m -e run%d.err -o run%d.out ./transfer %f %d %d' % (jnum, jnum, beta, S, R), cwd=cur_path, shell=True)
            jnum += 1

if __name__=="__main__":
    main()
