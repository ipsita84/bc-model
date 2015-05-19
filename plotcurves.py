import matplotlib
#matplotlib.use('agg')

from matplotlib import pyplot as plt
import numpy as np

def main():
    fnames = 'S2_%03d.dat'
    sizes = [8,16,32,64]
    for S in sizes:
        d = np.loadtxt(fnames % S)
        x = [i*1.0/S for i in range(S+1)]
        plt.errorbar(x, d[:,0]/S, d[:,1]/S)
    plt.show()

if __name__ == "__main__":
    main()
