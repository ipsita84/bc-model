import numpy as np
import matplotlib.pyplot as plt
import glob

def main():
    files = glob.glob("S2*")
    filedict = {}
    Lset = []
    for i in files:
        n = i.split('_')
        L = int(n[1][1:])
        if L not in Lset:
            Lset.append(L)
        beta = float(n[2][1:])
        filedict[(L,beta)] = i
    Lset.sort()
    for L in Lset:
        Bvals = []
        for i in filedict:
            if i[0] == L:
                if i[1] not in Bvals:
                    Bvals.append(i[1])
        Bvals.sort()
        # Now we have a sorted set of betas for the given system size, and we cna plot
        middata = []
        for beta in Bvals:
            data = np.loadtxt(filedict[(L,beta)])
            # Line L/2 contains the data we want
            MI = data[L/2,0] + data[L/2,0] - data[L,0]
            MIe = data[L,1]
            middata.append([MI,MIe])
        middata = np.array(middata)
        plt.errorbar(Bvals, middata[:,0]/L, middata[:,1]/L)
    plt.show()

if __name__ == "__main__":
    main()
