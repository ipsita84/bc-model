"""
Correct temps for Blume Capel
T_c/J = to be found
D/J = to be found

We will measure the ratio trick at T_c and 2T_c
or Bc and Bc/2
"""

def main():
    Bc = 1./0.9
    N = 20
    betas = [1./N*Bc*i for i in xrange(1,N+1)]
    meas = [0 for i in xrange(1,N+1)]
    betas = [(i,j) for i,j in zip(betas,meas)]
    N = 100
    minB = 1./0.85
    maxB = 1./0.75
    mb = [minB + (maxB - minB)*i*1./N for i in range(N+1)]
    meas = [1 for i in range(N+1)]
    mb += [(minB + (maxB - minB)*i*1./N)/2. for i in range(N+1)]
    meas += [1 for i in range(N+1)]
    betas += [(i,j) for i,j in zip(mb,meas)]
    betas.sort(key = lambda X:X[0])
    fout = open('betas.dat','w')
    for i in betas:
        fout.write('%1.8f %d\n' % (i[0], i[1]))
    fout.close()

if __name__ == "__main__":
    main()
