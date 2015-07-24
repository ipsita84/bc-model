"""
Correct temps for Blume Capel
T_c/J = 0.609
D/J = 1.965

We will measure the ratio trick at T_c and 2T_c
or Bc and Bc/2
"""

def main():
    Bc = 1./0.78
    N = 20
    betas = [1./N*Bc*i for i in xrange(1,N+1)]
    meas = [0 if (2*i % N) != 0 else 1 for i in xrange(1,N+1)]
    betas = [(i,j) for i,j in zip(betas,meas)]
    betas.append((1./0.77,1))
    betas.append((0.5/0.77,1))
    betas.append((1./0.79,1))
    betas.append((0.5/0.79,1))
    betas.append((1./0.80,1))
    betas.append((0.5/0.80,1))
    betas.sort(key = lambda X:X[0])
    fout = open('betas.dat','w')
    for i in betas:
        fout.write('%1.8f %d\n' % (i[0], i[1]))
    fout.close()

if __name__ == "__main__":
    main()
