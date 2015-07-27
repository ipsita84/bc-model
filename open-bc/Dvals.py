"""
Correct temps for Blume Capel
T_c/J = 0.609
D/J = 1.965

We will measure the ratio trick at T_c and 2T_c
or Bc and Bc/2
"""

def main():
    N = 8
    Dmin = 1.7
    Dmax = 2.5
    Dvals = [(Dmax-Dmin)*i*1./N + Dmin for i in xrange(N+1)]
    Dmin = 1.8
    Dmax = 2.0
    Dvals += [(Dmax-Dmin)*i*1./N + Dmin for i in xrange(N+1)]
    Dvals += [1.965]
    Dvals.sort()
    fout = open('dvals.dat','w')
    for i in xrange(len(Dvals)):
        if i>1 and (Dvals[i] - Dvals[i-1]) < 1e-4:
            continue
        fout.write('%1.8f\n' % Dvals[i])
    fout.close()

if __name__ == "__main__":
    main()
