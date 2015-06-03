import glob
#from matplotlib import pyplot

def main():
    files = glob.glob('S2*.dat')
    print files

if __name__ == "__main__":
    main()
