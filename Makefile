# Creates 4 binaries by compiling 3 different source files:
# normal: E-print-spin.cc
# replicaA: EmA-vs-beta.cc
# mutualinfo: Mutual-info-vs-beta.cc

# Example:
# make normal will compile E-vs-beta-Normal.cc and generate the
# executable "normal", which can then be run as ./normal and so on.

# Special commands:
# "make clean" will remove all generated binaries and .o object
# files;
# "make all" will compile all 3 source files and generate all
# 4 binaries in one step.

all: normal replicaA replicaB mutualinfo

CFLAGS = -I/home/ipsita/usr/include
LDFLAGS = -L/home/ipsita/usr/lib

normal: E-print-spin.cc
	icpc -Wall ${CFLAGS} ${LDFLAGS} -O3 E-print-spin.cc -o normal
	
replicaA: E-replicaA-l-print-spin.cc
	icpc -Wall ${CFLAGS} ${LDFLAGS} -O3 E-replicaA-l-print-spin.cc -o replicaA
	
replicaB: E-replicaB-l-print-spin.cc
	icpc -Wall ${CFLAGS} ${LDFLAGS} -O3 E-replicaB-l-print-spin.cc -o replicaB

transfer: Transfer-method.cc
	icpc -Wall ${CFLAGS} ${LDFLAGS} -O3 $< -o $@

mutualinfo: MI-partition-l.cc
	icpc -Wall -O3 ${CFLAGS} ${LDFLAGS} \
	`pkg-config --cflags --libs gsl tabdatrw-0.4 interp2dpp` \
	MI-partition-l.cc -o mutualinfo
	
.PHONY: clean

clean:
	rm -f normal replicaA replicaB mutualinfo transfer *.o 
