#!/usr/bin/python3
# vim: set ts=4 sw=4 tw=80:

import glob
import os.path
import numpy as nm
import sys

if len(sys.argv) != 2:
	print("\033[1;31mError\033[0m"
	      + ": Expecting exactly one (integer) input, got "
	      + str(len(sys.argv)-1) + '.')
	sys.exit(-1);

idx = int(sys.argv[1])

# Gets the name of all files in current dir matching 'I2*.dat'
# and stores as a list of strings in variable lof
lof = glob.glob("./I2*.dat")

val1 = []
val2 = []
files = {}

for name in lof:
#	Strip the leading './' from filenames indicating current dir
	name = name.strip('./')
	
#	Store the filename as a string with its extension ('.dat') removed
	name2   = os.path.splitext(name)[0]

#	Split the name2 string at every occurrence of '-'
	splstr = name2.split('-')

#	val1 is the array of all IDX1, without duplication
	if int(splstr[1]) not in val1:
		val1.append(int(splstr[1]))

#	val2 stores all the [IDX1, IDX2] pairs
	val2.append([int(splstr[1]), int(splstr[2])])
	
#	Dictionary to get filenames given (IDX1,IDX2) as a tuple
	files[(int(splstr[1]), int(splstr[2]))] = name

val1 = nm.array(val1)
val2 = nm.array(val2)

resarr = [];

print("# FILE IDX1 IDX2 IDX1/IDX2 Col0 Col1 Str.Val Sub.Value")
for matidx in val1:
	hffile = ''
	strval = 0;
# This for loop scans over val2 to search for IDX2 = IDX1/2
	for i in val2:
		if i[0] == matidx:
			if (i[1] * 2) == i[0]:
				hffile = files[tuple(i)]
				strval  = nm.loadtxt(hffile)[idx, 1]

#	Loop again over val2 to compute the subtracted value for each matching file
	for i in val2:
		if i[0] == matidx:
			fdata = nm.loadtxt(files[tuple(i)])

#			Print formatted data to terminal
			print(files[tuple(i)], format(i[0], ' 5d'),
			      format(i[1], '5d'),
			      format(i[1] * 1.0/i[0],' 12.4E'),
			      format(fdata[idx,0], '12.4E'),
			      format(fdata[idx,1], '12.4E'),
			      format(strval, '12.4E'),
			      format(fdata[idx,1]-strval, '12.4E')
			     )
			resarr.append([i[0], i[1], i[1] * 1.0/i[0], fdata[idx,0],
			               fdata[idx,1], strval, fdata[idx,1]-strval])

res = nm.reshape(resarr,(-1, 7))
# Write data to file using numpy.savetxt
nm.savetxt("SCANDATA.dat", res,
           fmt='%3i %3i % 15.4E % 15.4E % 15.4E % 15.4E % 15.4E',
           header="I1\tI2\tI2/I1\tCOL0\tCOL1\tSTRVAL\tSUBVAL")
