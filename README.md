Instructions for running the code on Sharcnet/orca

Make the executable with the command
    make transfer
this produces the executable named "transfer".

Copy this file to the directory "test"-- all other files are already there in the directory "test".

To run the code, type
    python submit_jobs.py
this will run the code, at the temperatures listed in "betas.dat".

Measurements of the Renyi entanglemtne entropy will be done for every temperature which has a "1" as the second column.

The "sizes = [...]" line in submit_jobs.py can be changed to affect the full lattice size that is being simulated.
It must also be changed later in analyze.py, but the plotting code will find the appropriate files automatically.

fit_all2.py scipt details

linFit() creates a linear extrapolation of the exact shape function from the file dedekindeta_07.dat
We then load the appropriate sizes needed for each cut from file into memory.
For each cut, we doa  linear fitting in 1/L (see def func) and we store the y-intercept and the uncertainty in the y-intercept
From this, we compare the extrapolated values vs. the exact values, weighting the different by 1/error.
In principle, we could also weight the difference by it's absolute magnitude.
This gives us the chi^2 for each set of parameters (D,beta)
Using this, we can find the minimum chi^2 for each D value, as a function of beta.
We can use the same set of data to find the minimum coefficient "c" for each (D,beta) pair as well.
