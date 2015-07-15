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
