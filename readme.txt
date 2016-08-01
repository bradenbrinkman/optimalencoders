readme:

Code listed in this folder is most up to date code as of 05/20/2015.

Folder should contain:

Simulation Code:
ppie_solver4b.c
functions_method7.h

Hyak files:
ppie_solver4_generatecommandsforhyak.c
Makefile


To compile, gsl v1.15 or greater is needed. Compile using command

gcc ppie_solver4b.c -o ppie_solver4 -lgsl -lgslcblas -lm

On hyak, compile the file by running the Makefile.

The program is set up to take in, in order, the stim variance (sigstim), upstream variance (sigin), downstream variance (sigout), upstream correlation coeff (vin) , downstream correlation coeff (vout), kappa, an integer seed that is passed to the random number generator, and a 0 or 1 to use an ON-ON initial guess (0) or an ON-OFF initial guess (1).

e.g., to run the program

./ppie_solver3 0.8 1.0 1.0 -0.64 0.0 0.0 123 1

runs the program with parameter values sigstim = 0.8, sigin = 1.0, sigout = 1.0, vin = -0.64, vout = 0.0, kappa = 0.0, a seed of 123, and uses an initial guess in which the two nonlinearities have opposite slopes. These parameters also set a benchmark case
for which the exact solution is known. The exact value for |D| is 0.16157863018260993. Program should give De and Di of 0.1615 or 0.1616 (tolerance of 1e-4), with either sign.

The program assumes it is in a folder called executables, which is contained in a parent folder that contains a separate folder called ‘data’. The program prints files into this data folder. i.e., the folder structure is

<parent>/executables (contains ppie_solver4b)
<parent>/data (data files printed to this folder)

The file ppie_solver4_generatecommandsforhyak.c, once compiled and run, generates a file called ’sweepfile’ which contains parameter sets to be passed to ppie_solver4 to be run on hyak.

To run the jobs on hyak, we use the parallel_sql database. (see here: https://sig.washington.edu/itsigs/Hyak_parallel-sql )

This submits the jobs listed in sweepfile to a database from which the tasks are pulled. To submit multiple jobs, run, e.g., 

for job in $(seq 1 20); do qsub -q bf parallel_sql_job.sh; done

This will submit 20 jobs to the backfill queue. Each job can pull tasks from the database until it is empty. Jobs which do not finish in the 4 hour time limit are killed and sent back to the list of available tasks. For this reason, the code has been set up to produce a file which records the current iteration and the decoding coefficients D1 and D2. If the task gets killed, upon being initialized again it will read in the previous completed iteration and D1 and D2 estimates to pick up where it was left off when killed.
