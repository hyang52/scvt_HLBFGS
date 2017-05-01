#!/bin/bash
#SBATCH -J scvt_fg_1           # job name
#SBATCH -o scvt_fg_1.o%j       # output and error file name (%j expands to jobID)
#SBATCH -n 32              # total number of mpi tasks requested
#SBATCH -p normal # queue (partition) -- normal, development, etc.
#SBATCH -t 01:00:00        # run time (hh:mm:ss) - 1.5 hours
#SBATCH --mail-user=lcgraham@fsu.edu
#SBATCH --mail-type=begin  # email me when the job starts
#SBATCH --mail-type=end    # email me when the job finishes
#SBATCH -A TG-EAR170003    # scvt 


module r scvt
cp parameters_fg_1 parameters
ibrun ./scvt_HLBFGS.exe           # run the MPI executable named a.out
