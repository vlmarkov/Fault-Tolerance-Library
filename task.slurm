#!/bin/bash

#SBATCH --nodes=1 --ntasks-per-node=8
#SBATCH --job-name=mpi_ft

cd $SLURM_SUBMIT_DIR

mpiexec ./primes