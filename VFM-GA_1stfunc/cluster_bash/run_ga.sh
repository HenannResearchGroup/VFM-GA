#!/bin/bash

# Specify a job name:
#SBATCH -J run
#SBATCH --time=48:00:00 # time for EACH simulation in the array
#SBATCH --ntasks=2
#SBATCH --cpu-freq=high
#SBATCH -N 1
#SBATCH --partition=batch

# Run a command
module load intel
./run_ga.exe

