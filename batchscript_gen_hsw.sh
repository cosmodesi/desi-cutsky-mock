#!/bin/bash
#SBATCH -n 1	       # Number of tasks
#SBATCH -J orig_1   # Name of the job
#SBATCH -q regular
#SBATCH -C haswell
#SBATCH -N 1           # number of nodes
#SBATCH -c 64          # number of cpus per tasks
#SBATCH --time=02:00:00
#SBATCH -o ./output/out.orig_1.out
#SBATCH -e ./output/err.orig_1.err

source activate desilightcone

srun python main.py
