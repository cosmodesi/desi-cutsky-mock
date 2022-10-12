#!/bin/bash
#SBATCH -n 1	       # Number of tasks
#SBATCH -J LRG4322   # Name of the job
#SBATCH -q regular
#SBATCH -C haswell
#SBATCH -N 1           # number of nodes
#SBATCH -c 64          # number of cpus per tasks
#SBATCH --time=00:30:00
#SBATCH -o ./output/out.LRG4322.out
#SBATCH -e ./output/err.LRG4322.err

source activate desilightcone
date
# time srun python main.py  --galaxy LRG --phase 4322
time srun python main.py  --phase 4322
date
