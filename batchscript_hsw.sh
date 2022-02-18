#!/bin/bash
#SBATCH -n 1		# Number of tasks
#SBATCH -J LRG998  # Name of the job
#SBATCH -q regular
#SBATCH -C haswell
#SBATCH -N 1          # number of nodes
#SBATCH -c 64          # number of cpus per tasks
#SBATCH --time=1:30:00
#SBATCH -o ./output/out.test_lrg998.out
#SBATCH -e ./output/err.test_lrg998.err

source activate desilightcone

srun python main.py --galaxy LRG --phase 998
