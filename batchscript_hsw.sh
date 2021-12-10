#!/bin/bash
#SBATCH -n 1		# Number of tasks
#SBATCH -J ab_elg0  # Name of the job
#SBATCH -q regular
#SBATCH -C haswell
#SBATCH -N 1          # number of nodes
#SBATCH -c 64          # number of cpus per tasks
#SBATCH --time=0:45:00
#SBATCH -o ./output/out.ab_elg0.out
#SBATCH -e ./output/err.ab_elg0.err

source activate desilightcone

srun python main_elg.py --phase 0
