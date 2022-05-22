#!/bin/bash
#SBATCH -n 1	       # Number of tasks
#SBATCH -J ELGPHPH   # Name of the job
#SBATCH -q regular
#SBATCH -C haswell
#SBATCH -N 1           # number of nodes
#SBATCH -c 64          # number of cpus per tasks
#SBATCH --time=01:00:00
#SBATCH -o ./output/out.ELGPHPH.out
#SBATCH -e ./output/err.ELGPHPH.err

source activate desilightcone

srun python main.py  --galaxy ELG --phase PHPH
