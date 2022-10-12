#!/bin/bash
#SBATCH -n 1	       # Number of tasks
#SBATCH -J GALGALPHPH   # Name of the job
#SBATCH -q regular
#SBATCH -C haswell
#SBATCH -N 1           # number of nodes
#SBATCH -c 64          # number of cpus per tasks
#SBATCH --time=00:30:00
#SBATCH -o ./output/out.GALGALPHPH.out
#SBATCH -e ./output/err.GALGALPHPH.err

source activate desilightcone
date
# time srun python main.py  --galaxy GALGAL --phase PHPH
time srun python main.py  --phase PHPH
date
