#!/bin/bash
#SBATCH -n 1		# Number of tasks
#SBATCH -J ez_qsoPHPH  # Name of the job
#SBATCH -q regular
#SBATCH -C haswell
#SBATCH -N 1          # number of nodes
#SBATCH -c 64          # number of cpus per tasks
#SBATCH --time=1:30:00
#SBATCH -o ./output/out.ez_qsoPHPH.out
#SBATCH -e ./output/err.ez_qsoPHPH.err

source activate desilightcone

srun python main_qso.py --phase PHPH
