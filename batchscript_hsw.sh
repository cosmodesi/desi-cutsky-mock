#!/bin/bash
#SBATCH -n 1		# Number of tasks
#SBATCH -J qso  # Name of the job
#SBATCH -q regular
#SBATCH -C haswell
#SBATCH -N 1          # number of nodes
#SBATCH -c 64          # number of cpus per tasks
#SBATCH --time=48:00:00
#SBATCH -o ./out.qso.out
#SBATCH -e ./err.qso.err

source activate desilightcone

srun python main_qso.py
