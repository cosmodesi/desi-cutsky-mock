#!/bin/bash
#SBATCH -N 1
#SBATCH -C cpu
#SBATCH --cpus-per-task=64
#SBATCH -t 10:00:00
#SBATCH --ntasks=1
#SBATCH -q regular
#SBATCH -A desi
#SBATCH --array=1-24
#SBATCH -J ic_array
#SBATCH -o logs/ic_%A_%a.out
#SBATCH -e logs/ic_%A_%a.err

source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
mkdir -p logs

srun --cpus-per-task=$SLURM_CPUS_PER_TASK --ntasks=1 python main_ic.py "$SLURM_ARRAY_TASK_ID"
