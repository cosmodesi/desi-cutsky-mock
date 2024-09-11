#!/bin/bash
#SBATCH --time=1:30:00
#SBATCH --qos=shared
#SBATCH --cpus-per-task=64
#SBATCH --constraint=cpu
#SBATCH --array=1-24
#SBATCH --account=desi
#SBATCH --ntasks=1
#SBATCH --mem=220G

#source /global/common/software/desi/desi_environment.sh main
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
srun python main_ic_new.py --phase $SLURM_ARRAY_TASK_ID

