#!/bin/bash

# Name of the job
#SBATCH --job-name=dmg_func
# Number of compute nodes
#SBATCH --nodes=1
# Number of tasks per node
#SBATCH --ntasks-per-node=1
# Number of CPUs per task
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4GB
# Walltime (job duration)
#SBATCH --time=8:00:00
# Email address
#SBATCH --mail-user=Alexander.R.Gottlieb.GR@dartmouth.edu
# Email notifications (comma-separated options: BEGIN,END,FAIL)
#SBATCH --mail-type=BEGIN,END,FAIL
# CMIG account
#SBATCH --partition=preemptable
#SBATCH --account=CMIG
#SBATCH --array=0-12

mv *.out /dartfs-hpc/rc/lab/C/CMIG/damages/slurm_out/

source /optnfs/common/miniconda3/etc/profile.d/conda.sh
conda activate arg24

PARAMS=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" product_combos.txt)
set -- $PARAMS  

# Run your Python script
python damage_func_bootstrap.py "$1" "$2" county
python damage_func_bootstrap.py "$1" "$2" state
python damage_func_bootstrap.py "$1" "$2" year

