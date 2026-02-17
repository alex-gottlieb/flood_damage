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
#SBATCH --time=36:00:00
# Email address
#SBATCH --mail-user=Alexander.R.Gottlieb.GR@dartmouth.edu
# Email notifications (comma-separated options: BEGIN,END,FAIL)
#SBATCH --mail-type=BEGIN,END,FAIL
# CMIG account
#SBATCH --account=CMIG
#SBATCH --array=0-9

mv *.out /dartfs-hpc/rc/lab/C/CMIG/damages/slurm_out/

source /optnfs/common/miniconda3/etc/profile.d/conda.sh
conda activate arg24

# Run your Python script
python attribution_bootstrap.py noflood $SLURM_ARRAY_TASK_ID 
python attribution_bootstrap.py both $SLURM_ARRAY_TASK_ID 
python attribution_bootstrap.py ppt $SLURM_ARRAY_TASK_ID 
python attribution_bootstrap.py tws $SLURM_ARRAY_TASK_ID 
