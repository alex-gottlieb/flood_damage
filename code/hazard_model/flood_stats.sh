#!/bin/bash

# Name of the job
#SBATCH --job-name=flood_stats
# Number of compute nodes
#SBATCH --nodes=1
# Number of tasks per node
#SBATCH --ntasks-per-node=1
# Number of CPUs per task
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4GB
# Walltime (job duration)
#SBATCH --time=2:00:00
# Email address
#SBATCH --mail-user=Alexander.R.Gottlieb.GR@dartmouth.edu
# Email notifications (comma-separated options: BEGIN,END,FAIL)
#SBATCH --mail-type=BEGIN,END,FAIL
# CMIG account
#SBATCH --account=CMIG
#SBATCH --partition=preemptable
#SBATCH --array=18-35
mv *.out /dartfs-hpc/rc/lab/C/CMIG/damages/slurm_out/

source /optnfs/common/miniconda3/etc/profile.d/conda.sh
conda activate arg24

PARAMS=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" product_combos.txt)
set -- $PARAMS  

# Run your Python script
python3 flood_stats.py "$1" "$2" 0 county
python3 flood_stats.py "$1" "$2" 0 huc4
