#!/bin/bash
#SBATCH --job-name sems1      # Set a name for your job. This is especially useful if you have multiple jobs queued.
#SBATCH --partition short     # Slurm partition to use
#SBATCH --ntasks 1          # Number of tasks to run. By default, one CPU core will be allocated per task
#SBATCH --time 0-10:00        # Wall time limit in D-HH:MM
#SBATCH --mem-per-cpu=500     # Memory limit for each tasks (in MB)
#SBATCH -o out/ivmbatch_%j.out    # File to which STDOUT will be written
#SBATCH -e err/ivmbatch_%j.err    # File to which STDERR will be written
#SBATCH --mail-type=ALL       # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=dag89@uw.edu # Email to which notifications will be sent
#SBATCH --array=1-20

R CMD BATCH 1.R
