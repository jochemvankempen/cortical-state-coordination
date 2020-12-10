#!/bin/bash
# Example SLURM job script for job arrays
#
# Tell SLURM which project's account to use:
#
#SBATCH -A atvc
#
# This job runs a command called myprog and takes input from a set of files called 
# nnn.input, where nnn is the array task IDs
#
# set up a job array with tasks numbered 1 to 100.
#SBATCH --array=1-83
# 
# give the array a single job name
#SBATCH -J analyses
#
# Standard output is saved in this file
#SBATCH -o analyses_%A_%a.out
#
# Standard error messages are saved in this file
#SBATCH -e analyses_%A_%a.err
#
# request n cores per task
#SBATCH -c 10
#
# request a run-time
#SBATCH -t 0-30:00:00
#
# SLURM defaults to the directory you were working in when you submitted the job.
# Output files are also put in this directory. To set a different working directory add:
#
# SBATCH –-workdir=/mnt/nfs/home/njv20/repositories/clustercomputing/outputFiles/
#
# Tell SLURM if you want to be emailed when your job starts, ends, etc.
# Currently mail can only be sent to addresses @ncl.ac.uk
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jochem.van-kempen@ncl.ac.uk
#

module load MATLAB
cd /mnt/nfs/home/njv20/repositories/cortical-state-coordination/analyses
matlab -nodisplay -r "batch_analyses"
