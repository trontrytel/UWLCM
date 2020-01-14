#!/bin/bash

# Submit this script with: sbatch file_name

#SBATCH --job-name=array_job_1m     # job name
#SBATCH --time=00:10:00             # time limit

#SBATCH --nodes=1                   # Use one node
#SBATCH --ntasks=1                  # Number of CPUs
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=1G
#SBATCH --partition=any

#SBATCH --mail-user=ajaruga@caltech.edu
#SBATCH --mail-type=ALL

#SBATCH --output=array_%A-%a.out    # Standard output and error log
#SBATCH --array=1-2                 # Array range

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module load singularity/3.2.1

BASE=/home/ajaruga/
CONTAINER=$BASE/singularities/sng_ubuntu_18_04_cuda_10_new.sif

srun singularity -d exec  $CONTAINER python python_1m_ensemble.py
