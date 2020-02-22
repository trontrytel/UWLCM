#!/bin/bash

# Submit this script with: sbatch file_name

#SBATCH --job-name=array_job_1m     # job name
#SBATCH --time=00:30:00             # time limit

#SBATCH --nodes=1                   # Use one node
#SBATCH --ntasks=1                  # Number of CPUs
#SBATCH --cpus-per-task=32          # number of OpenMP threads

#SBATCH --mail-user=ajaruga@caltech.edu
#SBATCH --mail-type=ALL

#SBATCH --output=array_%A-%a.out    # Standard output and error log
#SBATCH --array=1-10                # Array range

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module load singularity/3.2.1

BASE=/home/ajaruga/
CONTAINER=$BASE/singularities/sng_ubuntu_18_04_cuda_10_new.sif

srun singularity -d exec  $CONTAINER python python_1m_ensemble.py
