#!/bin/bash

# Submit this script with: sbatch this_script_name.sh

#SBATCH --job-name=array_of_sds     # job name
#SBATCH --time=10:00:00             # time limit

#SBATCH --nodes=1                   # Use one node
#SBATCH --ntasks=1                  # Number of CPUs
#SBATCH --mem-per-cpu=1G
#SBATCH --gres=gpu:1                # request 1 GPU

#SBATCH --mail-user=ajaruga@caltech.edu
#SBATCH --mail-type=ALL

#SBATCH --output=array_%A-%a.out    # Standard output and error log
#SBATCH --array=1-5                 # Array range

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module load singularity/3.2.1

BASE=/home/ajaruga/
CONTAINER=$BASE/singularities/sng_ubuntu_18_04_cuda_10_new.sif

srun singularity -d exec --nv $CONTAINER python python_lgrngn.py
