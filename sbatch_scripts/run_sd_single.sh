#!/bin/bash

# Submit this script with: sbatch run_array_central.sh

#SBATCH --job-name=l_base_17        # job name
#SBATCH --time=09:30:00             # time limit

#SBATCH --nodes=1                   # Use one node
#SBATCH --ntasks=1                  # Number of CPUs
#SBATCH --mem-per-cpu=1G
#SBATCH --gres=gpu:1                # request 1 GPU

#SBATCH --mail-user=ajaruga@caltech.edu
#SBATCH --mail-type=ALL

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module load singularity/3.2.1

BASE=/home/ajaruga/
CONTAINER=$BASE/singularities/sng_ubuntu_18_04_cuda_10_new.sif

srun singularity -d exec --nv $CONTAINER python python_lgrngn.py
