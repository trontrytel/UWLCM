#!/bin/bash

# Submit this script with: sbatch run_array_central.sh

#SBATCH --job-name=eki_init_job_1m           # job name
#SBATCH --time=00:10:00                      # time limit

#SBATCH --nodes=1                            # Use one node
#SBATCH --ntasks=1                           # Number of CPUs
#SBATCH --ntasks-per-node=1

#SBATCH --mail-user=ajaruga@caltech.edu
#SBATCH --mail-type=ALL

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module load singularity/3.2.1

BASE=/home/ajaruga/
CONTAINER=$BASE/singularities/sng_ubuntu_18_04_cuda_10_new.sif

echo "Starting at `date`"
echo "Running on hosts: $SLURM_NODELIST"
echo "Running on hosts: $SLURM_JOB_NODELIST"
echo "Running on $SLURM_NNODES nodes."
echo "Running on $SLURM_NPROCS processors."
echo "Current working directory is `pwd`"

srun singularity -d exec  $CONTAINER python python_1m_eki_setup.py
