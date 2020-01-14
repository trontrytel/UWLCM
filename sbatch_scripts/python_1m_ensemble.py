import os
import sys
import pickle
import numpy as np

# slurm array id to tell the jobs apart
slurm_id = os.environ['SLURM_ARRAY_TASK_ID']
print "slurm_id = " + slurm_id

# read in the input parameters
with open("enki_class_dump.pkl") as pkl_file:
    eki = pickle.load(pkl_file)
r_c0_arr = eki.u[-1]
print r_c0_arr

# construct the Dycoms run command
case = "dycoms_rf02"
nx = "129"
ny = "0"
nz = "301"
dt = "1"
spinup = "10"
nt = "20"
micro = "blk_1m"
outfreq = "10"
backend = "serial"

r_c0 = str(r_c0_arr[int(slurm_id) - 1][0])
outdir = "/home/ajaruga/clones/UWLCM/build/out_test_blk_1m_rc0_"+r_c0

print " "
print " inside run enseble saving to: ", outdir

cmd = "OMP_NUM_THREADS=1 /home/ajaruga/clones/UWLCM/build/src/bicycles"+\
      " --outdir="+outdir+" --case="+case+\
      " --nx="+nx+" --ny=0 --nz="+nz+" --dt="+dt+" --spinup="+spinup+\
      " --nt="+nt+" --micro="+micro+" --outfreq="+outfreq+\
      " --backend="+backend+" --r_c0="+r_c0+" --rng_seed=42"

print "running " + cmd
sys.stdout.flush()

# run model
os.system(cmd)
