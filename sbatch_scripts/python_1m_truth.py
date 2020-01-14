import os

outdir = "/home/ajaruga/clones/UWLCM/build/out_test_blk_1m_truth_long"

case = "dycoms_rf02"
nx = "129"
ny = "0"
nz = "301"
dt = "1"
spinup = "3600"
nt = "25200"
micro = "blk_1m"
outfreq = "900"
backend = "serial"


cmd = "OMP_NUM_THREADS=32 /home/ajaruga/clones/UWLCM/build/src/bicycles"+\
      " --outdir="+outdir+" --case="+case+\
      " --nx="+nx+" --ny=0 --nz="+nz+" --dt="+dt+" --spinup="+spinup+\
      " --nt="+nt+" --micro="+micro+" --outfreq="+outfreq+\
      " --backend="+backend+" --rng_seed=42"

print "running " + cmd
os.system(cmd)
