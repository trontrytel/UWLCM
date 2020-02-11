import os

outdir = "/home/ajaruga/clones/UWLCM/build/out_test_blk_1m_truth"

case = "dycoms_rf02"
nx = "129"
ny = "0"
nz = "301"

dt = "0.5"
outfreq = "10"
spinup = "900"
nt = "1800"

#dt = "0.5"
#outfreq = "1800"
#spinup = "7200"
#nt = "50400"

#nt = "25200"
#spinup = "3600"
#outfreq = "900"

micro = "blk_1m"
backend = "serial"

prs_tol = "6e-7"
rng_seed = "42"

cmd = "OMP_NUM_THREADS=32 /home/ajaruga/clones/UWLCM/build/src/bicycles"+\
      " --outdir="+outdir+" --case="+case+\
      " --nx="+nx+" --ny=0 --nz="+nz+" --dt="+dt+" --spinup="+spinup+\
      " --nt="+nt+" --micro="+micro+" --outfreq="+outfreq+\
      " --backend="+backend+" --rng_seed="+rng_seed+" --prs_tol="+prs_tol

print "running " + cmd
os.system(cmd)
