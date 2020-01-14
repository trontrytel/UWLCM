import numpy as np
import os

def gen_wet_bins(left, right, nbins):
    """
    Generate bin ranges for the histogram plot of droplets and rain drops
    input:
    - left - leftmost edge
    - right - rightmost edge
    - nbins - number of bins
    """
    bin_str = ""
    moms = [0,]
    moms_str = ','.join(str(x) for x in moms)
    str_fmt = "{:.2e}"
    left_edges = np.geomspace(left, right, num=nbins, endpoint=False, dtype=float)
    for i in np.arange(nbins):
        bin_str += str_fmt.format(left_edges[i])
        bin_str += ":"
        if i < nbins-1:
            bin_str += str_fmt.format(left_edges[i+1])
            bin_str += "|"+moms_str+";"
        else:
            bin_str += str_fmt.format(right)
            bin_str += "|"+moms_str
    return bin_str

case = "dycoms_rf02"
nx = "129"
ny = "0"
nz = "301"
dt = "1"
nt = "25200"
spinup = "3600"
outfreq = "900"
backend = "CUDA"

#slurm_id = os.environ['SLURM_ARRAY_TASK_ID']
#rng_seed = str(int(slurm_id) - 1 + 14)
#outdir = "out_test_lgrngn_"+rng_seed

rng_seed = "48"
outdir = "/home/ajaruga/clones/UWLCM/build/out_test_lgrngn_hangup_"+rng_seed

micro = "lgrngn"
sd_conc = "512"
sstp_cond = "10"
sstp_coal = "10"

left  = 25e-7
right = 2500e-6
nbins = 100
bins_wet_str = gen_wet_bins(left, right, nbins)

#      " --turb_coal=1 --sgs=true --gccn=1"+\

cmd = "OMP_NUM_THREADS=1 /home/ajaruga/clones/UWLCM/build/src/bicycles"+\
      " --outdir="+outdir+" --case="+case+\
      " --nx="+nx+" --ny=0 --nz="+nz+" --dt="+dt+" --spinup="+spinup+\
      " --nt="+nt+" --micro="+micro+" --outfreq="+outfreq+\
      " --backend="+backend+" --rng_seed="+rng_seed+" --sd_conc="+sd_conc+\
      " --sstp_cond="+sstp_cond+" --sstp_coal="+sstp_coal+\
      " --gccn=1"+\
      " --out_wet=\"25e-7:25e-6|0,1,2,3,6;25e-6:25e-4|0,1,2,3,6;25e-8:25e-4|0,1,2,3,6;"+bins_wet_str+"\""

print "running " + cmd
os.system(cmd)
