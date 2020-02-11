import h5py
import pickle
import os
import numpy as np
import enki as enki

with open("enki_class_dump.pkl") as pkl_file:
    eki = pickle.load(pkl_file)

# read in the input parameters
r_c0 = eki.u[-1]
print " "
print " inside after run enseble "
print "r_c0 : ", r_c0

# read data from all enseble runs
params_length = r_c0.shape[0]
obs_length    = eki.g_t.shape[0]
print params_length, obs_length
g_arr = np.zeros([params_length, obs_length])

for param_id in range(params_length):

    outdir = "/home/ajaruga/clones/UWLCM/build/out_test_blk_1m_rc0_"+str(r_c0[param_id][0])
    print "reding from: ", outdir

    my_file = h5py.File(outdir + '/timestep0000001800.h5', 'r')
    rr_data = np.array(my_file["rr"])

    # vertical profile of rr
    rr_profile = np.average(rr_data, axis=0)
    print rr_profile

    g_arr[param_id, :] = rr_profile

    # remove output files
    os.system("rm -r "+outdir)

eki.EnKI(g_arr)

print "parameters   : ", eki.u
print "observations : ", eki.g

with open("enki_class_dump.pkl", "w") as pkl_file:
    pickle.dump(eki, pkl_file)
