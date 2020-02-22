import h5py
import pickle
import os
import numpy as np
import enki as enki

with open("enki_class_dump.pkl") as pkl_file:
    eki = pickle.load(pkl_file)

with open("enki_cloud_range.pkl") as pkl_file:
    cloud_range = pickle.load(pkl_file)

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
    print "reading from: ", outdir

    my_file_1 = h5py.File(outdir + '/timestep0000050400.h5', 'r')
    my_file_2 = h5py.File(outdir + '/timestep0000048600.h5', 'r')
    my_file_3 = h5py.File(outdir + '/timestep0000046800.h5', 'r')
    my_file_4 = h5py.File(outdir + '/timestep0000045000.h5', 'r')

    rr_profile = np.average(my_file_1["rr"], axis=0) + np.average(my_file_2["rr"], axis=0) +\
                 np.average(my_file_3["rr"], axis=0) + np.average(my_file_4["rr"], axis=0)

    rc_profile = np.average(my_file_1["rc"], axis=0) + np.average(my_file_2["rc"], axis=0) +\
                 np.average(my_file_3["rc"], axis=0) + np.average(my_file_4["rc"], axis=0)


    rc_profile_slice = rc_profile[cloud_range[0] : cloud_range[1]]
    rr_profile_slice = rc_profile[0 : cloud_range[1]]

    g_arr[param_id, :] = np.append(rr_profile_slice, rc_profile_slice)

    # remove output files
    os.system("rm -r "+outdir)

eki.EnKI(g_arr)

print "parameters   : ", eki.u
print "error        : ", eki.error
print "converged    : ", eki.converged
#print "observations : ", eki.g

with open("enki_class_dump.pkl", "w") as pkl_file:
    pickle.dump(eki, pkl_file)
