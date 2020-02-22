import pickle
import h5py
import numpy as np
import enki as enki

# truth: 5e-4

number_of_ensembles = 10

parameters = np.zeros([number_of_ensembles, 1])
parameters[:,0] = np.linspace(1e-6, 8e-4, number_of_ensembles)

# 1 hour average
fname_1 = "/home/ajaruga/clones/UWLCM/build/out_test_blk_1m_truth/timestep0000050400.h5"
fname_2 = "/home/ajaruga/clones/UWLCM/build/out_test_blk_1m_truth/timestep0000048600.h5"
fname_3 = "/home/ajaruga/clones/UWLCM/build/out_test_blk_1m_truth/timestep0000046800.h5"
fname_4 = "/home/ajaruga/clones/UWLCM/build/out_test_blk_1m_truth/timestep0000045000.h5"
my_file_1 = h5py.File(fname_1, 'r')
my_file_2 = h5py.File(fname_2, 'r')
my_file_3 = h5py.File(fname_3, 'r')
my_file_4 = h5py.File(fname_4, 'r')

# rain and cloud profiles
rr_profile = np.average(my_file_1["rr"], axis=0) + np.average(my_file_2["rr"], axis=0) +\
             np.average(my_file_3["rr"], axis=0) + np.average(my_file_4["rr"], axis=0)

rc_profile = np.average(my_file_1["rc"], axis=0) + np.average(my_file_2["rc"], axis=0) +\
             np.average(my_file_3["rc"], axis=0) + np.average(my_file_4["rc"], axis=0)

# cloud base and cloud top indices
cb_idx = np.min(np.where(rc_profile > 1e-8))
ct_idx = np.max(np.where(rc_profile > 1e-8))
ext_idx = 5

# slice through rain from cloud top to the ground
# and through cloud from cloud top to cloud base
# +/- 5 vertical levels
rr_profile_slice = rr_profile[0 : ct_idx + ext_idx]
rc_profile_slice = rc_profile[cb_idx - ext_idx : ct_idx + ext_idx]

# covariance cannot be zero
# for places where the truth profiles are zero assign a small number
idx_zero_rr = np.where(rr_profile_slice == 0)
idx_zero_rc = np.where(rc_profile_slice == 0)
rr_profile_slice[idx_zero_rr] = 1e-11
rc_profile_slice[idx_zero_rc] = 1e-11

# construct covariance matrix
truth = np.append(rr_profile_slice, rc_profile_slice)
cov = np.zeros([truth.shape[0], truth.shape[0]])
np.fill_diagonal(cov, 0.1)
cov = cov * truth

eki = enki.EKI(parameters, truth, cov)

with open("enki_class_dump.pkl", "w") as pkl_file:
    pickle.dump(eki, pkl_file)

with open("enki_cloud_range.pkl", "w") as pkl_file:
    pickle.dump([cb_idx - ext_idx, ct_idx + ext_idx], pkl_file)
