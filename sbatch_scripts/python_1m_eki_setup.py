import pickle
import h5py
import numpy as np
import enki as enki

number_of_ensembles = 10

parameters = np.zeros([number_of_ensembles, 1])
parameters[:,0] = np.linspace(1e-6, 8e-4, number_of_ensembles)

fname = "/home/ajaruga/clones/UWLCM/build/out_test_blk_1m_truth/timestep0000001800.h5"
my_file = h5py.File(fname, 'r')
rr_profile = np.average(my_file["rr"], axis=0)

truth = rr_profile
cov = np.zeros([rr_profile.shape[0], rr_profile.shape[0]])
np.fill_diagonal(cov, 0.1)
cov = cov * truth

eki = enki.EKI(parameters, truth, cov)

with open("enki_class_dump.pkl", "w") as pkl_file:
    pickle.dump(eki, pkl_file)
