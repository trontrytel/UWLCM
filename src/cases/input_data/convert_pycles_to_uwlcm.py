"""
Convert from PyCLES grid to UWLCM grid
"""

import os
import h5py
import hickle as hkl
import pickle as pkl
import numpy as np

#fname_hdf = "tmp/1.hdf"
fname_pkl = "tmp/10000001.pkl"

#f_hdf = h5py.File(fname_hdf, 'r')
f_pkl = open(fname_pkl, 'rb')

data = pkl.load(f_pkl)
#print("Keys: %s" % data.keys())
#Keys: ['v', 'u', 'temperature', 'w', 'qv']

# Convert horizontal velocity
pycles_v = np.array(data['v'])
x_dim    = pycles_v.shape[0]
z_dim    = pycles_v.shape[1]
uwlcm_v  = np.zeros((x_dim+1, z_dim+1))
# linear interpolation from the two neighbour velocity points from pycles
for idx_i in range(1, x_dim+1):
    for idx_k in range(1, z_dim):
        uwlcm_v[idx_i, idx_k] = 0.5 * (pycles_v[idx_i-1, idx_k-1] + pycles_v[idx_i-1, idx_k])
# linear extrapolation at the bottom
for idx_i in range(1, x_dim+1):
    uwlcm_v[idx_i, 0] = 2 * pycles_v[idx_i-1, 0] - uwlcm_v[idx_i, 1]    
# linear extrapolation at the top
for idx_i in range(1, x_dim+1):
    uwlcm_v[idx_i, -1] = 2 * pycles_v[idx_i-1,-1] - uwlcm_v[idx_i, -2]    
# cyclic 
uwlcm_v[0, :] = uwlcm_v[-1, :] 

#Convert verticl velocity
pycles_w = np.array(data['w'])
x_dim    = pycles_v.shape[0]
z_dim    = pycles_v.shape[1]
uwlcm_w  = np.zeros((x_dim+1, z_dim+1))
# linear interpolation from the two neighbour velocity points from pycles
for idx_i in range(1, x_dim):
    for idx_k in range(1, z_dim+1):
        uwlcm_w[idx_i, idx_k] = 0.5 * (pycles_w[idx_i-1, idx_k-1] + pycles_w[idx_i, idx_k-1])
# linear extrapolation at the right
for idx_k in range(1, z_dim+1):
    uwlcm_w[-1, idx_k] = 2 * pycles_w[-1, idx_k-1] - uwlcm_w[-2, idx_k]    
# cyclic at the left
for idx_k in range(1, z_dim+1):
    uwlcm_w[0, idx_k] = uwlcm_w[-1, idx_k]    
# wall at the bottom
uwlcm_w[:, 0] = 0.

#Convert scalars (T, qv)
for scl in ['temperature', 'qv']:
    pycles_scl = np.array(data[scl])
    x_dim      = pycles_scl.shape[0]
    z_dim      = pycles_scl.shape[1]
    uwlcm_scl  = np.zeros((x_dim+1, z_dim+1))

    # linear interpolation from the four neighbours
    for idx_i in range(1, x_dim):
        for idx_k in range(1, z_dim):
             uwlcm_scl[idx_i, idx_k] = 0.25 * (pycles_scl[idx_i-1, idx_k-1] + pycles_scl[idx_i, idx_k-1] + 
                                               pycles_scl[idx_i-1, idx_k  ] + pycles_scl[idx_i, idx_k  ])
    # linear extrapolation at the bottom
    for idx_i in range(1, x_dim):
        uwlcm_scl[idx_i,0] = pycles_scl[idx_i-1, 0] + pycles_scl[idx_i,0] - uwlcm_scl[idx_i,1]
    # linear extrapolation at the top
    for idx_i in range(1, x_dim):
        uwlcm_scl[idx_i,-1] = pycles_scl[idx_i-1, -1] + pycles_scl[idx_i,-1] - uwlcm_scl[idx_i,-2]
    # linear extrapolation at the right
    for idx_k in range(1, z_dim):
        uwlcm_scl[-1, idx_k] = pycles_scl[-1, idx_k-1] + pycles_scl[-1,idx_k]  - uwlcm_scl[-2,idx_k]
    # even more linear extrapolation in the two corners
    uwlcm_scl[-1,-1] = 2*pycles_scl[-1,-1] - uwlcm_scl[-2,-2]
    uwlcm_scl[-1, 0] = 2*pycles_scl[-1, 0] - uwlcm_scl[-2, 1]
    # cyclic at the left
    uwlcm_scl[0,:] = uwlcm_scl[-1,:]

    if (scl == 'temperature'):
        uwlcm_T = np.zeros((x_dim+1, z_dim+1))
        uwlcm_T = np.copy(uwlcm_scl) 
    if (scl == 'qv'):
        uwlcm_qv = np.zeros((x_dim+1, z_dim+1))
        uwlcm_qv = np.copy(uwlcm_scl) 

#Convert density TODO - for initial profiles

# dump data to hdf
#hkl.dump(data, fdir+'../../slices/'+str(it)+'.hdf', mode='w')
#it+=1
# Get the data
#data = list(f[a_group_key])
