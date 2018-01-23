"""
Convert from PyCLES grid to UWLCM grid
  - initial condition
  - velocity for each time step
"""
import os
import h5py
import hickle as hkl
import pickle as pkl
import numpy as np

import sys
sys.path.append('/USers/ajaruga/clones/libcloudphxx/build/bindings/python/')
from libcloudphxx import common as cm

def read_horizontal_velocity(data, x_dim, z_dim):
    """ Convert horizontal velocity from PyCLES to UWLCM grid """
    pycles_v = np.array(data['v'])
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
    #import matplotlib.pyplot as plt
    #plt.subplot(1,2,1)
    #plt.imshow(pycles_v[:,:], vmin=-0.0005, vmax=0.0005, aspect='auto')
    #plt.title('pycles_v')
    #plt.subplot(1,2,2)
    #plt.imshow(uwlcm_v[:,:],  vmin=-0.0005, vmax=0.0005, aspect='auto')
    #plt.title('uwlcm_v')
    #plt.colorbar()
    #plt.show()
    return uwlcm_v

def read_vertical_velocity(data, x_dim, z_dim):
    """ Convert vertical velocity from PyCLES to UWLCM grid """
    pycles_w = np.array(data['w'])
    uwlcm_w  = np.zeros((x_dim+1, z_dim+1))
    # linear interpolation from the two neighbour velocity points from pycles
    #    TODO - interpolation smoothes the velocity field, think of something else?
    for idx_i in range(1, x_dim):
        for idx_k in range(1, z_dim+1):
            uwlcm_w[idx_i, idx_k] = 0.5 * (pycles_w[idx_i-1, idx_k-1] + pycles_w[idx_i, idx_k-1])
    # linear extrapolation at the right and left
    for idx_k in range(1, z_dim+1):
        uwlcm_w[-1, idx_k] = 0.5 * (pycles_w[-1, idx_k-1] + pycles_w[-1, idx_k-1] - uwlcm_w[-2, idx_k])
        uwlcm_w[ 0, idx_k] = 0.5 * (pycles_w[ 0, idx_k-1] + pycles_w[ 0, idx_k-1] - uwlcm_w[ 1, idx_k])
    # cyclic at the left
    #for idx_k in range(1, z_dim+1):
    #    uwlcm_w[0, idx_k] = uwlcm_w[-1, idx_k]    
    # wall at the bottom
    uwlcm_w[:, 0] = 0.
    return uwlcm_w
    
def read_scalars(data, x_dim, z_dim, scalar):
    """ Convert scalars (T, qv) from PyCLES to UWLCM grid """
    pycles_scl = np.array(data[scalar])
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
    return uwlcm_scl   
 
def read_profiles(fname_stats, x_dim, z_dim, profile):
    """ Convert initial profiles of density and pressure from PyCLES to 2D UWLCM """
    f_hdf      = h5py.File(fname_stats, 'r')
    pycles_prf = np.array(f_hdf['reference'][profile])
    uwlcm_prf  = np.zeros((x_dim+1, z_dim+1))
    uwlcm_prf[:, 1:] = pycles_prf
    uwlcm_prf[:, 0]  = 2 * uwlcm_prf[:, 1] - uwlcm_prf[:,2]
    return uwlcm_prf   



# file with initial profiles
#fname_stats = "/Users/ajaruga/clones/UWLCM/src/cases/input_data/dycoms/Output.DYCOMS_RF01.04b3a/stats/Stats.DYCOMS_RF01.nc"
fname_stats = "dycoms/Output.DYCOMS_RF01.04b3a/stats/Stats.DYCOMS_RF01.nc"
# file with simulation data from t=0
fname_pkl   = "dycoms/Output.DYCOMS_RF01.04b3a/Visualization/10000000.pkl"
# folder with simulation data
folder_pkl = "dycoms/Output.DYCOMS_RF01.04b3a/Visualization/"

# ------------ initial condition --------------

print "reading initial condition from ", fname_pkl 

# get simulation dimensions
f_pkl      = open(fname_pkl, 'rb')
data       = pkl.load(f_pkl)
pycles_v   = np.array(data['v'])
x_dim      = pycles_v.shape[0]
z_dim      = pycles_v.shape[1]

# velocity
uwlcm_v0    = read_horizontal_velocity(data, x_dim, z_dim)
uwlcm_w0    = read_vertical_velocity(data, x_dim, z_dim)
# water vapor mixing ratio
uwlcm_qv0   = read_scalars(data, x_dim, z_dim, 'qv')
# Convert qv to rv 
uwlcm_rv0   = uwlcm_qv0 / (1. - uwlcm_qv0)
# dry air density
uwlcm_p0    = read_profiles(fname_stats, x_dim, z_dim, 'p0_full')
uwlcm_T0    = read_scalars(data, x_dim, z_dim, 'temperature')
uwlcm_pd0   = uwlcm_p0 * (1. - uwlcm_rv0 / (uwlcm_rv0 + cm.eps))
uwlcm_rhod0 = uwlcm_pd0 / cm.R_d / uwlcm_T0
# dry air potential temperature
p_1000      = 100000.
uwlcm_thd0  = uwlcm_T0 * (p_1000 / uwlcm_pd0)**(cm.R_d/cm.c_pd)

# TODO - process velocity field to make it non-divergent

# save t=0 data
with h5py.File('dycoms/dycoms_init.h5', 'w') as init_hdf:
    init_hdf.create_dataset("uwlcm_rv0",   data=uwlcm_rv0)
    init_hdf.create_dataset("uwlcm_v0",    data=uwlcm_v0)
    init_hdf.create_dataset("uwlcm_w0",    data=uwlcm_w0)
    init_hdf.create_dataset("uwlcm_thd0",  data=uwlcm_thd0)
    init_hdf.create_dataset("uwlcm_rhod0", data=uwlcm_rhod0)

np.set_printoptions(threshold=np.inf)

#print "uwlcm_rv0 shape",    uwlcm_rv0.shape[0], " ", uwlcm_rv0.shape[1]
##print uwlcm_rv0 
#print "uwlcm_v0 shape",     uwlcm_v0.shape[0], " ", uwlcm_v0.shape[1]
##print uwlcm_v0
#print "uwlcm_w0 shape",     uwlcm_w0.shape[0], " ", uwlcm_w0.shape[1]
##print uwlcm_w0
#print "uwlcm_thd0 shape",   uwlcm_thd0.shape[0], " ", uwlcm_thd0.shape[1]
##print uwlcm_thd0
#print "uwlcm_rhod0 shape",  uwlcm_rhod0.shape[0], " ", uwlcm_rhod0.shape[1]
##print uwlcm_rhod0

# ------------ t=1,...,n --------------

# create the output hdf5 file
output_hdf = h5py.File("dycoms/dycoms_velocity.h5", 'w')
v = output_hdf.create_group(u'v')
v.attrs[u'units'] = u'm/s'
w = output_hdf.create_group(u'w')
w.attrs[u'units'] = u'm/s'

# read velocity from all timesteps
it = 0
for fname in os.listdir(folder_pkl):
    if fname.endswith(".pkl"):
        print "reading data from ", fname
        f_pkl   = open(folder_pkl+fname, 'rb')
        data    = pkl.load(f_pkl)

        # TODO - process velocity field to make it non-divergent

        uwlcm_v = read_horizontal_velocity(data, x_dim, z_dim)
        uwlcm_w = read_vertical_velocity(data, x_dim, z_dim)

        v.create_dataset(str(it), data=uwlcm_v)
        w.create_dataset(str(it), data=uwlcm_w)

        it+=1

output_hdf.close()
