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

import scipy
from scipy import ndimage

import matplotlib.pyplot as plt

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
    # linear extrapolation at the bottom and at the top
    for idx_i in range(1, x_dim+1):
        uwlcm_v[idx_i, 0] = 2 * pycles_v[idx_i-1, 0] - uwlcm_v[idx_i, 1]    
        uwlcm_v[idx_i,-1] = 2 * pycles_v[idx_i-1,-1] - uwlcm_v[idx_i,-2]    
    # cyclic
    uwlcm_v[0, :] = uwlcm_v[-1, :] 
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
    # linear extrapolation at the right (cyclic b-cond)
    for idx_k in range(1, z_dim+1):
        uwlcm_w[-1, idx_k] = 0.5 * (pycles_w[-1, idx_k-1] - pycles_w[0, idx_k-1])
    # cyclic at the left
    uwlcm_w[0, :] = uwlcm_w[-1, :]    
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
    #TODO - add cyclic assumption
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

def div_h_cyc(U_comp, W_comp, dx, dz): 
    """ Calculate divergence assuming cyclic boundary condition"""
    tmp_grad_x = np.gradient(U_comp, dx)[0]
    tmp_grad_z = np.gradient(W_comp, dz)[1]
    # cyclic horizontal b-cond
    tmp_grad_x[0,:] = (U_comp[1, :] - U_comp[-2, :]) / 2. / dx
    tmp_grad_x[-1,:] = tmp_grad_x[0,:]
    return tmp_grad_x + tmp_grad_z


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

## velocity
#uwlcm_v0    = read_horizontal_velocity(data, x_dim, z_dim)
#uwlcm_w0    = read_vertical_velocity(data, x_dim, z_dim)
## water vapor mixing ratio
#uwlcm_qv0   = read_scalars(data, x_dim, z_dim, 'qv')
## Convert qv to rv 
#uwlcm_rv0   = uwlcm_qv0 / (1. - uwlcm_qv0)
## dry air density
#uwlcm_p0    = read_profiles(fname_stats, x_dim, z_dim, 'p0_full')
#uwlcm_T0    = read_scalars(data, x_dim, z_dim, 'temperature')
#uwlcm_pd0   = uwlcm_p0 * (1. - uwlcm_rv0 / (uwlcm_rv0 + cm.eps))
#uwlcm_rhod0 = uwlcm_pd0 / cm.R_d / uwlcm_T0
## dry air potential temperature
#p_1000      = 100000.
#uwlcm_thd0  = uwlcm_T0 * (p_1000 / uwlcm_pd0)**(cm.R_d/cm.c_pd)
#
## TODO - process velocity field to make it non-divergent
#
## save t=0 data
#with h5py.File('dycoms/dycoms_init.h5', 'w') as init_hdf:
#    init_hdf.create_dataset("uwlcm_rv0",   data=uwlcm_rv0)
#    init_hdf.create_dataset("uwlcm_v0",    data=uwlcm_v0)
#    init_hdf.create_dataset("uwlcm_w0",    data=uwlcm_w0)
#    init_hdf.create_dataset("uwlcm_thd0",  data=uwlcm_thd0)
#    init_hdf.create_dataset("uwlcm_rhod0", data=uwlcm_rhod0)

# ------------ t=1,...,n --------------

## create the output hdf5 file
#output_hdf = h5py.File("dycoms/dycoms_velocity.h5", 'w')
#v = output_hdf.create_group(u'v')
#v.attrs[u'units'] = u'm/s'
#w = output_hdf.create_group(u'w')
#w.attrs[u'units'] = u'm/s'
#
## read velocity from all timesteps
#it = 0
#for fname in os.listdir(folder_pkl):
#    if fname.endswith(".pkl"):
#        print "reading data from ", fname
#        f_pkl   = open(folder_pkl+fname, 'rb')
#        data    = pkl.load(f_pkl)
#
#        uwlcm_v = read_horizontal_velocity(data, x_dim, z_dim)
#        uwlcm_w = read_vertical_velocity(data, x_dim, z_dim)
#
#        v.create_dataset(str(it), data=uwlcm_v)
#        w.create_dataset(str(it), data=uwlcm_w)
#
#        it+=1
#
#output_hdf.close()

# ---------- non-divergent velocity --------------------
# TODO
#for fname in os.listdir(folder_pkl):
#    if fname.endswith(".pkl"):

fname = "10007200.pkl"
print "reading data from ", fname
f_pkl   = open(folder_pkl+fname, 'rb')
data    = pkl.load(f_pkl)

uwlcm_v = read_horizontal_velocity(data, x_dim, z_dim)
uwlcm_w = read_vertical_velocity(data,   x_dim, z_dim)

# TODO - read in the dx and dz from file
dx=35.   #35
dz=5.    #5
# calculate the diveregence of u (our RHS)
div = div_h_cyc(uwlcm_v, uwlcm_w, dx, dz)

err=44.
it = 0
# first guess for the scalar potential
scl_pot     = np.zeros((x_dim+1, z_dim+1))
# boundary condition
scl_pot[ :,  0] = uwlcm_w[ :,  0] * dz # bottom
scl_pot[ :, -1] = uwlcm_w[ :, -1] * dz # top
scl_pot[ 0,  :] = uwlcm_v[ 0,  :] * dx # left
scl_pot[-1,  :] = uwlcm_v[-1,  :] * dx # right

# iterative search for scalar potential
#while (err >= 1e-7):
for it in range(1000):
    for idx_i in range(0, x_dim+1):
        for idx_k in range(1, z_dim):
            if (idx_i ==0):
                scl_pot[0, idx_k] = \
                   ( \
                       (scl_pot[-1, idx_k  ] + scl_pot[1 , idx_k  ]) * dz**2 + \
                       (scl_pot[idx_i, idx_k-1] + scl_pot[idx_i, idx_k+1]) * dx**2 - \
                       dx**2 * dz**2 * div[idx_i, idx_k] \
                   ) / 2. / (dx**2 + dz**2)
            elif (idx_i > 0 and idx_i < x_dim):
                scl_pot[idx_i, idx_k] = \
                    ( \
                        (scl_pot[idx_i-1, idx_k  ] + scl_pot[idx_i+1, idx_k  ]) * dz**2 + \
                        (scl_pot[idx_i,   idx_k-1] + scl_pot[idx_i,   idx_k+1]) * dx**2 - \
                        dx**2 * dz**2 * div[idx_i, idx_k] \
                    ) / 2. / (dx**2 + dz**2)
            elif (idx_i == x_dim):
                scl_pot[idx_i, idx_k] = \
                    ( \
                        (scl_pot[-2,   idx_k  ] + scl_pot[0,    idx_k  ]) * dz**2 + \
                        (scl_pot[idx_i,idx_k-1] + scl_pot[idx_i,idx_k+1]) * dx**2 - \
                        dx**2 * dz**2 * div[idx_i, idx_k] \
                    ) / 2. / (dx**2 + dz**2)
            else:
                assert(false)

    # calc u_div = grad(scalar_pot)
    u_div = np.gradient(scl_pot, dx, dz)
    u_div[0][0,:]  = (scl_pot[1,:] - scl_pot[-2,:]) / 2. / dx
    u_div[0][-1,:] = u_div[0][0,:]
    
    # calc u_rot = u - u_div
    divless_v = uwlcm_v - u_div[0]
    divless_w = uwlcm_w - u_div[1]
   
    div_check = div_h_cyc(divless_v, divless_w, dx, dz)

    min_v = min(np.min(divless_v), np.min(uwlcm_v))
    max_v = max(np.max(divless_v), np.max(uwlcm_v))
    min_w = min(np.min(divless_w), np.min(uwlcm_w))
    max_w = max(np.max(divless_w), np.max(uwlcm_w))

    zoom = -1

    plt.figure(figsize=(15,10))
    plt.subplot(2,3,1)
    #plt.imshow(scipy.ndimage.rotate(div_check[:,:zoom],90), vmin=np.min(div_check), vmax=np.max(div_check), aspect='auto')
    plt.imshow(div_check[:,:zoom], vmin=-0.01, vmax=0.01, aspect='auto')
    plt.title('div_check')
    plt.colorbar()
    plt.subplot(2,3,4)
    #plt.imshow(scipy.ndimage.rotate(scl_pot[:,:zoom], 90), vmin=np.min(scl_pot), vmax=np.max(scl_pot), aspect='auto')
    plt.imshow(scl_pot[:,:zoom], vmin=np.min(scl_pot), vmax=np.max(scl_pot), aspect='auto')
    plt.title('scl_pot')
    plt.colorbar()
    plt.subplot(2,3,2)
    #plt.imshow(scipy.ndimage.rotate(divless_v[:,:zoom], 90), vmin=min_v, vmax=max_v, aspect='auto')
    plt.imshow(divless_v[:,:zoom], vmin=min_v, vmax=max_v, aspect='auto')
    plt.title('divless_v')
    plt.colorbar()
    plt.subplot(2,3,5)
    #plt.imshow(scipy.ndimage.rotate(uwlcm_v[0:10,:zoom], 90), vmin=min_v, vmax=max_v, aspect='auto')
    plt.imshow(uwlcm_v[:,:zoom], vmin=min_v, vmax=max_v, aspect='auto')
    plt.title('uwlcm_v')
    plt.colorbar()
    plt.subplot(2,3,3)
    #plt.imshow(scipy.ndimage.rotate(divless_w[:,:zoom], 90), vmin=min_w, vmax=max_w, aspect='auto')
    plt.imshow(divless_w[:,:zoom], vmin=min_w, vmax=max_w, aspect='auto')
    plt.title('divless_w')
    plt.colorbar()
    plt.subplot(2,3,6)
    #plt.imshow(scipy.ndimage.rotate(uwlcm_w[:,:zoom], 90), vmin=min_w, vmax=max_w, aspect='auto')
    plt.imshow(uwlcm_w[:,:zoom], vmin=min_w, vmax=max_w, aspect='auto')
    plt.title('uwlcm_w')
    plt.colorbar()
    plt.tight_layout()
    plt.savefig(str(it)+".png")
    it +=1
    plt.close()

    #err = np.max(np.abs(div_check))
    err = np.sum(np.abs(div_check))
    print err

# save the new velocity

# TODO -add plotting pycles vs uwlcm ws divless version
# TODO add ploting the divergence field
# TODO at what error should we stop
# TODO what norm to calculate error?


#plt.subplot(1,2,1)
#plt.imshow(uwlcm_v[:,:], vmin=min_v, vmax=max_v, aspect='auto')
#plt.title('uwlcm_v')
#plt.subplot(1,2,2)
#plt.imshow(divless_v[:,:], vmin=min_v, vmax=max_v, aspect='auto')
#plt.title('divless_v')
#plt.colorbar()
#plt.show()
#
#plt.subplot(1,2,1)
#plt.imshow(uwlcm_w[:,:], vmin=min_w, vmax=max_w, aspect='auto')
#plt.title('uwlcm_w')
#plt.subplot(1,2,2)
#plt.imshow(divless_w[:,:], vmin=min_w, vmax=max_w, aspect='auto')
#plt.title('divless_w')
#plt.colorbar()
#plt.show()

#min_v = min(np.min(divless_v), np.min(uwlcm_v))
#max_v = max(np.max(divless_v), np.max(uwlcm_v))
#min_w = min(np.min(divless_w), np.min(uwlcm_w))
#max_w = max(np.max(divless_w), np.max(uwlcm_w))
# 
#plt.subplot(1,1,1)
#plt.imshow(div_check[:,:], vmin=np.min(div_check), vmax=np.max(div_check), aspect='auto')
#plt.title('div_check')
#plt.colorbar()
#plt.show()
