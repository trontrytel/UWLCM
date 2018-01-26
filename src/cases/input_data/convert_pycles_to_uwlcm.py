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


def make_non_divergent(U_comp, W_comp, dx, dz, eps, fname):
    """ Input:  v and w components of uwlcm velocity field
        Output: v and w velocity components such that d/dxv + d/dzw = 0
        Helmholtz decomposition using Jacobi method and cyclic horizontal boundary condition """

    x_dim = np.shape(U_comp)[0]
    z_dim = np.shape(U_comp)[1]
       
    # calculate the diveregence of u (our RHS)
    div = div_h_cyc(U_comp, W_comp, dx, dz)

    # first guess for the scalar potential
    scl_pot     = np.zeros((x_dim, z_dim))
    # boundary condition
    scl_pot[ :,  0] = W_comp[ :,  0] * dz # bottom
    scl_pot[ :, -1] = W_comp[ :, -1] * dz # top

    err=44
    it = 0
    # iterative search for scalar potential
    while (err >= eps):
        for idx_i in range(0, x_dim):
            for idx_k in range(1, z_dim-1):
                if (idx_i == 0):
                    scl_pot[0, idx_k] = \
                       ( \
                           (scl_pot[-1, idx_k  ] + scl_pot[1, idx_k  ]) * dz**2 + \
                           (scl_pot[ 0, idx_k-1] + scl_pot[0, idx_k+1]) * dx**2 - \
                           dx**2 * dz**2 * div[0, idx_k] \
                       ) / 2. / (dx**2 + dz**2)
                elif (idx_i > 0 and idx_i < x_dim - 1):
                    scl_pot[idx_i, idx_k] = \
                        ( \
                            (scl_pot[idx_i-1, idx_k  ] + scl_pot[idx_i+1, idx_k  ]) * dz**2 + \
                            (scl_pot[idx_i,   idx_k-1] + scl_pot[idx_i,   idx_k+1]) * dx**2 - \
                            dx**2 * dz**2 * div[idx_i, idx_k] \
                        ) / 2. / (dx**2 + dz**2)
                elif (idx_i == x_dim - 1):
                    scl_pot[idx_i, idx_k] = \
                        ( \
                            (scl_pot[-2,   idx_k  ] + scl_pot[0,    idx_k  ]) * dz**2 + \
                            (scl_pot[idx_i,idx_k-1] + scl_pot[idx_i,idx_k+1]) * dx**2 - \
                            dx**2 * dz**2 * div[idx_i, idx_k] \
                        ) / 2. / (dx**2 + dz**2)
                else:
                    print idx_i
                    assert(False)

        # calc u_div = grad(scalar_pot)
        u_div = np.gradient(scl_pot, dx, dz)
        u_div[0][0,:]  = (scl_pot[1,:] - scl_pot[-2,:]) / 2. / dx
        u_div[0][-1,:] = u_div[0][0,:]
        
        # calc u_rot = u - u_div
        divless_v = U_comp - u_div[0]
        divless_w = W_comp - u_div[1]
   
        div_check = div_h_cyc(divless_v, divless_w, dx, dz)
        err = np.sum(np.abs(div_check)) / x_dim / z_dim
        print "iter = ", it, "  L1(error)/nx/nz = ", err
        it +=1

    min_v = min(np.min(divless_v), np.min(uwlcm_v))
    max_v = max(np.max(divless_v), np.max(uwlcm_v))
    min_w = min(np.min(divless_w), np.min(uwlcm_w))
    max_w = max(np.max(divless_w), np.max(uwlcm_w))
    zoom = -1
    plt.figure(figsize=(15,10))
    plt.subplot(2,3,1)
    plt.imshow(div_check[:,:zoom], aspect='auto')
    plt.title('div_check')
    plt.colorbar()
    plt.subplot(2,3,4)
    plt.imshow(scl_pot[:,:zoom], vmin=np.min(scl_pot), vmax=np.max(scl_pot), aspect='auto')
    plt.title('scl_pot')
    plt.colorbar()
    plt.subplot(2,3,2)
    plt.imshow(divless_v[:,:zoom], vmin=min_v, vmax=max_v, aspect='auto')
    plt.title('divless_v')
    plt.colorbar()
    plt.subplot(2,3,5)
    plt.imshow(U_comp[:,:zoom], vmin=min_v, vmax=max_v, aspect='auto')
    plt.title('uwlcm_v')
    plt.colorbar()
    plt.subplot(2,3,3)
    plt.imshow(divless_w[:,:zoom], vmin=min_w, vmax=max_w, aspect='auto')
    plt.title('divless_w')
    plt.colorbar()
    plt.subplot(2,3,6)
    plt.imshow(W_comp[:,:zoom], vmin=min_w, vmax=max_w, aspect='auto')
    plt.title('uwlcm_w')
    plt.colorbar()
    plt.tight_layout()
    plt.savefig(fname+"_"+str(it)+".png")
    #plt.show()
    plt.close()

    return divless_v, divless_w


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

# create the output hdf5 file
output_hdf = h5py.File("dycoms/dycoms_velocity.h5", 'w')
v = output_hdf.create_group(u'v')
v.attrs[u'units'] = u'm/s'
w = output_hdf.create_group(u'w')
w.attrs[u'units'] = u'm/s'
v_nodiv = output_hdf.create_group(u'v_nodiv')
v_nodiv.attrs[u'units'] = u'm/s'
w_nodiv = output_hdf.create_group(u'w_nodiv')
w_nodiv.attrs[u'units'] = u'm/s'

# read velocity from all timesteps
it  = 0
dx  = 35.     # TODO - read in the dx and dz from file
dz  = 5.
eps = 4.4e-4  # http://culture.pl/en/article/mickiewicz-unravelled-a-little-known-fact-about-polands-best-known-bard

for fname in os.listdir(folder_pkl):
    if fname.endswith(".pkl"):
        print "reading data from ", fname
        f_pkl   = open(folder_pkl+fname, 'rb')
        data    = pkl.load(f_pkl)

        uwlcm_v = read_horizontal_velocity(data, x_dim, z_dim)
        uwlcm_w = read_vertical_velocity(data, x_dim, z_dim)

        uwlcm_v_nodiv, uwlcm_w_nodiv = make_non_divergent(uwlcm_v, uwlcm_w, dx, dz, eps, fname)

        v.create_dataset(str(it), data=uwlcm_v)
        w.create_dataset(str(it), data=uwlcm_w)
        v_nodiv.create_dataset(str(it), data=uwlcm_v_nodiv)
        w_nodiv.create_dataset(str(it), data=uwlcm_w_nodiv)

        it+=1

output_hdf.close()


