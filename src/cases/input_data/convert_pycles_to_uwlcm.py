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

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import sys
sys.path.append('/home/ajaruga/clones/libcloudphxx/build/bindings/python')
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
    # TODO - add open/closed top bottom boundary
    return uwlcm_v

def read_vertical_velocity(data, x_dim, z_dim):
    """ Convert vertical velocity from PyCLES to UWLCM grid """
    pycles_w = np.array(data['w'])
    uwlcm_w  = np.zeros((x_dim+1, z_dim+1))
    # linear interpolation from the two neighbour velocity points from pycles
    for idx_i in range(1, x_dim):
        for idx_k in range(1, z_dim+1):
            uwlcm_w[idx_i, idx_k] = 0.5 * (pycles_w[idx_i-1, idx_k-1] + pycles_w[idx_i, idx_k-1])
    # linear interpolation at the right (cyclic b-cond)
    for idx_k in range(1, z_dim+1):
        uwlcm_w[-1, idx_k] = 0.5 * (pycles_w[-1, idx_k-1] + pycles_w[0, idx_k-1])
    # cyclic at the left
    uwlcm_w[0, :] = uwlcm_w[-1, :]    
    # wall at the bottom
    uwlcm_w[:, 0] = 0.
    # TODO - add open/closed top bottom boundary
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
    # linear interplation at the right (cyclic boundary condition)
    for idx_k in range(1, z_dim):
        uwlcm_scl[-1, idx_k] =  0.25 * (pycles_scl[-1, idx_k-1] + pycles_scl[0, idx_k-1] + 
                                        pycles_scl[-1, idx_k  ] + pycles_scl[0, idx_k  ])
    # linear extrapolation at the bottom
    for idx_i in range(1, x_dim):
        uwlcm_scl[idx_i,0] = pycles_scl[idx_i-1, 0] + pycles_scl[idx_i,0] - uwlcm_scl[idx_i,1]
    # linear extrapolation at the top
    for idx_i in range(1, x_dim):
        uwlcm_scl[idx_i,-1] = pycles_scl[idx_i-1, -1] + pycles_scl[idx_i,-1] - uwlcm_scl[idx_i,-2]
    # top and bottom right corners (cyclic boundary condition)
    uwlcm_scl[-1, 0] = pycles_scl[-1, 0] + pycles_scl[ 0, 0] - uwlcm_scl[-1, 1]
    uwlcm_scl[-1,-1] = pycles_scl[-1,-1] + pycles_scl[ 0,-1] - uwlcm_scl[-1,-2]
    # cyclic at the left
    uwlcm_scl[0, :] = uwlcm_scl[-1, :] 
    # TODO - add open/closed top bottom boundary
    return uwlcm_scl   
 
def read_profiles(fname_stats, x_dim, z_dim, profile):
    """ Convert initial profiles of density and pressure from PyCLES to 2D UWLCM """
    f_hdf      = h5py.File(fname_stats, 'r')
    pycles_prf = np.array(f_hdf['reference'][profile])
    uwlcm_prf  = np.zeros((x_dim+1, z_dim+1))
    uwlcm_prf[:, 1:] = pycles_prf
    uwlcm_prf[:, 0]  = 2 * uwlcm_prf[:, 1] - uwlcm_prf[:,2]
    # TODO - add open/closed top bottom boundary
    return uwlcm_prf   

def fill_halo_cyclic(field, halo):
    """ Fill cyclic left - right halo region """
    x_dim   = np.shape(field)[0]
    z_dim   = np.shape(field)[1]
    h_field = np.zeros((x_dim + 2*halo, z_dim))

    h_field[halo:-halo,:] = field

    for idx in range(halo):
        h_field[idx,       :] = h_field[-halo-halo+idx, :]
        h_field[-halo+idx, :] = h_field[halo+idx,       :]

    return h_field

def fill_halo_top(field, halo):
    """ Open boundary at the top """
    x_dim   = np.shape(field)[0]
    z_dim   = np.shape(field)[1]
    h_field = np.zeros((x_dim, z_dim+halo))

    h_field[:,0:-halo] = field

    for idx in range(halo):
        h_field[:,-halo+idx] = h_field[:, -halo-1]

    return h_field

def fill_halo_bottom_closed(field, halo):
    """ closed boundary at the bottom """

    x_dim   = np.shape(field)[0]
    z_dim   = np.shape(field)[1]
    h_field = np.zeros((x_dim, z_dim+halo))

    h_field[:, halo:] = field

    for idx in range(halo):
        h_field[:, halo-idx-1] = -1 * h_field[:, halo + idx + 1]

    return h_field

def fill_halo_bottom_open(field, halo):
    """ closed boundary at the bottom """

    x_dim   = np.shape(field)[0]
    z_dim   = np.shape(field)[1]
    h_field = np.zeros((x_dim, z_dim+halo))

    h_field[:, halo:] = field

    for idx in range(halo):
        h_field[:, halo-idx-1] = -1 * h_field[:, halo]

    return h_field

def grad_h_cyc(scl, dx, dz):
    """ Calculate gradient operator of a scalar assuming cyclic boundary condition"""
    x_dim  = np.shape(scl)[0]
    z_dim  = np.shape(scl)[1]
    grad_x = np.zeros((x_dim, z_dim))
    grad_z = np.zeros((x_dim, z_dim))
   
    for idx_i in range(1,x_dim-1):
        grad_x[idx_i,:] = (scl[idx_i+1,:] - scl[idx_i-1, :]) / 2. / dx
    grad_x[0,  :] = (scl[1,:] - scl[-1,:]) / 2. / dx
    grad_x[-1, :] = (scl[0,:] - scl[-2,:]) / 2. / dx

    for idx_k in range(1, z_dim-1):
        grad_z[:, idx_k] = (scl[:, idx_k+1] - scl[:, idx_k-1]) / 2. / dz
    grad_z[:,  0] = (scl[:,  1] - scl[:,  0]) / dz
    grad_z[:, -1] = (scl[:, -1] - scl[:, -2]) / dz

    return grad_x, grad_z

def grad_dual_cyc(scl, dx, dz):
    x_dim = np.shape(scl)[0]
    z_dim = np.shape(scl)[1]
    grad_x = np.zeros([x_dim+1, z_dim])
    grad_z = np.zeros([x_dim , z_dim+1])

    hlo = 2

    h_scl = fill_halo_cyclic(scl, hlo)
    h_scl = fill_halo_top(h_scl, hlo)
    h_scl = fill_halo_bottom_open(h_scl, hlo)

    for idx_i in range(0, x_dim+1):
        for idx_k in range(0, z_dim):
            grad_x[idx_i, idx_k] = (h_scl[idx_i+hlo, idx_k+hlo] - h_scl[idx_i+hlo-1, idx_k+hlo]) / dx
    for idx_k in range(0, z_dim+1):
        for idx_i in range(0,x_dim):
            grad_z[idx_i, idx_k] = (h_scl[idx_i+hlo, idx_k+hlo] - h_scl[idx_i+hlo, idx_k+hlo-1]) / dz

    return grad_x, grad_z

def div_h_cyc(U_comp, W_comp, dx, dz): 
    """ Calculate divergence of a vector assuming cyclic boundary condition"""
    # U_comp and W_comp passed here are one column shorter than the uwlcm velocity fields
    # This is to avoid doubling the one column that is cyclic (i.e. the left edge of the velocity field matrix
    #    is exactly equal to its right edge
    tmp_grad_x = np.gradient(U_comp, dx)[0]
    tmp_grad_z = np.gradient(W_comp, dz)[1]
    # cyclic horizontal b-cond
    tmp_grad_x[ 0,:] = (U_comp[1, :] - U_comp[-1, :]) / 2. / dx
    tmp_grad_x[-1,:] = (U_comp[0, :] - U_comp[-2, :]) / 2. / dx
    # closed top bottom boundary
    tmp_grad_z[:, 0] = 0.
    tmp_grad_x[:,-1] = 0.

    return tmp_grad_x + tmp_grad_z

def div_dual(dual_u, dual_w, dx, dz):

    x_dim = dual_w.shape[0]
    z_dim = dual_u.shape[1]

    div = np.zeros([x_dim, z_dim])

    for idx_i in range(0, x_dim):
        for idx_k in range(0, z_dim):
            div[idx_i, idx_k] = (dual_u[idx_i+1, idx_k] - dual_u[idx_i, idx_k]) / dx + (dual_w[idx_i, idx_k+1] - dual_w[idx_i, idx_k]) / dz

    return div

def lap_h_cyc(scl, dx, dz): 
    """ Calculate Laplace operator of a scalar assuming cyclic boundary condition"""
    x_dim = np.shape(scl)[0]
    z_dim = np.shape(scl)[1]
    lap   = np.zeros((x_dim, z_dim))

    hlo = 4
    h_scl = fill_halo_cyclic(scl, hlo)
    h_scl = fill_halo_top(h_scl, hlo)
    h_scl = fill_halo_bottom_open(h_scl, hlo)

    for idx_i in range(0, x_dim):
        for idx_k in range(0, z_dim):
            lap[idx_i, idx_k] = (h_scl[idx_i+1+hlo, idx_k+hlo  ] + h_scl[idx_i-1+hlo, idx_k+hlo  ] - 2*h_scl[idx_i+hlo, idx_k+hlo]) / dx**2 + \
                                (h_scl[idx_i+hlo  , idx_k+1+hlo] + h_scl[idx_i+hlo,   idx_k-1+hlo] - 2*h_scl[idx_i+hlo, idx_k+hlo]) / dz**2

    #for idx_i in range(1, x_dim-1):
    #    for idx_k in range(1, z_dim-1):
    #        lap[idx_i, idx_k] = (scl[idx_i+1, idx_k  ] + scl[idx_i-1, idx_k  ] - 2*scl[idx_i, idx_k]) / dx**2 + \
    #                            (scl[idx_i  , idx_k+1] + scl[idx_i,   idx_k-1] - 2*scl[idx_i, idx_k]) / dz**2
    #for idx_i in range(1, x_dim-1):
    #    lap[idx_i,  0] =  (scl[idx_i+1, 0] + scl[idx_i-1, 0] - 2*scl[idx_i, 0]) / dx**2 - 2*(scl[idx_i, 0]) / dx**2
    #    lap[idx_i, -1] =  (scl[idx_i+1,-1] + scl[idx_i-1,-1] - 2*scl[idx_i,-1]) / dx**2 - 2*(scl[idx_i,-1]) / dz**2
    #for idx_k in range(1, z_dim-1):
    #    lap[0,  idx_k] = (scl[ 1, idx_k  ] + scl[-1, idx_k  ] - 2*scl[ 0, idx_k]) / dx**2 + \
    #                     (scl[ 0, idx_k+1] + scl[ 0, idx_k-1] - 2*scl[ 0, idx_k]) / dz**2
    #    lap[-1, idx_k] = (scl[ 0, idx_k  ] + scl[-2, idx_k  ] - 2*scl[-1, idx_k]) / dx**2 + \
    #                     (scl[-1, idx_k+1] + scl[-1, idx_k-1] - 2*scl[-1, idx_k]) / dz**2
    #lap[ 0, 0] = (scl[1, 0] + scl[-1, 0] - 2*scl[ 0, 0]) / dx**2 - (2*scl[ 0, 0]) / dz**2
    #lap[ 0,-1] = (scl[1,-1] + scl[-1,-1] - 2*scl[ 0,-1]) / dx**2 - (2*scl[ 0,-1]) / dz**2
    #lap[-1, 0] = (scl[0, 0] + scl[-2, 0] - 2*scl[-1, 0]) / dx**2 - (2*scl[-1, 0]) / dz**2
    #lap[-1,-1] = (scl[0,-1] + scl[-2,-1] - 2*scl[-1,-1]) / dx**2 - (2*scl[-1,-1]) / dz**2
    return lap

def make_non_divergent_conjugate_residual(U_comp, W_comp, dx, dz, eps, max_it, fname, file_counter):
    """ Input:  v and w components of uwlcm velocity field
        Output: v and w velocity components such that d/dxv + d/dzw = 0
        Conjugate residual method assumic cyclic boundary in left-right and closed boundary in top-bottom """

    tmp_U_comp = np.copy(U_comp[0:-1,:])
    tmp_W_comp = np.copy(W_comp[0:-1,:])
    x_dim      = tmp_U_comp.shape[0]
    z_dim      = tmp_U_comp.shape[1]

    dual_u = np.zeros([x_dim+1, z_dim]  )
    dual_w = np.zeros([x_dim  , z_dim+1])
    hlo = 4
    h_u = fill_halo_cyclic(tmp_U_comp, hlo)
    h_u = fill_halo_top(h_u, hlo)
    h_u = fill_halo_bottom_open(h_u, hlo)
    h_w = fill_halo_cyclic(tmp_W_comp, hlo)
    h_w = fill_halo_top(h_w, hlo)
    h_w = fill_halo_bottom_closed(h_w, hlo)
    for idx_i in range(0, x_dim+1):
        for idx_k in range(0, z_dim):
            dual_u[idx_i, idx_k] = 0.5 *(h_u[idx_i+hlo-1, idx_k+hlo] + h_u[idx_i+hlo, idx_k+hlo])
    for idx_k in range(0, z_dim+1):
        for idx_i in range(0, x_dim):
            dual_w[idx_i, idx_k] = 0.5*(h_w[idx_i+hlo, idx_k+hlo-1] + h_w[idx_i+hlo, idx_k+hlo])
 
    # initail guess
    scl_pot    = np.zeros([x_dim, z_dim])
    #err        = -1 * div_h_cyc(tmp_U_comp, tmp_W_comp, dx, dz)
    err        = -1 * div_dual(dual_u, dual_w, dx, dz)
    div_L0_old = 1 
    div_L1_old = 1

    iter_helper   = np.array([])
    L0_div_helper = np.array([])
    L1_div_helper = np.array([])
    L0_err_helper = np.array([])
    L1_err_helper = np.array([])
 
    tmp_iter = 0
    flag = True
    while (flag):

        tmp_iter += 1

        lap_err = lap_h_cyc(err, dx, dz)
        tmp_den = np.sum(lap_err * lap_err)

        if (tmp_den != 0):
            beta = -1 * np.sum(err * lap_err) / tmp_den

        scl_pot  += beta * err
        err      += beta * lap_err

        #u_div, w_div = grad_h_cyc(scl_pot, dx, dz)
        u_div, w_div = grad_dual_cyc(scl_pot, dx, dz)
       
        # calc u_rot = u - u_div
        #divless_v = tmp_U_comp - u_div
        #divless_w = tmp_W_comp - w_div
        divless_v = dual_u - u_div
        divless_w = dual_w - w_div
   
        #div_check     = div_h_cyc(divless_v, divless_w, dx, dz)
        #div_check_ini = div_h_cyc(tmp_U_comp, tmp_W_comp, dx, dz)
        div_check     = div_dual(divless_v, divless_w, dx, dz)
        div_check_ini = div_dual(dual_u, dual_w, dx, dz)
 

        div_L1 = np.longdouble(np.sum(np.abs(div_check)))
        div_L0 = np.longdouble(np.max(np.abs(div_check)))
        print "it= ", tmp_iter,\
              " L0(div) vs L0(err) ", div_L0, " vs ", np.max(np.abs(err)),\
              " L1(div) vs L1(err) ", div_L1, " vs ", np.sum(np.abs(err))

        if (tmp_iter % 100 == 0):
            iter_helper   = np.append(iter_helper,   tmp_iter)
            L0_div_helper = np.append(L0_div_helper, div_L0)
            L1_div_helper = np.append(L1_div_helper, div_L1)
            L0_err_helper = np.append(L0_err_helper, np.max(np.abs(err)))
            L1_err_helper = np.append(L1_err_helper, np.sum(np.abs(err)))

        if (div_L0 < eps or tmp_iter > max_it):
            flag = False
        #if (div_L0 == div_L0_old and div_L1 == div_L1_old):  
        ##    flag = False

        div_L0_old = div_L0
        div_L1_old = div_L1


    if (file_counter % 100 == 0):
        # every 200 files plot whats happening to the velocity fields
        min_v = min(np.min(divless_v), np.min(U_comp))
        max_v = max(np.max(divless_v), np.max(U_comp))
        min_w = min(np.min(divless_w), np.min(W_comp))
        max_w = max(np.max(divless_w), np.max(W_comp))
        zoom = -1
        plt.figure(figsize=(15,10))

        plt.subplot(2,4,1)
        plt.imshow(div_check_ini[:,:zoom], aspect='auto')
        plt.title('initial divergence')
        plt.colorbar()

        plt.subplot(2,4,2)
        plt.imshow(err[:,:zoom], aspect='auto')
        plt.title('solver error')
        plt.colorbar()

        plt.subplot(2,4,3)
        plt.imshow(divless_v[:,:zoom], aspect='auto')
        plt.title('divless_v')
        plt.colorbar()

        plt.subplot(2,4,4)
        plt.imshow(divless_w[:,:zoom], aspect='auto')
        plt.title('divless_w')
        plt.colorbar()

        plt.subplot(2,4,5)
        plt.imshow(div_check[:,:zoom], aspect='auto')
        plt.title('final divergence')
        plt.colorbar()

        plt.subplot(2,4,6)
        plt.imshow(scl_pot[:,:zoom], aspect='auto')
        plt.title('scalar potential')
        plt.colorbar()

        plt.subplot(2,4,7)
        plt.imshow(u_div[:,:zoom], aspect='auto')
        plt.title('div_v')
        plt.colorbar()

        plt.subplot(2,4,8)
        plt.imshow(w_div[:,:zoom], aspect='auto')
        plt.title('div_w')
        plt.colorbar()

        plt.tight_layout()
        plt.savefig("output_conversion_plots/"+str(file_counter)+"_"+str(tmp_iter)+".png")
        plt.close()
  
    U_comp_fin = np.zeros([x_dim+1, z_dim])
    W_comp_fin = np.zeros([x_dim+1, z_dim])
 
    U_comp_fin[0:-1, :] = np.copy(divless_v[1:,:])
    W_comp_fin[0:-1, :] = np.copy(divless_w[:,1:])
    U_comp_fin[-1,   :] = U_comp[0,:]
    W_comp_fin[-1,   :] = W_comp[0,:]

    plt.figure(figsize=(15,10))
    plt.subplot(221)
    plt.title('L0 div')
    plt.plot(iter_helper , L0_div_helper, '.')
    plt.subplot(222)
    plt.title('L1 div')
    plt.plot(iter_helper , L1_div_helper, '.')
    plt.subplot(223)
    plt.title('L0 err')
    plt.plot(iter_helper , L0_err_helper, '.')
    plt.subplot(224)
    plt.title('L1 err')
    plt.plot(iter_helper, L1_err_helper, '.')
    #plt.tight_layout()
    #plt.show()
    plt.savefig("output_conversion_plots/error_convergence_"+str(file_counter)+"_"+str(tmp_iter)+".png")
    plt.close()

    return U_comp_fin, W_comp_fin

# ------------ Here it all begins --------------

# file with initial profiles
fname_stats = "dycoms/Output.DYCOMS_RF01.04b3a/stats/Stats.DYCOMS_RF01.nc"
# file with simulation data from t=0
fname_pkl   = "dycoms/Output.DYCOMS_RF01.04b3a/Visualization/10000010.pkl"
# folder with simulation data
folder_pkl = "dycoms/Output.DYCOMS_RF01.04b3a/Visualization/"

dx           = 35.     # TODO - read in the dx and dz from file
dz           = 5.
eps          = 1e-8
max_it       = 1000
file_counter = 0
halo         = 10

# ------------ initial condition --------------
print "reading initial condition from ", fname_pkl 

# get simulation dimensions
f_pkl      = open(fname_pkl, 'rb')
data       = pkl.load(f_pkl)
pycles_v   = np.array(data['v'])
x_dim      = pycles_v.shape[0]
z_dim      = pycles_v.shape[1]
scl_pot    = np.zeros((x_dim, z_dim+1))

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
## make sure that initial velocity is non-divergent
#uwlcm_v0_nodiv, uwlcm_w0_nodiv = make_non_divergent_conjugate_residual(uwlcm_v0, uwlcm_w0, dx, dz, eps, max_it, fname_pkl, 0)
#
## save t=0 data
#with h5py.File('dycoms/dycoms_init.h5', 'w') as init_hdf:
#    init_hdf.create_dataset("uwlcm_rv0",      data=uwlcm_rv0)
#    init_hdf.create_dataset("uwlcm_v0",       data=uwlcm_v0)
#    init_hdf.create_dataset("uwlcm_w0",       data=uwlcm_w0)
#    init_hdf.create_dataset("uwlcm_v0_nodiv", data=uwlcm_v0_nodiv)
#    init_hdf.create_dataset("uwlcm_w0_nodiv", data=uwlcm_w0_nodiv)
#    init_hdf.create_dataset("uwlcm_thd0",     data=uwlcm_thd0)
#    init_hdf.create_dataset("uwlcm_rhod0",    data=uwlcm_rhod0)
#
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
for fname in os.listdir(folder_pkl):
    if fname.endswith(".pkl"):
        #print file_counter
        #if (file_counter <= 2001):
        #if (file_counter in [7200]):
        if (file_counter == 1400):
            print "reading data from ", fname
            f_pkl   = open(folder_pkl+fname, 'rb')
            data    = pkl.load(f_pkl)

            uwlcm_v = read_horizontal_velocity(data, x_dim, z_dim)
            uwlcm_w = read_vertical_velocity(data,   x_dim, z_dim)

            h_uwlcm_v = fill_halo_cyclic(uwlcm_v[0:-1,:], halo)
            h_uwlcm_v = fill_halo_top(h_uwlcm_v, halo)
            h_uwlcm_v = fill_halo_bottom_open(h_uwlcm_v, halo)
            h_uwlcm_w = fill_halo_cyclic(uwlcm_w[0:-1,:], halo)
            h_uwlcm_w = fill_halo_top(h_uwlcm_w, halo)
            h_uwlcm_w = fill_halo_bottom_closed(h_uwlcm_w, halo)

            sigma_y = 5  #3.0
            sigma_x = 2  #1.0
            sigma = [sigma_x, sigma_y]
            filtered_v_tmp = ndimage.filters.gaussian_filter(h_uwlcm_v, sigma, mode='wrap')
            filtered_w_tmp = ndimage.filters.gaussian_filter(h_uwlcm_w, sigma, mode='wrap')
            filtered_v = np.zeros(uwlcm_v.shape)
            filtered_w = np.zeros(uwlcm_w.shape)
            filtered_v[0:-1,:] = np.copy(filtered_v_tmp[halo:-halo, halo:-halo])
            filtered_w[0:-1,:] = np.copy(filtered_w_tmp[halo:-halo, halo:-halo])
            filtered_v[-1,:] = filtered_v[0,-1]
            filtered_w[-1,:] = filtered_w[0,-1]

            x_dim = uwlcm_v.shape[0]
            z_dim = uwlcm_v.shape[1]
            lin_v = np.zeros(uwlcm_v.shape)
            lin_w = np.zeros(uwlcm_w.shape)
            for idx_i in range(0, x_dim):
                for idx_k in range(0, z_dim):
                    lin_v[idx_i, idx_k] = .25*(h_uwlcm_v[idx_i+1+halo, idx_k+halo] + h_uwlcm_v[idx_i-1+halo, idx_k+halo] + h_uwlcm_v[idx_i+halo, idx_k+1+halo] + h_uwlcm_v[idx_i+halo, idx_k-1+halo])
                    lin_w[idx_i, idx_k] = .25*(h_uwlcm_w[idx_i+1+halo, idx_k+halo] + h_uwlcm_w[idx_i-1+halo, idx_k+halo] + h_uwlcm_w[idx_i+halo, idx_k+1+halo] + h_uwlcm_w[idx_i+halo, idx_k-1+halo])

            uwlcm_v_nodiv, uwlcm_w_nodiv = make_non_divergent_conjugate_residual(\
                uwlcm_v, uwlcm_w,\
                #filtered_v, filtered_w,\
                #lin_v, lin_w,
                dx, dz, eps, max_it, fname, file_counter\
            )
 
            if (file_counter % 200 == 0):
 
                zoom = -1
                tmp_div_ini = div_h_cyc(uwlcm_v, uwlcm_w, dx, dz)
                tmp_div_fin = div_h_cyc(uwlcm_v_nodiv, uwlcm_w_nodiv, dx, dz)

                diff_v = uwlcm_v - uwlcm_v_nodiv
                diff_w = uwlcm_w - uwlcm_w_nodiv

                print "max difference in horizontal velocity: ",  np.max(np.abs(diff_v))
                print "max difference in vertical velocity: ",    np.max(np.abs(diff_w))

                plt.figure(figsize=(15,10))
                plt.subplot(2,3,1)
                plt.imshow(tmp_div_ini[:,:zoom], aspect='auto')
                plt.title('initial diveregence')
                plt.colorbar()
                plt.subplot(2,3,4)
                plt.imshow(tmp_div_fin[:,:zoom], aspect='auto')
                plt.title('final divergence')
                plt.colorbar()
                plt.subplot(2,3,2)
                plt.imshow(diff_v[:,:zoom], aspect='auto')
                plt.title('uwlcm_v - uwlcm_v_nodiv')
                plt.colorbar()
                plt.subplot(2,3,5)
                plt.imshow(diff_w[:,:zoom], aspect='auto')
                plt.title('uwlcm_w - uwlcm_w_nodiv')
                plt.colorbar()
                plt.subplot(2,3,3)
                plt.imshow(uwlcm_v_nodiv[:,:zoom], aspect='auto')
                plt.title('uwlcm_v_nodiv')
                plt.colorbar()
                plt.subplot(2,3,6)
                plt.imshow(uwlcm_w_nodiv[:,:zoom], aspect='auto')
                plt.title('uwlcm_w_nodiv')
                plt.colorbar()
                plt.tight_layout()
                plt.savefig("output_conversion_plots/div_check_"+str(file_counter)+"_"+str(max_it)+".png")
                plt.close()
 
            v.create_dataset(str(file_counter), data=uwlcm_v)
            w.create_dataset(str(file_counter), data=uwlcm_w)
            v_nodiv.create_dataset(str(file_counter), data=uwlcm_v_nodiv)
            w_nodiv.create_dataset(str(file_counter), data=uwlcm_w_nodiv)

        file_counter += 1

output_hdf.close()
