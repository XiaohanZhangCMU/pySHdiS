# coding: utf-8

import time
import numpy as np
from scipy.misc import factorial
import sys
if sys.version_info[0] < 3:
  sys.path.append('../../../SHTOOLS-4.0')
else:
  sys.path.append('../../../SHTOOLS-4.1')
import pyshtools


# ## Implementation of ShElastic Solution
# 
# First we import the libraries we developed:

#import numpy as np
import pyshtools
from SHUtil import SHCilmToVector
from SHBV import spec_J2lmk, spec_lmk2J, subCmat, print_SH_mode, create_Cmat, visualize_Cmat, stress_solution
from ShElastic import T_mode, S_mode


# Then we generate the meshgrid on the void surface for boundary conditions. 

#import matplotlib.pyplot as plt
#from matplotlib import cm, colors
#from mpl_toolkits.mplot3d import Axes3D
import scipy.sparse as spm
from scipy.io import loadmat, savemat

def void_mesh(Ngrid, mu, nu, a, b, t):
    theta = np.arange(0, np.pi, np.pi/Ngrid)
    phi = np.arange(0, 2*np.pi, np.pi/Ngrid)
    THETA, PHI = np.meshgrid(theta, phi)
    Zm = a*np.cos(THETA)
    Ym = a*np.sin(THETA)*np.sin(PHI)
    Xm = a*np.sin(THETA)*np.cos(PHI)
    return (Xm,Ym,Zm)

def stress_soln_void_dislocation(Ngrid, mu, nu, a, b, t, XX,YY,ZZ, X, Y, Z, sigma_inf):
#### Willis, Hayns & Bullough (1972) ####
    mu0 = mu; mu = 1;
    r0 = a; b0 = b; x01 = X; x02 = Y; x03 = Z; sigma_0 = sigma_inf;
    a = a/r0; b = b0/r0; X = x01/r0; Y = x02/r0; Z = x03/r0; sigma_inf = sigma_0/mu0
#    r0 = a; b0 = b; a = a/r0; b = b0/r0; X /=r0; Y /=r0; Z /=r0; sigma_inf /=mu0


    #x3 = np.linspace(0, 3, 30)
    lmax_full = 15
    lmax_mode = 12
    Cmat_file = 'Cmat_lmax%d_mode%d_mu%f_nu%f.mat' % (lmax_full, lmax_mode, mu, nu)

    N = -np.stack([XX/r0, YY/r0, ZZ/r0], axis=-1)
        
    T_inf = N*0
    for i,x in np.ndenumerate(XX):
        T_inf[i] = np.dot(sigma_inf[i], N[i])
    T_usr_mesh = T_inf
    # Then we expand the traction boundary conditions to spherical harmonic modes:
    
    T_usr_vec = np.array([])
    lmax_sub = np.int(Ngrid/2) - 1
    mode_sub = lmax_sub - 3
    np.set_printoptions(suppress=True)
    
    mk = '*x.'
    for k in range(3):
        T_usr_cilm = pyshtools.expand.SHExpandDHC(T_usr_mesh[:,:,k].T, sampling=2, lmax_calc=lmax_sub)
        T_usr_mode = pyshtools.SHCoeffs.from_array(T_usr_cilm)
        T_usr_grid = T_usr_mode.expand()
        T_usr_vec = np.hstack((T_usr_vec, SHCilmToVector(T_usr_cilm, lmax = lmax_sub)))
    
    
    m_max = 3
    Cmat = create_Cmat(lmax_full, lmax_mode, mu, nu, m_max=m_max, Cmat_file=Cmat_file, recalc=False, etol=1e-10)
    #Cmat = loadmat(Cmat_file)['Cmat']
    #savemat(Cmat_file, {'Cmat': Cmat})
    
    Csub = subCmat(Cmat, lmax_full, lmax_mode, lmax_sub, mode_sub, m_max=m_max)
    #plt.figure(figsize=(24,24))
    #visualize_Cmat(Csub, precision=1e-8, m_max=m_max)
    #plt.show()
    #savemat(Cmat_file, {'Cmat': Cmat})
    
    Csub_copy = spm.lil_matrix(Csub.shape, dtype=np.complex)
    etol = 1e-8
    real_idx = np.abs(np.real(Csub)) > etol
    imag_idx = np.abs(np.imag(Csub)) > etol
    Csub_copy[real_idx] = Csub[real_idx]
    Csub_copy[imag_idx] = Csub[imag_idx]
    Csub = Csub_copy
    
    tic = time.time()
    A = spm.linalg.lsqr(Csub, T_usr_vec.transpose())
    toc = time.time()
    A_sol = A[0]
    index_sol = print_SH_mode(A_sol, m_dir=3)
    
    #X, Y, Z = np.meshgrid([t, ], [0,], x3)
    sigma_tot = stress_solution(index_sol, X, Y, Z, MU=mu, NU=nu, lmax=lmax_sub, recalc=False).real * mu0

    return sigma_tot
