import numpy as np
from scipy.io import loadmat, savemat
import scipy.sparse as spm
import pyshtools
from SHUtil import SHCilmToVector
from SHBV import spec_J2lmk, spec_lmk2J, subCmat, print_SH_mode, create_Cmat, stress_solution
from ShElastic import T_mode, S_mode

#### load the image stress problem setting
print('load the image stress problem setting')
data = loadmat('Susr.mat')
sigma_inf = np.swapaxes(data['S'], 0, 1)
a = data['R']
mu = data['MU']
nu = data['NU']
b = 1
b1 = data['b1']
x1 = data['x1']
x2 = data['x2']

#### generate meshing on the void surface
print('generate meshing on the void surface')
Ngrid = sigma_inf.shape[1]
theta = np.arange(0, np.pi, np.pi/Ngrid)
phi = np.arange(0, 2*np.pi, np.pi/Ngrid)
THETA, PHI = np.meshgrid(theta, phi)

Z = a*np.cos(THETA)
Y = a*np.sin(THETA)*np.sin(PHI)
X = a*np.sin(THETA)*np.cos(PHI)
N = -np.stack([X/a, Y/a, Z/a], axis=-1)

#### decompose traction boundary condition
print('decompose traction boundary condition')
T_inf = N*0
for i,x in np.ndenumerate(THETA):
    T_inf[i] = np.dot(sigma_inf[i], N[i])
T_usr_mesh = T_inf

T_usr_vec = np.array([])
lmax_sub = np.int(Ngrid/2) - 1
mode_sub = lmax_sub - 3

mk = '*x.'
for k in range(3):
    T_usr_cilm = pyshtools.expand.SHExpandDHC(T_usr_mesh[:,:,k].T, sampling=2, lmax_calc=lmax_sub)
    T_usr_mode = pyshtools.SHCoeffs.from_array(T_usr_cilm)
    T_usr_grid = T_usr_mode.expand()
    T_usr_vec = np.hstack((T_usr_vec, SHCilmToVector(T_usr_cilm, lmax = lmax_sub)))

#### load traction mode matrix
print('load traction mode matrix')

lmax_full = 25
lmax_mode = 22
#Cmat_file = 'Cmat_lmax%d_mode%d_mu%f_nu%f.mat' % (lmax_full, lmax_mode, mu, nu)
Csub_file = 'Csub_lmax%d_mode%d_mu%f_nu%f.mat' % (lmax_sub,  mode_sub,  mu, nu)
m_max = 3
#Cmat = create_Cmat(lmax_full, lmax_mode, mu, nu, m_max=m_max, Cmat_file=Cmat_file, recalc=False, etol=1e-15)
#Csub = subCmat(Cmat, lmax_full, lmax_mode, lmax_sub, mode_sub, m_max=m_max, Csub_file=Csub_file, verbose=False)
Csub = loadmat(Csub_file)['Cmat']

#### determine the mode coefficients
print('determine the mode coefficients')
A = spm.linalg.lsqr(Csub, T_usr_vec.transpose())
A_sol = A[0]
index_sol = print_SH_mode(A_sol, m_dir=3, etol=1e-10, verbose=False)

#### evaluate stress solution and nodal force
print('evaluate stress solution and nodal force')
nseg, m = b1.shape
f1 = np.zeros(x1.shape)
f2 = np.zeros(x2.shape)
xi = x2 - x1
sigma_1 = stress_solution(index_sol, x1[:, 0], x1[:, 1], x1[:, 2], MU=mu, NU=nu, lmax=lmax_sub, recalc=False, verbose=False).real
sigma_2 = stress_solution(index_sol, x2[:, 0], x2[:, 1], x2[:, 2], MU=mu, NU=nu, lmax=lmax_sub, recalc=False, verbose=False).real

for i in range(nseg):
    f1[i, :] = np.cross(np.dot(sigma_1[i, :, :], b1[i, :]), xi[i, :])
    f2[i, :] = np.cross(np.dot(sigma_2[i, :, :], b1[i, :]), xi[i, :])

savemat('fimg.mat', {'fi0': f1, 'fi1': f2})
