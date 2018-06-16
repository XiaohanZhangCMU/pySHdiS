import numpy as _np
import scipy as _sp
import scipy.sparse as _spm
import matplotlib.pyplot as _plt
import pyshtools as _psh
import ShElastic as _she
from SHUtil import SHCilmToVector
import time as _time
import sys, os

###  spec_lmk2J, spec_J2lmk ###
# Given lmax, transform between (l,m,k) degrees and column number K

def spec_lmk2J(l, m, k, lmax=100):
    nk = (lmax+1)**2
    return k*nk + l**2 + (m + l)

def spec_J2lmk(J, lmax=100):
    nk = (lmax+1)**2
    k = _np.int(J/nk)
    l = _np.int(_np.sqrt(J - k*nk))
    m = J - k*nk - l**2 - l
    return (l, m, k)

def save_mode(filename, mode):
    print('saving mode:', filename)
    mode_arr = _np.empty((3,3)+mode[0][0].to_array().shape, dtype=_np.complex)
    for i in range(3):
        for j in range(3):
            mode_arr[i, j] = mode[i][j].to_array()
    _np.save(filename, mode_arr)

def load_mode(filename, axis=2):
    data = _np.load(filename)
    mode_arr = data['mode_arr']
    #mode = _np.empty(mode_arr.shape[:axis])
    #for i, _ in _np.ndenumerate(mode):
    #    mode[i] = _psh.SHCoeffs.from_array(mode_arr[i,:])
    data.close()
    return mode_arr

# create full matrix C
def create_Cmat(lmax, mode_lmax, MU=100, NU=0.25, m_max=3, Cmat_file=None, recalc=True, etol=1e-8):
    if not (Cmat_file is None):
        if recalc or (not os.path.exists(Cmat_file)):

            assert 0, 'python2 is not tested. Cmat from python3 is required.'

            print('Integrating traction modes to a matrix')
            size_row = 3*(lmax+1)**2
            size_col = m_max*(mode_lmax+1)**2
            print(size_row, size_col)
            Cmat = _spm.lil_matrix( (size_row, size_col), dtype=_np.complex_)
            c_omega = 1.0

            for k in range(m_max):
                tic_k = _time.time()
                for l in range(mode_lmax+1):
                    tic_l = _time.time()
                    for m in range(-l, l+1):
                        tic = _time.time()
                        Tlm_vec = _np.empty(3*(lmax+1)**2)*0j
                        T_K, S_K = _she.T_mode(l, m, k, MU, NU, c_omega=c_omega, shtype='irr',\
                                               lmax=lmax, stress=True, recalc=recalc, etol=etol)
                        for i in range(3):
                            Tlm_vec_i = SHCilmToVector(T_K[i].to_array(), lmax=lmax)
                            Tlm_vec[i*(lmax+1)**2:(i+1)*(lmax+1)**2] = Tlm_vec_i
                        Cmat[:, spec_lmk2J(l, m, k, lmax=mode_lmax)] = _np.matrix(Tlm_vec).T
                        toc = _time.time()
                        print(k, l, m, toc-tic)
                    toc_l = _time.time()
                    print (k, l, toc_l-tic_l)
                toc_k = _time.time()
                print ('k =', k, toc_k-tic_k)
            _sp.io.savemat(Cmat_file, {'Cmat': Cmat})
        else:
            Cmat = _sp.io.loadmat(Cmat_file)['Cmat']
    else:
        print('create_Cmat: Must provide matrix file name!')
        Cmat = -1

    return Cmat

# sub matrix of C
def subCmat(Cmat, lmax, mode_lmax, lmax_sub, mode_sub, m_max=3, Csub_file=None, verbose=False):
    if (not (Csub_file is None)) and (os.path.exists(Csub_file)):
        Csub = _sp.io.loadmat(Csub_file)['Cmat']
    else:
        size_row_sub = 3*(lmax_sub+1)**2
        size_col_sub = m_max*(mode_sub+1)**2
        Csub = _spm.lil_matrix( (size_row_sub, size_col_sub), dtype=_np.complex_ )

        for m_dir in range(m_max):
            for l in range(mode_sub+1):
                for m in range(-l, l+1):
                    if verbose:
                        print(m_dir, l, m)
                    for k in range(3):
                        nk = (lmax_sub+1)**2
                        NK = (lmax+1)**2
                        J_full = spec_lmk2J(l, m, m_dir, lmax=mode_lmax)
                        J_sub = spec_lmk2J(l, m, m_dir, lmax=mode_sub)
                        Tlm = Cmat[k*NK:k*NK+nk, J_full]
                        Csub[k*nk:(k+1)*nk, J_sub] = Tlm
        if not (Csub_file is None):
            _sp.io.savemat(Csub_file, {'Cmat': Csub})
    return Csub

def visualize_Cmat(Csub, precision=0, m_max=4):
    _plt.spy(Csub, precision=precision, markersize = 3)
    mode_sub = _np.int(_np.sqrt(Csub.shape[1]/m_max))-1
    lmax_sub = _np.int(_np.sqrt(Csub.shape[0]/3))-1
    print(mode_sub, lmax_sub)
    lsy = _np.arange(0, lmax_sub+1)
    lsx = _np.arange(0, mode_sub+1)

    x_range = [0, Csub.shape[1]]
    y_range = [0, Csub.shape[0]]
    # traction direction dividing line:
    for i in range(1, 3+1):
        y_divide = y_range[1]*i/3
        _plt.plot(x_range, _np.ones(2)*y_divide-0.5)
        y_text = y_divide - y_range[1]/m_max/2
        _plt.text(-mode_sub-5, y_text, '$T_'+str(i)+'$', fontsize=24)
    # mode dividing line:
    for i in range(1, m_max+1):
        x_divide = x_range[1]*i/m_max
        _plt.plot(_np.ones(2)*x_divide-0.5, y_range)
        x_text = x_divide - x_range[1]/3/2 - 1
        _plt.text(x_text, -lmax_sub, '$\\psi_'+str(i)+'$', fontsize=24)

    # ticks for different traction directions
    ticks_y = _np.array([])
    for i in range(3):
        ticks_Ti = i*(lmax_sub+1)**2+lsy**2
        ticks_y = _np.hstack((ticks_y, ticks_Ti))
    _plt.yticks(ticks_y, _np.tile(lsy, 3))
    # ticks for different modes
    ticks_x = _np.array([])
    for i in range(m_max):
        ticks_Psi_i = i*(mode_sub+1)**2+lsx**2
        ticks_x = _np.hstack((ticks_x, ticks_Psi_i))
    _plt.xticks(ticks_x, _np.tile(lsx, m_max))
    if m_max == 4:
        _plt.title('Solution A+B', fontsize=24)
    else:
        _plt.title('Solution B', fontsize=24)

def print_SH_mode(vec, m_dir=4, etol=1e-8, verbose=False):
    # vec is a *complex* spherical harmonic vector
    # with m different directions
    idx_type = [('index', '<i4', 3), ('coeff', _np.complex_)]
    idx_mode = _np.array([], dtype=idx_type)
    lmax = _np.int(_np.sqrt(vec.size / m_dir)) - 1
    idx = [spec_J2lmk(i, lmax) for i in range(vec.size)]
    for i in range(vec.size):
        if (_np.abs(_np.real(vec[i])) > etol or \
            _np.abs(_np.imag(vec[i])) > etol):
            if verbose:
                print('index:',i, idx[i], 'coeff:', vec[i])
            new_idx = _np.array((idx[i], vec[i]), dtype=idx_type)
            idx_mode = _np.append(idx_mode, new_idx)
    return idx_mode

def stress_solution(index_sol, X, Y, Z, MU=100, NU=0.25, lmax=100, recalc=False, verbose=False):
    R = _np.sqrt(X**2+Y**2+Z**2)
    assert R.all()>0, (R, X, Y, Z)
    THETA = _np.arccos(Z/R)
    PHI = _np.arctan2(Y, X)
    sigma_tot = _np.zeros(X.shape+(3,3))*(0j)
    Slm_dir = 'Slm_SH_mu'+str(MU)+'_nu'+str(NU)
    if verbose:
        print('the stress solution includes: ')
        print('l m k coeff')
    for idx_sol in index_sol:
        l, m, k = idx_sol['index']
        c_omega = idx_sol['coeff']
        if verbose:
            print(l, m, k, c_omega)
        sigma = _np.zeros(X.shape+(3,3))*(0j)
        Slm = _she.S_mode(l, m, k, MU, NU, c_omega=c_omega, shtype='irr', lmax=lmax, recalc=recalc)
        for idx in _np.ndindex(X.shape):
            lat_d = 90-THETA[idx]/_np.pi*180
            lon_d = PHI[idx]/_np.pi*180
            for u in range(3):
                for v in range(3):
                    Slm_uv = Slm[u][v].to_array()
                    sigma_uv = _psh.expand.MakeGridPointC(Slm_uv, lat_d, lon_d)
                    sigma[idx+(u,v)] = sigma_uv
            if k == 3:
                sigma[idx] /= R[idx]**(l+3)
            else:
                sigma[idx] /= R[idx]**(l+2)
        sigma_tot = sigma_tot + sigma # sigma(r,theta,phi)
    return sigma_tot
