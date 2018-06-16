import numpy as _np
import scipy as _sp
import matplotlib.pyplot as _plt
import pyshtools as _psh

def CartCoord_to_SphCoord(X, Y, Z):
    # translate Cartesian coordinates into spherical coordinates
    # make sure X, Y, Z have same dimensions
    R = _np.sqrt(X**2 + Y**2 + Z**2)
    THETA = _np.arccos(Z/R)
    PHI = _np.arctan2(Y, X)
    return (R, THETA, PHI)

def SphCoord_to_CartCoord(R, THETA, PHI):
    # translate spherical coordinates into Cartesian coordinates
    # make sure R, THETA, PHI have same dimensions
    Z = R * _np.cos(THETA)
    X = R * _np.sin(THETA) * _np.cos(PHI)
    Y = R * _np.sin(THETA) * _np.sin(PHI)
    return (X, Y, Z)

def LM_list(lmax):
    l_list = _np.empty((lmax)**2, dtype=_np.int)
    m_list = _np.empty((lmax)**2, dtype=_np.int)
    for l in range(lmax):
        l_list[l**2:(l+1)**2] = l
        m_list[l**2:(l+1)**2] = _np.arange(-l, l+1)
    return (l_list, m_list)

def l_coeffs(lmax):
    # This function to create the l degrees corresponding to the coefficient matrix shape
    l_list = _np.arange(lmax + 1)
    l_coeffs_T = _np.broadcast_to(l_list, (2, lmax+1, lmax+1))
    return _np.swapaxes(l_coeffs_T, 1, 2)
    # l_coeffs = _np.swapaxes(l_coeffs_T, 1, 2)
    # return l_coeffs

###  eval_GridC ###
# The function will take a set of points to evaluate the Spherical harmonic functions. 
# coeff: SHCoeffs class
# latin, lonin: latitude and longitude coordinates of the points to evaluate, in DEGREE
# rin: 
#     scalar value - evaluate on sphere with r; 
#     array like - evaluate on 3D points (must have same dimensions as latin and lonin)
# lmax_calc: the max l degree to calculate. If not specified, evaluate entire coeff
# norm: normalization convention. If not specified, use the convention same as coeff
#     1 - Geodesy 4-pi; 2 - Schimidt Semi-orthogonalized; 3 - Unnormalized; 4 - Orthonormal.

def eval_GridC(coeff, latin, lonin, rin=1.0, lmax_calc=None, norm=None, shtype=None):
    if shtype == None:
        C_n = 0.0
    elif shtype == 'irr':
        C_n = -l_coeffs(coeff.lmax)-1
    elif shtype == 'reg':
        C_n = l_coeffs(coeff.lmax)
    elif type(shtype) is int:
        C_n = shtype
    else:
        print('invalid shtype (irr, reg, float)')
    if norm == None: # normalization not given, use the coefficient norm convention
        if coeff.normalization == '4pi':
            norm = 1
        elif coeff.normalization == 'schmidt':
            norm = 2
        elif coeff.normalization == 'ortho':
            norm = 4
    if lmax_calc == None:
        lmax_calc = coeff.lmax
    if (type(rin) is _np.ndarray) or (type(rin) is list):
        values = _np.empty_like(latin, dtype=_np.complex)
        for v, latitude, longitude, radius in _np.nditer([values, latin, lonin, rin], op_flags=['readwrite']):
            C_r = radius**C_n
            v[...] = _psh.expand.MakeGridPointC(C_r*coeff.to_array(),
                                                lat=latitude, lon=longitude, 
                                                lmax=lmax_calc, norm=norm, 
                                                csphase=coeff.csphase)
    else:
        values = _np.empty_like(latin, dtype=_np.complex)
        C_r = rin**C_n
        for v, latitude, longitude in _np.nditer([values, latin, lonin], op_flags=['readwrite']):
            v[...] = _psh.expand.MakeGridPointC(C_r*coeff.to_array(),
                                                lat=latitude, lon=longitude, 
                                                lmax=lmax_calc, norm=norm, 
                                                csphase=coeff.csphase)
    return values

def VSH_vector_add(v1, s1, v2, s2):
    # v = s1*v1 + s2*v2
    return [_psh.SHCoeffs.from_array(v1[i].to_array()*s1 + v2[i].to_array()*s2) for i in range(len(v1))]

def SHCilmToVector(cilm, lmax):
    vec = _np.zeros((lmax+1)**2)*0j
    #print(vec.shape)
    for l in range(lmax+1):
        #print(cilm[1,l,l:0:-1])
        vec[(l**2):(l**2+l)] = cilm[1, l, l:0:-1]
        #print(cilm[0,l,:l+1])
        #print(l, l**2+l+1, (l+1)**2)
        vec[(l**2+l):((l+1)**2)] = cilm[0, l, 0:l+1]
    return vec

########### savemode, loadmode: save and load vectorized spherical harmonic mode (l, m, value)

def sparse_mode(M1, etol=1e-8):
    _np.set_printoptions(suppress=True)
    #print(M1.to_array())
    l_mat = _psh.SHCoeffs.from_zeros(M1.lmax)
    m_mat = _psh.SHCoeffs.from_zeros(M1.lmax)
    l_list, m_list = LM_list(M1.lmax+1)
    l_mat.set_coeffs(l_list, l_list, m_list)
    m_mat.set_coeffs(m_list, l_list, m_list)
    mode_mat = M1.to_array()
    mode_idx = _np.logical_or(_np.abs(_np.real(mode_mat)) > etol, _np.abs(_np.imag(mode_mat)) > etol)
    values = mode_mat[mode_idx]
    ms = m_mat.to_array()[mode_idx]
    ls = l_mat.to_array()[mode_idx]
    sparseM = _np.vstack((ls, ms, values))
    #print(sparseM)
    return sparseM

def savemode(M, d=1, etol=1e-8):
    if d==0:
        mode = sparse_mode(M, etol)
    elif d==1:
        raw_mode = [None for _ in range(3)]
        for i in range(3):
            mode_i = sparse_mode(M[i], etol)
            idx_i= _np.ones(mode_i.shape[1])*i
            raw_mode[i] = _np.vstack((idx_i, mode_i))
        mode = _np.hstack(raw_mode)
    elif d==2:
        raw_mode = [None for _ in range(9)]
        for i in range(3):
            for j in range(3):
                mode_j = sparse_mode(M[i][j], etol)
                idx_i  = _np.ones(mode_j.shape[1])*i
                idx_j  = _np.ones(mode_j.shape[1])*j
                raw_mode[i*3+j] = _np.vstack((idx_i, idx_j, mode_j))
        mode = _np.hstack(raw_mode)
    else:
        print('savemode: invalid dimension d=', d)
        return 0

    return mode

def loadmode(mode, lmax):
    d = mode.shape[0] - 3
    if d==0:
        M = _psh.SHCoeffs.from_zeros(lmax, kind='complex')
        l_list = _np.array(_np.real(mode[d+0, :]), dtype=_np.int)
        m_list = _np.array(_np.real(mode[d+1, :]), dtype=_np.int)
        values = mode[d+2, :]
        M.set_coeffs(values, l_list, m_list)
    elif d==1:
        M = [None for _ in range(3)]
        for i in range(3):
            M[i] = _psh.SHCoeffs.from_zeros(lmax, kind='complex')
            colidx = (_np.array(_np.real(mode[0, :]), dtype=_np.int) == i)
            l_list = _np.array(_np.real(mode[d+0, colidx]), dtype=_np.int)
            m_list = _np.array(_np.real(mode[d+1, colidx]), dtype=_np.int)
            values = mode[d+2, colidx]
            M[i].set_coeffs(values, l_list, m_list)
    elif d==2:
        M = [[None for _ in range(3)] for _ in range(3)]
        for i in range(3):
            for j in range(3):
                M[i][j] = _psh.SHCoeffs.from_zeros(lmax, kind='complex')
                colidx = _np.logical_and((_np.array(_np.real(mode[0, :]), dtype=_np.int) == i),\
                                         (_np.array(_np.real(mode[1, :]), dtype=_np.int) == j))
                l_list = _np.array(_np.real(mode[d+0, colidx]), dtype=_np.int)
                m_list = _np.array(_np.real(mode[d+1, colidx]), dtype=_np.int)
                values = mode[d+2, colidx]
                M[i][j].set_coeffs(values, l_list, m_list)
    else:
        print('loadmode: invalid dimension d=', d)
        return 0
    
    return M, d