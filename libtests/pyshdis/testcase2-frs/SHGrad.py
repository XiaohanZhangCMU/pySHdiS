import numpy as _np
import scipy as _sp
from scipy.misc import factorial
import matplotlib.pyplot as _plt
import pyshtools as _psh
from SHUtil import l_coeffs, LM_list

###  SHMultiplyC ###
# The function will take two sets of spherical harmonic coefficients, 
# multiply the functions in the space domain, and expand the resulting
# field in spherical harmonics using SHExpandGLQ. The spherical harmonic
# bandwidth of the resulting field is lmax, which is also the bandwidths
# of the input fields.
# The employed spherical harmonic normalization ('4pi','ortho','schmidt')
# and Condon-Shortley phase (1 or -1)
# conventions can be set by the optional arguments norm and csphase;
# if not set, the default is to use geodesy 4-pi normalized harmonics
# that exclude the Condon-Shortley phase of (-1)^m.

def SHMultiplyC(sh1, sh2, lmax=100, norm='4pi', csphase=1):
    SH1 = _psh.SHCoeffs.from_array(sh1, normalization=norm, csphase=csphase)
    SH2 = _psh.SHCoeffs.from_array(sh2, normalization=norm, csphase=csphase)
    grid1 = SH1.expand('GLQ')
    grid2 = SH2.expand('GLQ')
    GridProd = grid1 * grid2
    SHProd = GridProd.expand()
    return SHProd.to_array()

###  VSH1, VSH2 ###
# The function will take a spherical harmonic coefficient class,
# multiply by the normal vector, and return the first/second vector
# spherical harmonic result.

def VSH1(clm):
    grid = clm.expand('GLQ')
    grid_data = grid.data
    latglq, longlq = _psh.expand.GLQGridCoord(clm.lmax)
    LON, LAT = _np.meshgrid(longlq, latglq)
    THETA = _np.deg2rad(90-LAT)
    PHI = _np.deg2rad(LON)

    grid.data = grid_data * _np.sin(THETA) * _np.cos(PHI)
    Y_x = grid.expand()
    grid.data = grid_data * _np.sin(THETA) * _np.sin(PHI)
    Y_y = grid.expand()
    grid.data = grid_data * _np.cos(THETA)
    Y_z = grid.expand()

    return (Y_x, Y_y, Y_z)

def DiffNormCoeffs(lmax, norm=None, csphase=1, shtype='irr'):
    
    if csphase == -1:
        print('csphase == -1 is not implemented!')

    # first we create a list of l, m degrees lower than lmax
    l_list, m_list = LM_list(lmax)
    if shtype == 'irr':
        l_d = l_list + 1
    elif shtype == 'reg':
        l_d = l_list - 1
    else:
        print('invalid shtype (irr, reg)')
        return (0, 0, 0)

    # calculate the normalization coefficient from the Ilm definition
    C_d_p = factorial(l_list + _np.abs(m_list))
    C_d_m = factorial(l_list - _np.abs(m_list))
    C_d_norm_4pi = _np.sqrt((2*l_list + 1)/(2*l_d + 1))
    C_dz_p = factorial(l_d + _np.abs(m_list))
    C_d1_p = factorial(l_d + _np.abs(m_list-1))
    C_d2_p = factorial(l_d + _np.abs(m_list+1))
    C_dz_m = factorial(l_d - _np.abs(m_list))
    C_d1_m = factorial(l_d - _np.abs(m_list-1))
    C_d2_m = factorial(l_d - _np.abs(m_list+1))

    # remove the zeros in the denominator (for regular solid harmonics only)
    C_dz_m_i = _np.zeros(C_dz_m.shape)
    C_d1_m_i = _np.zeros(C_d1_m.shape)
    C_d2_m_i = _np.zeros(C_d2_m.shape)
    C_dz_m_i[C_dz_m != 0.0] = 1.0/C_dz_m[C_dz_m != 0.0]
    C_d1_m_i[C_d1_m != 0.0] = 1.0/C_d1_m[C_d1_m != 0.0]
    C_d2_m_i[C_d2_m != 0.0] = 1.0/C_d2_m[C_d2_m != 0.0]

    C_dz_norm = _np.sqrt(C_dz_p*C_dz_m_i*C_d_m/C_d_p)
    C_d1_norm = _np.sqrt(C_d1_p*C_d1_m_i*C_d_m/C_d_p)
    C_d2_norm = _np.sqrt(C_d2_p*C_d2_m_i*C_d_m/C_d_p)

    # calculate the normalization coefficients of the spherical harmonics
    if csphase==1:
        if norm == '4pi' or norm == 'ortho':
            C_dz_list = C_dz_m/C_d_m * C_dz_norm * C_d_norm_4pi
            C_d1_list = C_d1_m/C_d_m * C_d1_norm * C_d_norm_4pi
            C_d2_list = C_d2_m/C_d_m * C_d2_norm * C_d_norm_4pi
        elif norm == 'schmidt':
            C_dz_list = C_dz_m/C_d_m * C_dz_norm
            C_d1_list = C_d1_m/C_d_m * C_d1_norm
            C_d2_list = C_d2_m/C_d_m * C_d2_norm
        else: # unnormalized
            C_dz_list = C_dz_m/C_d_m
            C_d1_list = C_d1_m/C_d_m
            C_d2_list = C_d2_m/C_d_m
    else: # not implemented
        print('DiffNormCoeffs: csphase = -1 not implemented')
        return (0, 0, 0)

    C_dz = _psh.SHCoeffs.from_zeros(lmax=lmax,kind='complex')
    C_d1 = _psh.SHCoeffs.from_zeros(lmax=lmax,kind='complex')
    C_d2 = _psh.SHCoeffs.from_zeros(lmax=lmax,kind='complex')
    C_dz.set_coeffs(-C_dz_list, l_list + 1, m_list)
    C_d1.set_coeffs(C_d1_list, l_list + 1, m_list - 1)
    C_d2.set_coeffs(-C_d2_list, l_list + 1, m_list + 1)

    return (C_d1.to_array(), C_d2.to_array(), C_dz.to_array())

def ISHgrad(clm, r=1.0):
# Calculate the gradient of irregular solid harmonics
    # normalize the coefficients by radius
    C_d1, C_d2, C_dz = DiffNormCoeffs(clm.lmax, norm=clm.normalization, csphase=clm.csphase)
    lp1 = l_coeffs(clm.lmax) + 1
    clm.coeffs = clm.to_array() / (r**lp1)

    # Dz(I_l^m) = C_dz * I_(l+1)^m
    clm_dz = _psh.SHCoeffs.from_zeros(lmax=clm.lmax,kind='complex',normalization=clm.normalization,csphase=clm.csphase)
    clm_dz.coeffs[:,1:,:] = clm.coeffs[:,:-1,:]
    clm_dz.coeffs = clm_dz.to_array() * C_dz

    # D1(I_l^m) = C_d1 * I_(l+1)^(m-1)
    clm_d1 = _psh.SHCoeffs.from_zeros(lmax=clm.lmax,kind='complex',normalization=clm.normalization,csphase=clm.csphase)
    clm_d1.coeffs[0,1:,:-1] = clm.coeffs[0,:-1,1:]
    clm_d1.coeffs[1,1:,1:] = clm.coeffs[1,:-1,:-1]
    clm_d1.coeffs[1,1:,1] = clm.coeffs[0,:-1,0]
    clm_d1.coeffs = clm_d1.to_array() * C_d1

    # D2(I_l^m) = C_d2 * I_(l+1)^(m+1)
    clm_d2 = _psh.SHCoeffs.from_zeros(lmax=clm.lmax,kind='complex',normalization=clm.normalization,csphase=clm.csphase)
    clm_d2.coeffs[0,1:,1:] = clm.coeffs[0,:-1,:-1]
    clm_d2.coeffs[1,1:,:-1] = clm.coeffs[1,:-1,1:]
    clm_d2.coeffs[0,1:,0] = clm.coeffs[1,:-1,1]
    clm_d2.coeffs = clm_d2.to_array() * C_d2

    # Dx = (D1 + D2)/2; Dy = (D2 - D1)/2j;
    clm_dx = _psh.SHCoeffs.from_zeros(lmax=clm.lmax,kind='complex',normalization=clm.normalization,csphase=clm.csphase)
    clm_dy = _psh.SHCoeffs.from_zeros(lmax=clm.lmax,kind='complex',normalization=clm.normalization,csphase=clm.csphase)
    clm_dx.coeffs = (clm_d1.to_array() + clm_d2.to_array())/2.0
    clm_dy.coeffs = (clm_d2.to_array() - clm_d1.to_array())/2.0j

    return (clm_dx, clm_dy, clm_dz)

def VSH2(clm):
    clm_dx, clm_dy, clm_dz = ISHgrad(clm)

    # the first spherical harmonic part
    l_list, m_list = LM_list(clm.lmax+1)
    lp1 = _np.swapaxes(_np.broadcast_to(_np.arange(clm.lmax+1)+1, (2, clm.lmax+1, clm.lmax+1)), 1,2)
    clm_lp1 = clm.copy()
    clm_lp1.coeffs = clm.to_array()*lp1
    lp1_Y_x, lp1_Y_y, lp1_Y_z = VSH1(clm_lp1)
    #print(clm_dz.to_array())
    #print(Y_z.to_array())
    #print(C_lp1.to_array())
    clm_dx.coeffs = clm_dx.to_array() + lp1_Y_x.to_array()
    clm_dy.coeffs = clm_dy.to_array() + lp1_Y_y.to_array()
    clm_dz.coeffs = clm_dz.to_array() + lp1_Y_z.to_array()

    return (clm_dx, clm_dy, clm_dz)
