import numpy as _np
import scipy as _sp
import matplotlib.pyplot as _plt
import pyshtools as _psh
import os
from SHUtil import VSH_vector_add, savemode, loadmode
from SHGrad import VSH1, VSH2

def U_mode(l, m, k, mu=0.5, nu=None, c_omega=1.0, r=1.0, shtype='irr', lmax=100, recalc=False, etol=1e-8):
    dir_Umode  = 'Umode'
    file_Umode = os.path.join(dir_Umode, 'U_l%d_m%d_k%d.npz' % (l, m, k))
    if os.path.exists(file_Umode) and (not recalc):
        U_nu, d = loadmode(_np.load(file_Umode)['load_Unu'], lmax)
        U_0,  d = loadmode(_np.load(file_Umode)['load_U0'] , lmax)
        if d != 1:
            print('U_mode: invalid dimension d=', d)
            return 0
    else:
        assert 0, 'python2 is not tested. U_mode from python3 is required.'
        if shtype == 'irr':
            C_r = r**(-l-1)
            C_vsh1 = -l
        elif shtype == 'reg':
            C_r = r**l
            C_vsh1 = l+1
        else:
            print('U_mode: invalid shtype (irr, reg)')
            return 0
        U_nu = [_psh.SHCoeffs.from_zeros(lmax=lmax, kind='complex') for _ in range(3)]
        C_U = C_r*c_omega
        U_nu[k].set_coeffs(C_U, l, m)
        Ylm = _psh.SHCoeffs.from_zeros(lmax, kind='complex')
        Ylm.set_coeffs(1.0, l, m)
        rdotPsi = VSH1(Ylm)[k]
        U_0 = VSH_vector_add(VSH1(rdotPsi), C_vsh1*C_U, VSH2(rdotPsi), C_U)

        save_Unu = savemode(U_nu, d=1, etol=etol)
        save_U0  = savemode(U_0,  d=1, etol=etol)
        if not os.path.exists(dir_Umode):
            os.makedirs(dir_Umode)
        _np.savez(file_Umode, load_Unu=save_Unu, load_U0=save_U0)

    if (nu == None):
        return (U_nu, U_0)
    else:
        C_nu = -4*(1-nu)
        return VSH_vector_add(U_nu, C_nu/2.0/mu, U_0, 1/2.0/mu)

def gradU(U, c1, c2):
    gradu = [None for _ in range(3)]
    for i in range(3):
        Ui_VSH1 = VSH1(U[i])
        Ui_VSH2 = VSH2(U[i])
        gradu[i] = VSH_vector_add(Ui_VSH1, c1, Ui_VSH2, c2)
        # gradu[i] = [_psh.SHCoeffs.from_array(c1*U_VSH1[j].to_array()+c2*U_VSH2[j].to_array()) for j in range(3)]

    return gradu
    

def S_mode(l, m, k, mu=0.5, nu=None, c_omega=1.0, r=1.0, shtype='irr', lmax=100, recalc=False, etol=1e-8):
    dir_Smode  = 'Smode'
    file_Smode = os.path.join(dir_Smode, 'S_l%d_m%d_k%d.npz' % (l, m, k))
    if os.path.exists(file_Smode)  and (not recalc):
        S_nu1, d = loadmode(_np.load(file_Smode)['load_Snu1'], lmax)
        S_nu2, d = loadmode(_np.load(file_Smode)['load_Snu2'], lmax)
        S_nu3, d = loadmode(_np.load(file_Smode)['load_Snu3'], lmax)
        S_0,   d = loadmode(_np.load(file_Smode)['load_S0'],   lmax)
        #print('S_nu1\n', _np.load(file_Smode)['load_Snu1'])
        #print('S_nu2\n', _np.load(file_Smode)['load_Snu2'])
        #print('S_nu3\n', _np.load(file_Smode)['load_Snu3'])
        #print('S_0\n', _np.load(file_Smode)['load_S0'])
        if d != 2:
            print('S_mode: invalid dimension d=', d)
            return 0
    else:
        assert 0, 'python2 is not tested. S_mode from python3 is required.'
        if shtype == 'irr':
            c1 = -(l+1)/r**(l+2)
            c2 = 1/r**(l+2)
        elif shtype == 'reg':
            c1 = l*r**(l-1)
            c2 = r**(l-1)
        else: # if (shtype != 'irr') and (shtype != 'reg'):
            print('S_mode: invalid shtype (irr, reg)')
            return 0

        S_nu1= [None for _ in range(3)]
        S_nu2= [None for _ in range(3)]
        S_nu3= [None for _ in range(3)]
        S_0  = [None for _ in range(3)]
        U_nu, U_0 = U_mode(l, m, k, mu=mu, nu=None, c_omega=1.0, r=r, shtype=shtype, lmax=lmax, recalc=recalc, etol=etol)
        gradu_nu = gradU(U_nu, c1, c2)
        gradu_0  = gradU(U_0, c1, c2)
        ukk_nu_arr = gradu_nu[0][0].to_array()+gradu_nu[1][1].to_array()+gradu_nu[2][2].to_array()
        ukk_0_arr  = gradu_0[0][0].to_array() +gradu_0[1][1].to_array() +gradu_0[2][2].to_array()
        for i in range(3):
            S_nu1[i] = [_psh.SHCoeffs.from_array(int(i==j)*ukk_nu_arr) for j in range(3)]
            S_nu2[i] = [_psh.SHCoeffs.from_array(0.5*(gradu_nu[i][j].to_array()+gradu_nu[j][i].to_array())) for j in range(3)]
            S_nu3[i] = [_psh.SHCoeffs.from_array(int(i==j)*ukk_0_arr) for j in range(3)]
            S_0[i]   = [_psh.SHCoeffs.from_array(0.5*(gradu_0[i][j].to_array() +gradu_0[j][i].to_array())) for j in range(3)]
        '''
            for j in range(3):
                print('i=%d, j=%d, Snu1, Snu2, Snu3, S0' % (i, j))
                fig, ax = S_nu1[i][j].plot_spectrum2d(vrange=(1e-8, 1))
                fig, ax = S_nu2[i][j].plot_spectrum2d(vrange=(1e-8, 1))
                fig, ax = S_nu3[i][j].plot_spectrum2d(vrange=(1e-8, 1))
                fig, ax = S_0[i][j].plot_spectrum2d(vrange=(1e-8, 1))
        '''

        save_Snu1 = savemode(S_nu1, d=2, etol=etol)
        save_Snu2 = savemode(S_nu2, d=2, etol=etol)
        save_Snu3 = savemode(S_nu3, d=2, etol=etol)
        save_S0   = savemode(S_0,   d=2, etol=etol)
        #print('S_nu1\n', save_Snu1)
        #print('S_nu2\n', save_Snu2)
        #print('S_nu3\n', save_Snu3)
        #print('S_0\n', save_S0)
        if not os.path.exists(dir_Smode):
            os.makedirs(dir_Smode)
        _np.savez(file_Smode, load_Snu1=save_Snu1, load_Snu2=save_Snu2, load_Snu3=save_Snu3, load_S0=save_S0)

    if (nu == None):
        return (S_nu1, S_nu2, S_nu3, S_0)
    else:
        c_nu1 = -4*(1-nu)
        c_nu2 = nu/(1-2*nu)
        S = [[None for _ in range(3)] for _ in range(3)]
        for i in range(3):
            for j in range(3):
                S[i][j] = _psh.SHCoeffs.from_array(c_omega*( \
                         c_nu1*c_nu2*S_nu1[i][j].to_array() + \
                         c_nu1*S_nu2[i][j].to_array() + \
                         c_nu2*S_nu3[i][j].to_array() + \
                         S_0[i][j].to_array() ))
        return S

def T_from_S(S_K):
    T_K = [None for _ in range(3)]
    for i in range(3):
        T_K[i] = _psh.SHCoeffs.from_array(VSH1(S_K[i][0])[0].to_array()+VSH1(S_K[i][1])[1].to_array()+VSH1(S_K[i][2])[2].to_array())
    return T_K

def T_mode(l, m, k, mu=0.5, nu=None, c_omega=1.0, r=1.0, shtype='irr', lmax=100, stress=False, recalc=False, etol=1e-8):
    if (nu == None):
        S_nu1, S_nu2, S_nu3, S_0 = S_mode(l, m, k, mu=mu, nu=nu, c_omega=c_omega, r=r, shtype=shtype, lmax=lmax, etol=etol)
        T_nu1 = T_from_S(S_nu1)
        T_nu2 = T_from_S(S_nu2)
        T_nu3 = T_from_S(S_nu3)
        T_0 = T_from_S(S_0)
        if (stress):
            return (T_nu1, T_nu2, T_nu3, T_0, S_nu1, S_nu2, S_nu3, S_0)
        else:
            return (T_nu1, T_nu2, T_nu3, T_0)
    else:
        S_K = S_mode(l, m, k, mu, nu, c_omega=c_omega, r=r, shtype=shtype, lmax=lmax, recalc=recalc, etol=etol)
        T_K = T_from_S(S_K)
        if (stress):
            return (T_K, S_K)
        else:
            return T_K
