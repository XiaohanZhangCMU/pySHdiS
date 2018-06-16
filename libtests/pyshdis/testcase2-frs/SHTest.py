import numpy as _np
import pyshtools as _psh
import scipy as _sp

def second_deriv_R( X, Y, Z ):
# Rij = 
# [ (y^2 + z^2)/r^3,      -(x*y)/r^3,      -(x*z)/r^3]
# [      -(x*y)/r^3, (x^2 + z^2)/r^3,      -(y*z)/r^3]
# [      -(x*z)/r^3,      -(y*z)/r^3, (x^2 + y^2)/r^3]

    R = _np.sqrt(X**2 + Y**2 + Z**2)
    B = [[(Y**2+Z**2)/R**3, -X*Y/R**3, -X*Z/R**3], \
         [-X*Y/R**3, (X**2+Z**2)/R**3, -Y*Z/R**3], \
         [-X*Z/R**3, -Y*Z/R**3, (X**2+Y**2)/R**3]]
    return B

def third_deriv_R(X, Y, Z):
    # B(ind,:) = int( d_ijk R, s=-inf..Z1)
    # where ind = (i-1)^2+(j-1)^2+(k-1)^2

    # sYms X Y Z a
    # R = sqrt(X^2+Y^2+Z^2+a^2)
    R = _np.sqrt(X**2 + Y**2 + Z**2)
    B = [_np.array([]) for X in range(13)]

    # X X X: [i-1 j-1 k-1]^2 = [0 0 0], ind-1 = 0
    # diff(diff(diff(R,X),X),X)
    # (3*X^3)/(a^2 + X^2 + Y^2 + Z^2)^(5/2) - (3*X)/(a^2 + X^2 + Y^2 + Z^2)^(3/2)
    B[0] = (3*X**3)/(R**5) - (3*X)/(R**3)

    # X X Y: [i-1 j-1 k-1]^2 = [0 0 1], ind-1 = 1
    #diff(diff(diff(R,X),X),Y)
    #(3*X^2*Y)/(a^2 + X^2 + Y^2 + Z^2)^(5/2) - Y/(a^2 + X^2 + Y^2 + Z^2)^(3/2)
    B[1] =  (3*X**2*Y)/(R**5) - Y/(R**3)

    # X Y Y: [i-1 j-1 k-1]^2 = [0 1 1], ind-1 = 2
    #diff(diff(diff(R,X),Y),Y)
    #(3*X*Y^2)/(a^2 + X^2 + Y^2 + Z^2)^(5/2) - X/(a^2 + X^2 + Y^2 + Z^2)^(3/2)
    B[2] = (3*X*Y**2)/(R**5) - X/(R**3)

    # Y Y Y: [i-1 j-1 k-1]^2 = [1 1 1], ind-1 = 3
    #diff(diff(diff(R,Y),Y),Y)
    #(3*Y^3)/(a^2 + X^2 + Y^2 + Z^2)^(5/2) - (3*Y)/(a^2 + X^2 + Y^2 + Z^2)^(3/2)
    B[3] = (3*Y**3)/(R**5) - (3*Y)/(R**3)

    # X X Z: [i-1 j-1 k-1]^2 = [0 0 4], ind-1 = 4
    #diff(diff(diff(R,X),X),Z)
    #
    B[4] = (3*X**2*Z)/(R**5) - Z/(R**3)

    # X Y Z: [i-1 j-1 k-1]^2 = [0 1 4], ind-1 = 5
    #diff(diff(diff(R,X),Y),Z)
    #(3*X*Y*Z)/(a^2 + X^2 + Y^2 + Z^2)^(5/2)
    B[5] = (3*X*Y*Z)/(R**5)

    # Y Y Z: [i-1 j-1 k-1]^2 = [1 1 4], ind-1 = 6
    #diff(diff(diff(R,Y),Y),Z)
    #
    B[6] = (3*Y**2*Z)/(R**5) - Z/(R**3)

    # X Z Z: [i-1 j-1 k-1]^2 = [0 4 4], ind-1 = 8
    #diff(diff(diff(R,X),Z),Z)
    #
    B[8] = (3*X*Z**2)/(R**5) - X/(R**3)

    # Y Z Z: [i-1 j-1 k-1]^2 = [1 4 4], ind-1 = 9
    #diff(diff(diff(R,Y),Z),Z)
    #
    B[9] = (3*Y*Z**2)/(R**5) - Y/(R**3)

    # Z Z Z: [i-1 j-1 k-1]^2 = [4 4 4], ind-1 = 12
    #diff(diff(diff(R,Z),Z),Z)
    #
    B[12] = (3*Z**3)/(R**5) - (3*Z)/(R**3)

    return B

######################### Analytical solution from Takahashi & Ghoniem (2008) ################

def Legendre_poly(N, Z, LegendreP,dl,dm, csphase=-1):
    nmax = N.max()
    dl = _np.array(dl)
    dm = _np.array(dm)
    p = _np.empty(Z.shape+dl.shape)
    for idx, _ in _np.ndenumerate(N):
        n, z = (N[idx], Z[idx])
        Pvalue = LegendreP(nmax+3, z, csphase=csphase)
        l = n + dl; m = n + dm;
        Pidx = _np.array((l*(l+1)/2+m), dtype=_np.int)
        p[idx] = Pvalue[Pidx]

    return p

def void_screw_disl(nmax, zs, ts, mu1, mu2, nu1=0.25, nu2=0.25, a=1, b=1):
    ns = _np.arange(nmax) + 1
    n, z, t = _np.meshgrid(ns, zs, ts)
    r = _np.sqrt(z**2 + t**2)
    z = z/r
    lambda1 = 2*mu1*nu1/(1-2*nu1)
    lambda2 = 2*mu2*nu2/(1-2*nu2)

    Kn = -(lambda1+mu1)/2/((n+2)*lambda1+(3*n+5)*mu1)
    EnI = 1/(2*n+1)*((n-1)*lambda1 - (n+4)*mu1)/((n+2)*lambda1+(3*n+5)*mu1)
    eta_n = b/4/_np.pi * ((-2)**n) *_sp.misc.factorial(n-1)/_sp.misc.factorial(2*n-1)

    Cab = (mu1-mu2)*eta_n/(mu1*(n+2)+mu2*(n-1))
    denom_ab = mu1*((n+2) + EnI*(2*n+1)*(n+1)) + mu2*2*n
    a_t_n = (a/t)**n
    alpha_n = Cab/denom_ab * (mu1*((n+2) - EnI*(2*n+1)*(n-1))) * a_t_n
    beta_n = -Cab/denom_ab * (mu1*((n+2) + EnI*(2*n+1)*(n-1)) + 2*mu2*(n-1)) * a_t_n

    Omega_n = 2*eta_n*(mu1-mu2)* a_t_n * (a**(n+1)) / denom_ab

    C_term1 = (-a**(n+1) * (2*mu1*alpha_n + mu1*beta_n) + 6*mu1*Kn*Omega_n*(2*n+1)/(2*n+5)) /r**(n+2)
    C_term2 = (-a**(n+1)/2 *mu1*beta_n + mu1*Kn*Omega_n*(2*n+3)/(2*n+5)) /r**(n+2)
    C_term3 = (2*mu1*(a**2/r**2 - 1)*Kn*Omega_n+mu1*Kn*Omega_n*2/(2*n+5)) /r**(n+2)
    C_term4 = (24*mu1*(a**2/r**2 - 1)*Kn*Omega_n+mu1*Kn*Omega_n*24/(2*n+5)) /r**(n+2)

    # After experiment, the implementation in Takahashi&Ghoniem(2008) is unnormalized Legendre Polynomial
    p = Legendre_poly(n, z, _psh.legendre.PLegendreA,dl = _np.array([1, 1, 3, 3]),dm = _np.array([-1,1, 1,-1]), csphase=-1)
    return _np.sum(C_term1*p[:,:,:,0] + C_term2*p[:,:,:,1] + C_term3*p[:,:,:,2] + C_term4*p[:,:,:,3], axis=1)