import numpy as np
from sympy.physics.wigner import wigner_d
from sympy import matrix2numpy, Ynm, Matrix

# rotation matrices (counter-clockwise)
def Rz(angle):
    return np.array([[np.cos(angle), -np.sin(angle), 0],
                     [np.sin(angle), np.cos(angle), 0],
                     [0, 0, 1]])


def Ry(angle):
    return np.array([[np.cos(angle), 0, np.sin(angle)],
                     [0, 1, 0],
                     [-np.sin(angle), 0, np.cos(angle)]])


def R(alpha, beta, gamma):
    return Rz(alpha) @ Ry(beta) @ Rz(gamma)


def cart2sph(x, y, z):
    '''
    Cartesian to spherical coordinates.

    '''
    r = np.sqrt(x**2 + y**2 + z**2)
    polar = np.arccos(z/r)
    azimuth = np.arctan2(y, x)
    return r, polar, azimuth


def sph2cart(r, polar, azimuth):
    '''
    Spherical to Cartesian coordinates.

    '''
    x = r * np.sin(polar) * np.cos(azimuth)
    y = r * np.sin(polar) * np.sin(azimuth)
    z = r * np.cos(polar)
    return x, y, z


def Ylm_real(l, m, polar, azimuth):
    '''
    Real spherical harmonics (with sign convention in accordance with ABACUS)

    '''
    if m > 0:
        return np.sqrt(2) * complex(Ynm(l, m, polar, azimuth)).real
    elif m < 0:
        return np.sqrt(2) * complex(Ynm(l, -m, polar, azimuth)).imag
    else: # m == 0
        return float(Ynm(l, 0, polar, azimuth))


def Yl_real(l, polar, azimuth, m_order='abacus'):
    '''
    An array of real spherical harmonics for a given l.

    Note
    ----
    ABACUS arranges m in the following order:

    0, 1, -1, 2, -2, ..., l, -l


    '''
    if m_order == 'abacus':
        return np.array([Ylm_real(l, (m+1)//2 * (-1)**(m+1), polar, azimuth) \
                for m in range(0, 2*l+1)])
    elif m_order == 'descend':
        return np.array([Ylm_real(l, m, polar, azimuth) for m in range(l,-l-1,-1)])
    elif m_order == 'ascend':
        return np.array([Ylm_real(l, m, polar, azimuth) for m in range(-l, l+1)])
    else:
        raise ValueError('m_order must be "abacus", "descend" or "ascend"')


def ToReal(l):
    '''
    Standard-to-real transform for spherical harmonics:

            |lm_real> = \sum_{m'} |lm'> U[m',m]

    where m is arranged in descending order (in order to use with sympy's wigner_d)

    '''
    U = np.zeros((2*l+1, 2*l+1), dtype=complex)
    U[l, l] = 1.0
    for m in range(1, l+1):
        U[l-m, l-m] = 1.0 / np.sqrt(2)
        U[l-m, l+m] = -1j / np.sqrt(2)
        U[l+m, l-m] = (-1)**m / np.sqrt(2)
        U[l+m, l+m] = 1j * (-1)**m / np.sqrt(2)

    return U


def wigner_d_real(l, alpha, beta, gamma):
    '''
    Wigner D-matrix for real spherical harmonics.

    This function returns the Wigner D-matrix for real spherical harmonics where
    m is arranged in accordance with ABACUS, and the rotations are counter-clockwise:

    D_real[m',m] = <lm'|exp(-i*Jz*alpha) * exp(-i*Jy*beta) * exp(-i*Jz*gamma))|lm>

    Notes
    -----
    sympy's wigner_d generates the Wigner D-matrix for standard spherical
    harmonics where m is arranged in descending order and the rotations
    are clockwise. See sympy's source code and documentation for details.

    '''
    D = matrix2numpy(wigner_d(l, -alpha, -beta, -gamma), dtype=complex)

    # transform to real spherical harmonics
    U = ToReal(l)
    D_real = np.real(U.T.conj() @ D @ U)

    # rearrange m to the ABACUS order from the descending order
    idx = np.zeros(2*l+1, dtype=int)
    idx[::2] = np.arange(l, 2*l+1)
    idx[1::2] = np.arange(l-1, -1, -1)
    # idx = [l, l-1, l+1, l-2, l+2, ..., 0, 2*l]

    D_real = D_real[:, idx]
    D_real = D_real[idx, :]

    return D_real


########################################################################
#                           Tests
########################################################################

import unittest

class _TestWignerReal(unittest.TestCase):

    def test_yreal(self):
        # random vector
        v = np.random.randn(3)
        r, polar, azimuth = cart2sph(*v)

        C = np.sqrt(0.75/np.pi)
        YReal_ref = [C*v[2]/r, -C*v[0]/r, -C*v[1]/r] # z/r, -x/r, -y/r
        YReal_1 = [Ylm_real(1, m, polar, azimuth) for m in [0, 1, -1]]
        YReal_2 = Yl_real(1, polar, azimuth)

        self.assertTrue(np.allclose(YReal_ref, YReal_1, rtol=1e-12, atol=1e-12))
        self.assertTrue(np.allclose(YReal_ref, YReal_2, rtol=1e-12, atol=1e-12))


    def test_toreal(self):
        l = 3

        # random vector
        v = np.random.randn(3)
        r, polar, azimuth = cart2sph(*v)

        Yl = np.array([complex(Ynm(l, m, polar, azimuth)) for m in range(l, -l-1, -1)])
        Yl_real_ref = Yl_real(l, polar, azimuth, m_order='descend')
        U = ToReal(l)

        # unitarity
        self.assertTrue(np.allclose(U.T.conj() @ U, np.eye(2*l+1), rtol=1e-12, atol=1e-12))

        # Yl_real_ref == Yl @ U
        self.assertTrue(np.allclose(Yl @ U, Yl_real_ref, rtol=1e-12, atol=1e-12))


    def test_wigner_d_real(self):
        l = 3

        # Yl_real along some random vector v
        v = np.random.randn(3)
        _, polar, azimuth = cart2sph(*v)
        Yl_real_v = Yl_real(l, polar, azimuth) # ABACUS order

        # Euler angles
        alpha = np.random.rand() * 2 * np.pi
        beta = np.random.rand() * np.pi
        gamma = np.random.rand() * 2 * np.pi

        # Yl_real along the rotated vector u = R @ v
        u = R(alpha, beta, gamma) @ v
        _, polar, azimuth = cart2sph(*u)
        Yl_real_u = Yl_real(l, polar, azimuth) # ABACUS order

        # Wigner D-matrix for real spherical harmonics with m in ABACUS order
        D_real = wigner_d_real(l, alpha, beta, gamma)

        # Yl_real_v == Yl_real_u @ D_real
        self.assertTrue( np.allclose(Yl_real_v, Yl_real_u @ D_real) )


if __name__ == '__main__':
    unittest.main()
