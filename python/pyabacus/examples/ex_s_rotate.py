from __future__ import annotations

# import pytest
from pyabacus import ModuleNAO as nao
from pyabacus import ModuleBase as base
import numpy as np

import copy
from indexmap import _index_map
from wigner_real import R, wigner_d_real

import matplotlib.pyplot as plt

orb_dir = '../../../tests/PP_ORB/'
file_list = ['H_gga_8au_100Ry_2s1p.orb', 'O_gga_10au_100Ry_2s2p1d.orb']

file_list = [orb_dir + orbfile for orbfile in file_list ]
nfile = len(file_list)

orb = nao.RadialCollection()
orb.build(nfile, file_list, 'o')
sbt = base.SphericalBesselTransformer()
orb.set_transformer(sbt)

rmax = orb.rcut_max * 2.0
dr = 0.01
nr = int(rmax / dr) + 1

orb.set_uniform_grid(True, nr, rmax, 'i', True)
S_intor = nao.TwoCenterIntegrator()
S_intor.tabulate(orb, orb, 'S', nr, rmax)
vR = np.array([1.0,0,0]).astype(np.float64)
#out = S_intor.calculate(1, 1, 0, 1, 0, 1, 0, 1, vR)[0]
#print(out)

comp2mu, mu2comp = _index_map(2, [2, 1], [1, 2], [[2,1], [2,2,1]])
#print(comp2mu, mu2comp)

nao = len(comp2mu.items())

coord = [
        [
            np.array([-12.081531451316582,16.463368531712373,10.304287878967891]),
            np.array([ -12.056180479123837,19.25408045699522,10.010554611214044])
            ],
        [
            np.array([-13.1930046246741,17.91132430713516,10.440289103003526])
            ]
        ]


S = np.zeros((nao, nao))
for mu in range(nao):
    itype1, iatom1, l1, zeta1, m1 = mu2comp[mu]
    for nu in range(nao):
        itype2, iatom2, l2, zeta2, m2 = mu2comp[nu]
        vR = coord[itype2][iatom2] - coord[itype1][iatom1]
        S[mu, nu] = S_intor.calculate(itype1, l1, zeta1, m1, itype2, l2, zeta2, m2, vR)[0]


S_prerot = np.copy(S)
coord_pre = copy.deepcopy(coord)

val_pre, vec_pre = np.linalg.eigh(S_prerot)

# random rotation
alpha = np.random.rand() * 2 * np.pi
beta = np.random.rand() * np.pi
gamma = np.random.rand() * 2 * np.pi
Rmat = R(alpha, beta, gamma)

for coord_type in coord:
    for i in range(len(coord_type)):
        coord_type[i] = Rmat @ coord_type[i]

D = np.zeros((nao, nao))
mu = 0
while mu < nao:
    itype, iatom, l, zeta, m = mu2comp[mu]
    D[mu:mu+2*l+1, mu:mu+2*l+1] = wigner_d_real(l, alpha, beta, gamma)
    mu += 2*l+1

S_e3 = D @ S_prerot @ D.T # overlap matrix after rotation (by angular momentum theory)

# numerical evaluation of S_e3 by two-center integrals
S = np.zeros((nao, nao))
for mu in range(nao):
    itype1, iatom1, l1, zeta1, m1 = mu2comp[mu]
    for nu in range(nao):
        itype2, iatom2, l2, zeta2, m2 = mu2comp[nu]
        vR = coord[itype2][iatom2] - coord[itype1][iatom1]
        S[mu, nu] = S_intor.calculate(itype1, l1, zeta1, m1, itype2, l2, zeta2, m2, vR)[0]

val, vec = np.linalg.eigh(S)

#print('S eigval diff: ', np.linalg.norm(val-val_pre, np.inf))
print('norm(S_e3 - S_numer) = ', np.linalg.norm(S_e3 - S, np.inf))


