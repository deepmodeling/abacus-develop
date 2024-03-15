def _index_map(ntype, natom, lmax, nzeta=None):
    '''
    Bijective map between the composite index (itype, iatom, l, zeta, m)
    and linearized orbital index mu.

    An atomic orbital is labeled by its type & atomic index, its angular
    momentum quantum numbers l & m, and possibly a zeta number. Suppose
    there's a total of N orbitals, each orbital can also be assigned a
    unique index mu \in [0, N-1].

    This function returns a pair of bijective maps between the composite
    index (itype, iatom, l, zeta, m) and the linearized orbital index mu.
    The composite index is linearized in C-style (lexicographic order).

    Parameters
    ----------
        ntype : int
            Number of atom types.
        natom : list of int
            Number of atoms for each type.
        lmax : list of int
            lmax[i] specifies the maximum angular momentum of type i.
        nzeta : list of list of int
            nzeta[i][l] specifies the number of zeta orbitals of the
            angular momentum l of type i.
            If None, nzeta is assumed to be 1 for all.

    Returns
    -------
        comp2mu : dict
            A dict that maps a composite index (itype, iatom, l, zeta, m)
            to its linearized index.
        mu2comp : dict
            A dict that maps a linearized index to its composite index.

    Notes
    -----
    The linearized index arranges m in accordance with ABACUS:

        0, 1, -1, 2, -2, 3, -3, ..., l, -l

    '''
    if nzeta is None:
        nzeta = [[1]*(lmax[itype]+1) for itype in range(ntype)]

    assert len(natom) == ntype
    assert len(lmax) == ntype
    assert lmax == [len(nzeta[itype])-1 for itype in range(ntype)]

    comp2mu = {}
    mu2comp = {}
    mu = 0
    for itype in range(ntype):
        for iatom in range(natom[itype]):
            for l in range(lmax[itype]+1):
                for zeta in range(nzeta[itype][l]):
                    '''
                    In ABACUS, real spherical harmonics Ylm arranges its m
                    in the following order:

                                  0, 1, -1, 2, -2, 3, -3, ..., l, -l
                    
                    (see module_base/ylm.cpp and module_base/math_ylmreal.cpp
                    for details)

                    '''
                    for mm in range(0, 2*l+1):
                        m = -mm // 2 if mm % 2 == 0 else (mm + 1) // 2
                        comp2mu[(itype, iatom, l, zeta, m)] = mu
                        mu2comp[mu] = (itype, iatom, l, zeta, m)
                        mu += 1

    return comp2mu, mu2comp


############################################################
#                           Test
############################################################
import unittest

class _TestIndexMap(unittest.TestCase):

    def test_index_map(self):
        ntype = 3
        natom = [2, 1, 3]
        lmax = [1, 2, 4]
        nzeta = [[2,3], [1,0,1], [1, 2, 2, 1, 3]]
        comp2mu, mu2comp = _index_map(ntype, natom, lmax, nzeta)

        # check the total number of orbitals
        nao = sum(sum( (2*l+1) * nzeta[itype][l] for l in range(lmax[itype]+1) ) \
                * natom[itype] for itype in range(ntype))
        self.assertEqual( len(mu2comp.items()), nao )

        # check bijectivity
        for mu in range(nao):
            self.assertEqual( comp2mu[mu2comp[mu]], mu )

        # check the first and the last
        self.assertEqual( mu2comp[0], (0, 0, 0, 0, 0) )
        self.assertEqual( mu2comp[nao-1], \
                (ntype-1, natom[-1]-1, lmax[-1], nzeta[-1][-1]-1, -lmax[-1]) )


if __name__ == '__main__':
    unittest.main()

