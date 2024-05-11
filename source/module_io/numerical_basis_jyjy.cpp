#include "module_io/numerical_basis_jyjy.h"
#include "module_basis/module_nao/two_center_integrator.h"

namespace NumericalBasis
{
#ifdef __LCAO
using index_t = std::tuple<int, int, int, int>;

std::vector<index_t> indexgen(const std::vector<int>& natom,
                              const std::vector<int>& lmax)
{
#ifdef __DEBUG
    assert(natom.size() == lmax.size());
#endif
    std::vector<index_t> index;
    for (size_t itype = 0; itype < natom.size(); ++itype)
    {
        for (int iatom = 0; iatom < natom[itype]; ++iatom)
        {
            for (int l = 0; l <= lmax[itype]; ++l)
            {
                for (int M = 0; M < 2*l+1; ++M)
                {
                    // convert the "abacus M" to the conventional m
                    const int m = (M % 2 == 0) ? -M/2 : (M+1)/2;
                    index.emplace_back(itype, iatom, l, m);
                }
            }
        }
    }
    return index;
}


ModuleBase::ComplexArray cal_overlap_Sq(
    const char type,
    const int lmax,
    const int nbes,
    const double rcut,
    const std::vector<std::vector<ModuleBase::Vector3<double>>>& R,
    const std::vector<index_t>& mu_index
)
{
    // allocate output array
    const int nlocal = mu_index.size();
    ModuleBase::ComplexArray overlap_Sq(nlocal, nlocal, nbes, nbes);
    overlap_Sq.zero_out();

    // build a RadialCollection of spherical Bessel functions
    const double dr = 0.005; // grid spacing for SphbesRadials; smaller for higher precision
    RadialCollection orb;
    orb.build(lmax, nbes, rcut, 0.0, dr);

    ModuleBase::SphericalBesselTransformer sbt;
    orb.set_transformer(sbt);

    const double rmax = orb.rcut_max() * 2.0;
    const int nr = static_cast<int>(rmax / dr) + 1;
    
    orb.set_uniform_grid(true, nr, rmax, 'i', true);

    // build the two-center integrator
    TwoCenterIntegrator intor;
    intor.tabulate(orb, orb, type, nr, rmax);

    // traverse the vector of composite index (itype, iatom, l, m)
    int t1 = 0, a1 = 0, l1 = 0, m1 = 0;
    int t2 = 0, a2 = 0, l2 = 0, m2 = 0;
    for (auto it1 = mu_index.cbegin(); it1 != mu_index.cend(); ++it1)
    {
        std::tie(t1, a1, l1, m1) = *it1;
        for (auto it2 = mu_index.cbegin(); it2 != mu_index.cend(); ++it2)
        {
            std::tie(t2, a2, l2, m2) = *it2;
            ModuleBase::Vector3<double> dR = R[t2][a2] - R[t1][a1];
            for (int zeta1 = 0; zeta1 < nbes; ++zeta1)
            {
                for (int zeta2 = 0; zeta2 < nbes; ++zeta2)
                {
                    double elem = 0.0;
                    intor.calculate(t1, l1, zeta1, m1,
                                    t2, l2, zeta2, m2,
                                    dR, &elem);
                    overlap_Sq(it1 - mu_index.begin(),
                               it2 - mu_index.begin(),
                               zeta1, zeta2) = elem;
                }
            }
        }
    }

    return overlap_Sq;
}
#endif

} // end of namespace NumericalBasis
