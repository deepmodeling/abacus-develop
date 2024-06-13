#include "module_io/numerical_basis_jyjy.h"

#include "module_base/matrix3.h"
#include "module_base/vector3.h"
#include "module_basis/module_nao/two_center_integrator.h"

namespace NumericalBasis
{
#ifdef __LCAO
using index_t = std::tuple<int, int, int, int>;
using Vec3 = ModuleBase::Vector3<double>;
using ModuleBase::cross;
using ModuleBase::dot;

std::vector<index_t> indexgen(const std::vector<int>& natom, const std::vector<int>& lmax)
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
                for (int M = 0; M < 2 * l + 1; ++M)
                {
                    // convert the "abacus M" to the conventional m
                    const int m = (M % 2 == 0) ? -M / 2 : (M + 1) / 2;
                    index.emplace_back(itype, iatom, l, m);
                }
            }
        }
    }
    return index;
}

ModuleBase::ComplexArray cal_overlap_Sq(const char type, const int lmax, const int nbes, const double rcut,
                                        const std::vector<std::vector<Vec3>>& R, const ModuleBase::Matrix3& latvec,
                                        const std::vector<index_t>& mu_index)
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
            std::vector<Vec3> dR_all = neighbor_vec(R[t2][a2] - R[t1][a1], latvec, 2 * rcut);
            for (auto& dR: dR_all)
            {
                for (int zeta1 = 0; zeta1 < nbes; ++zeta1)
                {
                    for (int zeta2 = 0; zeta2 < nbes; ++zeta2)
                    {
                        double elem = 0.0;
                        intor.calculate(t1, l1, zeta1, m1, t2, l2, zeta2, m2, dR, &elem);
                        overlap_Sq(it1 - mu_index.cbegin(), it2 - mu_index.cbegin(), zeta1, zeta2) += elem;
                    }
                }
            }
        }
    }

    return overlap_Sq;
}

std::vector<Vec3> neighbor_vec(const Vec3& d0, const ModuleBase::Matrix3& latvec, const double r)
{
    // convert the format of lattice vectors from Matrix3 to a vector of
    // Vector3<double> before searching for neighboring periodic images
    std::vector<Vec3> a{{latvec.e11, latvec.e12, latvec.e13},
                        {latvec.e21, latvec.e22, latvec.e23},
                        {latvec.e31, latvec.e32, latvec.e33}};

    // unit vectors normal to the three faces of the parallelepiped
    std::vector<Vec3> u(3);
    for (int i = 0; i < 3; ++i)
    {
        u[i] = cross(a[(i + 1) % 3], a[(i + 2) % 3]);
        u[i] = u[i] / u[i].norm();
    }

    // range of supercell coordinates
    //
    // Suppose n0, n1, n2 can take any real numbers (rather than just integers),
    // given any n0, one can always find n1 & n2 such that
    //
    //          d0 + n0*a0 + n1*a1 + n2*a2
    //
    // is perpendicular to the plane spanned by a1 and a2. Therefore, the
    // min/max values of n0 should satisfy
    //
    //          dot(d0 + n0*a0, u0) = r
    //
    std::vector<int> n_min(3), n_max(3);
    for (int i = 0; i < 3; ++i)
    {
        double n_min_tmp = (-r - dot(d0, u[i])) / dot(a[i], u[i]);
        double n_max_tmp = (r - dot(d0, u[i])) / dot(a[i], u[i]);
        if (n_min_tmp > n_max_tmp) // in case dot(a[i], u[i]) < 0
        {
            std::swap(n_min_tmp, n_max_tmp);
        }
        // cast to int truncates towards 0
        // slightly widens the range for floating point precision error
        n_min[i] = static_cast<int>(n_min_tmp - 1e-6);
        n_max[i] = static_cast<int>(n_max_tmp + 1e-6);
    }

    // loop over the neighboring supercells
    std::vector<Vec3> d;
    for (int n0 = n_min[0]; n0 <= n_max[0]; ++n0)
    {
        for (int n1 = n_min[1]; n1 <= n_max[1]; ++n1)
        {
            for (int n2 = n_min[2]; n2 <= n_max[2]; ++n2)
            {
                Vec3 d_tmp = d0 + static_cast<double>(n0) * a[0] + static_cast<double>(n1) * a[1]
                             + static_cast<double>(n2) * a[2];

                if (d_tmp.norm() < r)
                {
                    d.push_back(d_tmp);
                }
            }
        }
    }

    return d;
}

#endif

} // end of namespace NumericalBasis
