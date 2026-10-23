#ifndef RADIAL_PROJECTION_H
#define RADIAL_PROJECTION_H

/**
 * @file radial_proj.h
 *
 * project any atom-centered function that has separable radial and angular parts
 * or any function that can be expanded with spherical harmonics onto the planewave basis,
 * although the latter will be somewhat cumbersome:
 * f(r) = sum_{l,m} f_{lm}(r) * Ylm(theta, phi)
 * F(q) = sum_{l,m} i^l * 4*pi/sqrt(omega) * Jl[f_{lm}](q) * Ylm(q)
 */

#include "source_base/vector3.h"
#include "source_base/cubic_spline.h"
#include "source_base/matrix.h"
#include "source_base/realarray.h"
#include <memory>
#include <vector>
#include <complex>
#include <map>
#include <tuple>

namespace RadialProjection
{
    /**
     * Notation of the following two functions:
     *
     * Given all the projectors are listed in a series, so the `iproj` is the index goes across
     * all atomtypes, which means if for the first type, the iproj goes from 0 to 4, then the
     * second atomtypes the iproj will start from 5, and so on...
     * However, there is also another convention, like numerical atomic orbitals, developer always
     * use "l" to index orbitals, here, in all output map, the `iproj` will start from 0, which
     * means in output the `iproj` is local index.
     * -----------------------------------------------------------------------------------------
     * First, the following lists should be prepared as early as possible,
     *
     * it2iproj: for given it, the index of atom type, return the list of index of projectors.
     *
     * iproj2l: for given iproj, the index of projectors, return the l of this projector. More
     * simply explaning, it is just the list of angular momentum of projectors.
     *
     * it2ia: just a list that stolen information from UnitCell, for given it, the index of atom
     * within the range of it. So this list is different from the it2iproj, iproj is the index
     * across type but ia is the index within the type. So for each it2ia[it], the ia, in principle
     * , always/can start from 0.
     *
     * One may question that does the indexing support one atom type with multiple projectors? The
     * answer is YES. Combining the it2iproj and it2ia, one can even support PART of atoms of one
     * type has multiple projectors.
     * -----------------------------------------------------------------------------------------
     * Then the returned lists,
     *
     * irow2it: for given `irow`, the index of row, return the `it`: the index of atom type.
     *
     * irow2iproj: for given `irow`, the index of row, return the `iproj`, the index of projectors,
     * note that this `iproj` is the local index.
     *
     * irow2m: for given irow, the index of row, return the m, the magnetic quantum number of this
     * projector.
     *
     * One may complain that cannot get `l` from the `irow`, but the truth is, not exactly. One can
     * get the `l` starting from `irow` by:
     * ```c++
     * const int iproj = irow2iproj[irow];
     * const int it = irow2it[irow];
     * const int iproj_g = it2iproj[it][iproj];
     * const int l = iproj2l[iproj_g];
     * ```
     */
    void build_backward_map(const std::vector<std::vector<int>>& it2iproj,
                            const std::vector<int>& iproj2l,
                            std::vector<int>& irow2it,
                            std::vector<int>& irow2iproj,
                            std::vector<int>& irow2m);

    void build_forward_map(const std::vector<std::vector<int>>& it2ia,
                           const std::vector<std::vector<int>>& it2iproj,
                           const std::vector<int>& iproj2l,
                           std::map<std::tuple<int, int, int, int>, int>& itiaiprojm2irow);

    /**
     * @brief make interpolation tables for the Spherical Bessel Transform of
     * type-wise radial projectors, in the (ntype, nprojmax*npol, nq) realArray
     * layout used by the nonlocal-operator kernels, plus the (it, ih) -> l map.
     *
     * @param nproj number of projectors for each atom type
     * @param r radial grids, shared by all radial functions
     * @param radials radial functions, each element is a radial function
     * @param l angular momentum quantum number for each radial function
     * @param nq number of q-points
     * @param dq space between q-points
     * @param omega cell volume, used in the prefactor 4*pi/sqrt(omega)
     * @param npol number of spinor components (for nspin 4)
     * @param tab [out] interpolation table, (ntype, nprojmax*npol, nq)
     * @param nhtol [out] map from (it, ih) to l, with ih the (l, m)-distinctive index
     */
    void build_sbt_tab(const std::vector<int>& nproj,
                       const std::vector<double>& r,
                       const std::vector<std::vector<double>>& radials,
                       const std::vector<int>& l,
                       const int nq,
                       const double dq,
                       const double omega,
                       const int npol,
                       ModuleBase::realArray& tab,
                       ModuleBase::matrix& nhtol);

    /**
     * @brief RadialProjector holds the interpolation table of the Spherical Bessel
     * Transform of a set of radial functions and evaluates the analytical Fourier
     * transform on arbitrary q vectors.
     *
     * Usage:
     *
     * reciprocal space integration
     * ```c++
     * RadialProjector rp;
     * const int nq = 1000;
     * const double dq = 0.01;
     * // given `r` is the real space grid and `radials` is the collection of radial
     * // functions, `l` is the angular momentum quantum number for each radial function
     * // then the interpolation table can be rapidly built by calling SphericalBesselTransformer
     * // and CubicSpline modules.
     * rp.build_sbt_tab(r, radials, l, nq, dq);
     * // then the set of q will used to calculate the Fourier transform
     * rp.sbtft(qs, out, 'r', omega, tpiba);
     * // in `out`, there will be the Fourier transform of the radial functions organized
     * // in the same way as the input `radials` and `qs`, as row and column respectively.
     * // but one should note for each radials, there are 2*l+1 components now instead of
     * // just one.
     * ```
     */
    class RadialProjector
    {
        public:
            RadialProjector() = default;
            ~RadialProjector() = default;

            /**
             * @brief make a interpolation table for the Spherical Bessel Transform of f(r)
             *
             * @param nr number of grid points, shared by all radial functions
             * @param r radial grids, shared by all radial functions
             * @param radials radial functions, each element is a radial function
             * @param l angular momentum quantum number for each radial function
             * @param nq number of q-points
             * @param dq space between q-points
             */
            void build_sbt_tab(const int nr,
                               const double* r,
                               const std::vector<double*>& radials,
                               const std::vector<int>& l,
                               const int nq,
                               const double dq);
            void build_sbt_tab(const std::vector<double>& r,
                               const std::vector<std::vector<double>>& radials,
                               const std::vector<int>& l,
                               const int nq,
                               const double dq);
            /**
             * @brief perform analytical version of the Fourier transform:
             * F(q) = int(f(r)*exp(-iq.r) d^3r)
             *      = 4*pi/sqrt(omega) * (-i)^l * Jl[f](q) * Ylm(q)
             * , where Ylm(q) is real spherical harmonic function, and Jl[f](q) is
             * the Spherial Bessel Transform of f(r):
             * Jl[f](q) = int(f(r)*j_l(q*r)*r^2 dr)
             * , where j_l(q*r) is the spherical Bessel function of the first kind.
             * . If use another notation, F(q) = <q|f>, this is denoted as type
             * "r" for ket |>, and "l" for bra <|.
             */

            void sbtft(const std::vector<ModuleBase::Vector3<double>>& qs,
                       std::vector<std::complex<double>>& out,
                       const char type = 'r',            // 'r' for ket |>, 'l' for bra <|
                       const double omega = 1.0,
                       const double tpiba = 1.0);

        private:
            std::unique_ptr<ModuleBase::CubicSpline> cubspl_;
            std::vector<int> l_;
    };

    /**
     * @brief get the mask function for SBFFT
     *
     * @param mask mask function
     */
    void mask_func(std::vector<double>& mask);
}

#endif // RADIAL_PROJECTION_H
