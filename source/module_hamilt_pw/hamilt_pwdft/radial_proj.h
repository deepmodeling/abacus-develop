/**
 * @file radial_projection.h
 * 
 * project any atom-centered function that has seperatable radial and angular parts
 * or any function can expanded with spherical harmonics onto the planewave basis,
 * although the latter will be somewhat cumbersome:
 * f(r) = sum_{l,m} f_{lm}(r) * Ylm(theta, phi)
 * F(q) = sum_{l,m} i^l * 4*pi/sqrt(omega) * Jl[f_{lm}](q) * Ylm(q)
 */

#ifndef RADIAL_PROJECTION_H
#define RADIAL_PROJECTION_H

#include "module_base/vector3.h"
#include "module_base/cubic_spline.h"
#include <memory>
#include <vector>
#include <complex>
#include <map>
#include <utility>
#include <algorithm>

#include "module_cell/unitcell.h"
#include "module_psi/psi.h"
#include "module_basis/module_pw/pw_basis_k.h"

namespace RadialProjection
{
    /**
     * @brief RadialProjector is for projecting a function who has seperatable radial
     * and angular parts:
     * f(r) = f(|r|) * Ylm(theta, phi)
     * onto planewave basis. 
     * 
     * Usage:
     * 
     * 1. classical way: reciprocal space integration
     * ```c++
     * RadialProjector rp(false);  
     * const int nq = 1000;  
     * const double dq = 0.01;  
     * // given `r` is the real space grid and `radials` is the collection of radial
     * // functions, `l` is the angular momentum quantum number for each radial function
     * // then the interpolation table can be rapidly built by calling SphericalBesselTransformer
     * // and CubicSpline modules.
     * rp._build_sbt_tab(r, radials, l, nq, dq);
     * // then the set of q will used to calculate the Fourier transform
     * rp.sbtft(qs, out, 'r', omega, tpiba);
     * // in `out`, there will be the Fourier transform of the radial functions organized
     * // in the same way as the input `radials` and `qs`, as row and column respectively.
     * // but one should note for each radials, there are 2*l+1 components now instead of
     * // just one.
     * 
     * // one may find it is not easy to maintain such a large table, so here, also provides
     * // a tool function to map the 2D index to 1D index, and vice versa. With the angular
     * // momentum used in function _build_sbt_tab, one can easily build a map from 
     * // [irad][im] to 1D index, and use two functions _irad_m_to_idx and _idx_to_irad_m
     * // to convert between 1D index and [irad][m] (instead of im!).
     * std::vector<std::vector<int>> map_;
     * rp._build_sbtft_map(l, map_);
     * ```
     * 
     * 2. SBFFT: small box fast-fourier-transform (not implemented yet)
     */
    class RadialProjector
    {
        public:
            /**
             * Notation of following two functions:
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
             * irow2ia: for given `irow`, the index of row, return the `ia`: the index of atom within the range
             * of `it`.
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
            
            static void _build_backward_map(const std::vector<std::vector<int>>& it2iproj,
                                            const std::vector<int>& iproj2l,
                                            std::vector<int>& irow2it,
                                            std::vector<int>& irow2iproj,
                                            std::vector<int>& irow2m);
            static void _build_forward_map(const std::vector<std::vector<int>>& it2ia,
                                           const std::vector<std::vector<int>>& it2iproj,
                                           const std::vector<int>& iproj2l,
                                           std::map<std::tuple<int, int, int, int>, int>& itiaiprojm2irow);
        public:
            /**
             * @brief Construct a new Radial Projector object
             * 
             * @param realspace if perform integration in real space rather than reciprocal space
             * , default is false
             * 
             * @attention Currectly only reciprocal space method is implemented
             */
            RadialProjector(const bool realspace = false) {}
            ~RadialProjector() {}

            // it is more feasible to build interpolation table. this function will tabulate
            // for those functions. This function will write the in-build tab_, to place the
            // values of Jl[f](q) for each q and l.
            // Here I provide two versions of tabulate, one for the case may be capable to 
            // avoid the memory copy operation, and the other one is for the case of 
            // std::vector<double> and std::vector<std::vector<double>>.
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
            void _build_sbt_tab(const int nr,
                                const double* r,
                                const std::vector<double*>& radials,
                                const std::vector<int>& l,
                                const int nq,                             //< GlobalV::DQ
                                const double& dq);                        //< GlobalV::NQX
            void _build_sbt_tab(const std::vector<double>& r,
                                const std::vector<std::vector<double>>& radials,
                                const std::vector<int>& l,
                                const int nq,                             //< GlobalV::DQ
                                const double& dq);                        //< GlobalV::NQX
            // compatibility concern: for FS_Nonlocal_tools. Will not call sbtft so need omega
            void _build_sbt_tab(const std::vector<int>& nproj,
                                const std::vector<double>& r,
                                const std::vector<std::vector<double>>& radials,
                                const std::vector<int>& l,
                                const int nq,                             //< GlobalV::DQ
                                const double& dq,                         //< GlobalV::NQX
                                const double& omega,
                                const int npol,
                                ModuleBase::realArray& tab,
                                ModuleBase::matrix& nhtol);
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
                       const char type = 'r',                                   // 'r' for ket |>, 'l' for bra <|
                       const double& omega = 1.0,
                       const double& tpiba = 1.0);                                  // 'n' for no gradient, 'x', 'y', 'z' for gradient in x, y, z direction
            
            void sbfft(); // interface for SBFFT


        private:
            std::unique_ptr<ModuleBase::CubicSpline> cubspl_;
            std::vector<int> l_;
    };
  
    /** ====================================================================================
     * 
     *                       Small box Fast-Fourier-Transform (SBFFT)
     * 
     * ====================================================================================
     * Small box FFT is a technique for quickly intergrating real-space localized functions
     * , or say perform FFT in a small box defined by a "mask function" with relatively low 
     * time complexity, will be out-performing in system where number of atoms are larger 
     * than ten. For details please refer to the work: 
     * Mask-function real-space implementations of nonlocal pseudopotentials
     * by Wang, L.-W., PHYSICAL REVIEW B, VOLUME64,201107(R)
     *
     * Following are the brief technical review of this technique. Given the function to
     * be transformed w(r):
     * 1. Generate the q-grid in range |q| < 2qmax - qc, in which the qmax is the one
     *    defined by ecutrho, and qc is the cutoff defined by ecutwfc.
     * 2. With mask function m(r) generated, make division of w(r) by m(r). The real space
     *    cutoff (or the radius that w(r) vanishes), r0, must be smaller than the cutoff
     *    of m(r) to ensure the convergence of the division. Denote wm(r) = w(r)/m(r).
     * 3. Make FT on wm(r) to get wm(q)
     * 4. Make inverse FT with only the q-grid in range |q| < 2qmax - qc to get the
     *    wm'(r).
     * 5. Perform real-space integration on function w'(r)*m(r)*exp(iqr).
     */

    /**
     * @brief get the mask function for SBFFT
     * 
     * @param mask mask function
     */
    void _mask_func(std::vector<double>& mask);

    /**
     * @brief do operation w(r)/m(r) on a radial function. The cutoff radius of w(r) 
     * is smaller than the cutoff radius of m(r). The m(r) has been rescaled so that 
     * r ranges from 0 to 1.
     * 
     * @param nr1 number of grid points of function to operate
     * @param r grid points of function to operate
     * @param in function to operate
     * @param nr2 number of grid points of mask function
     * @param mask mask function
     * @param out output value
     */
    void _do_mask_on_radial(const int nr1,
                            const double* r,
                            const double* in,
                            const int nr2,
                            const double* mask,
                            double* out);
}

#endif // RADIAL_PROJECTION_H