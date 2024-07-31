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
// project any atom-centered function that has seperatable radial and angular parts
// onto the planewave basis
namespace RadialProjection
{
    /**
     * @brief RadialProjector is for projecting function who has seperatable radial
     * and angular parts:
     * f(r) = f(|r|) * Ylm(theta, phi)
     * onto planewave basis. Each instance of RadialProjector can project any of
     * this function onto a fixed set of q-points (or say one well-defined G-grid
     * points)
     * 
     * This class is designed based on the fact that ABACUS will always calculate
     * all Ylm values from 0 to lmax in-one-shot, therefore it is not wise to always
     * recalculate it, instead this class will cache at the very beginning, with one
     * specified lmax.
     * Additionally, the SphericalBesselTransformer is also cached at the very
     * beginning.
     * 
     * The "language" here will be always: "do something on fixed-bundled q-grids"
     */
    class RadialProjector
    {
        public:
            RadialProjector() {}
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

            /**
             * @brief perform analytical version of the Fourier transform:
             * F(q) = int(f(r)*exp(-iq.r) d^3r)
             *      = 4*pi/sqrt(omega) * i^l * Jl[f](q) * Ylm(q)
             * , where Ylm(q) is real spherical harmonic function, and Jl[f](q) is 
             * the Spherial Bessel Transform of f(r):
             * Jl[f](q) = int(f(r)*j_l(q*r)*r^2 dr)
             * , where j_l(q*r) is the spherical Bessel function of the first kind.
             * 
             */
            
            void sbtft(const std::vector<ModuleBase::Vector3<double>>& qs,
                       std::vector<std::complex<double>>& out,
                       const double& omega,
                       const double& tpiba,
                       const char type = 'r'); // 'r' for ket |>, 'l' for bra <|
            
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