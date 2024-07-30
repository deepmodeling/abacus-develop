#ifndef RADIAL_PROJECTION_H
#define RADIAL_PROJECTION_H

#include "module_base/vector3.h"
#include "module_base/matrix.h"
#include "module_base/math_ylmreal.h"
#include "module_base/spherical_bessel_transformer.h"
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
    class RadialProjector
    {
        public:
            RadialProjector(const int lmax, 
                            const std::vector<ModuleBase::Vector3<double>>& qs,
                            const double& omega = 1.0,
                            const double& tpiba = 1.0)
            {
                // can cache the Ylm and SphericalBesselTransformer at the very beginning
                const int total_lm = std::pow(lmax+1, 2);
                const int npw = qs.size();
                ylm_.create(total_lm, npw);
                ModuleBase::YlmReal::Ylm_Real(total_lm, npw, qs.data(), ylm_);
                sbt_ = std::unique_ptr<ModuleBase::SphericalBesselTransformer>(
                    new ModuleBase::SphericalBesselTransformer(true));
                qnorm_.resize(npw);
                std::transform(qs.begin(), qs.end(), qnorm_.begin(), 
                [tpiba](const ModuleBase::Vector3<double>& q){return tpiba * q.norm();});
                omega_ = omega;
            }
            ~RadialProjector() {}
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
            void sbtfft(const int nr,
                        const double* r,
                        const double* in,
                        const int l,
                        std::complex<double>* out) const;
            void sbtfft(const std::vector<double>& r,
                        const std::vector<double>& in,
                        const int l,
                        std::vector<std::complex<double>>& out) const;

        private:
            ModuleBase::matrix ylm_;
            std::unique_ptr<ModuleBase::SphericalBesselTransformer> sbt_;
            std::vector<double> qnorm_;
            double omega_;
    };

    /**
     * @brief build a map from [it][l][zeta] to vectorized index, together with its reverse
     * 
     * @param ntype number of atom types
     * @param lmax maximal angular momentum for each kind
     * @param nzeta number of zeta for each l for each kind
     * @param map [out] map from [it][l][zeta] to vectorized index
     * @param rmap [out] reverse map from vectorized index to [it][l][zeta]
     */
    void _radial_indexing(const int ntype,                                  //< from UnitCell::ntype
                          const std::vector<int>& lmax,
                          const std::vector<std::vector<int>>& nzeta,
                          std::map<std::tuple<int, int, int>, int>& map,
                          std::vector<std::tuple<int, int, int>>& rmap);
    
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
     * @brief initialize a q-grid for SBFFT
     * 
     * @param ecutwfc 
     * @param ecutrho 
     * @param R 
     * @return ModulePW::PW_Basis 
     */
    // static ModulePW::PW_Basis _init_qgrid(const double& ecutwfc,
    //                                       const double& ecutrho,
    //                                       const ModuleBase::Matrix3& R); // do I need to change to std::vector?

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