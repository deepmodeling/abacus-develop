#ifndef RADIAL_PROJECTION_H
#define RADIAL_PROJECTION_H
// project any atom-centered function that has seperatable radial and angular parts
// onto the planewave basis
#include <vector>
#include <complex>
#include <map>
#include <utility>
#include "module_base/vector3.h"
#include "module_basis/module_pw/pw_basis.h"

namespace RadialProjection
{
    /**
     * ====================================================================================
     * 
     *              Basic functions for performing analytical Fourier transform
     * 
     * ====================================================================================
     * 
     */
    /**
     * @brief perform analytical version of the Fourier transform:
     * F(q) = int(f(r)*exp(-iq.r) d^3r)
     *      = 4*pi/sqrt(omega) * i^l * Jl[f](q) * Ylm(q)
     * , where Ylm(q) is real spherical harmonic function, and Jl[f](q) is 
     * the Spherial Bessel Transform of f(r):
     * Jl[f](q) = int(f(r)*j_l(q*r)*r^2 dr)
     * , where j_l(q*r) is the spherical Bessel function of the first kind.
     * 
     * @note
     * One may need the opposite phase, in such case, simply do out.conj().
     * This function also has overloaded version for std::vector input, and
     * cached version for ylm_q.
     * 
     * @param nr number of real space grid
     * @param r real space grid (radial)
     * @param in radial part of the input function
     * @param l angular momentum of the angular part
     * @param m magnetic quantum number of the angular part
     * @param q wavevector
     * @param omega cell volume
     * @param tpiba 2*pi/lat0
     * @param out output value
     */
    static void _sbtfft(const int nr,
                        const double* r,
                        const double* in,
                        const int l,
                        const int m,
                        const ModuleBase::Vector3<double>& q,
                        const double& omega,
                        const double& tpiba,
                        std::complex<double>& out);
    // for cached call of ylm_q. Details on param list, see the above function
    static void _sbtfft(const int nr,
                        const double* r,
                        const double* in,
                        const int l,
                        const double& qnorm,
                        const double& ylm_q,
                        const double& omega,
                        std::complex<double>& out);
    // std::vector version
    static void _sbtfft(const std::vector<double>& r,
                        const std::vector<double>& in,
                        const int l,
                        const int m,
                        const ModuleBase::Vector3<double>& q,
                        const double& omega,
                        const double& tpiba,
                        std::complex<double>& out);
    // for cached call of ylm_q, std::vector version
    static void _sbtfft(const std::vector<double>& r,
                        const std::vector<double>& in,
                        const int l,
                        const double& qnorm,
                        const double& ylm_q,
                        const double& omega,
                        std::complex<double>& out);
    /**
     * @brief project any atom-centered function that has seperatable radial and angular parts
     * onto the planewave basis (based on Spherical Bessel Transformer)
     * 
     * @param nr number of real space grid
     * @param r real space grid
     * @param in input function
     * @param l angular momentum
     * @param qs wavevectors
     * @param omega cell volume
     * @param tpiba 2*pi/lat0
     * @param out output value
     */
    static void _sbtfft_for_each_q(const int nr,
                                   const double* r,
                                   const double* in,
                                   const int l,
                                   const std::vector<ModuleBase::Vector3<double>>& qs, // do I need to change to std::vector?
                                   const double& omega,
                                   const double& tpiba,
                                   std::complex<double>* out);
    // std::vector version
    static void _sbtfft_for_each_q(const std::vector<double>& r,
                                   const std::vector<double>& in,
                                   const int l,
                                   const std::vector<ModuleBase::Vector3<double>>& qs, // do I need to change to std::vector?
                                   const double& omega,
                                   const double& tpiba,
                                   std::vector<std::complex<double>>& out);
    /**
     * @brief project a collection of atom-centered functions that have seperatable radial and angular parts
     * onto the planewave basis (based on Spherical Bessel Transformer). In principle the collection
     * is organized as flattened [it][l][zeta], and l is organized as [i], the index of the input function
     * 
     * @param r uniformed radial grid
     * @param in_clkxn collection(clkxn) of input functions
     * @param l_clkxn collection of angular momentum for each input function
     * @param qs wavevectors
     * @param omega cell volume
     * @param tpiba 2*pi/lat0
     * @param out output value
     */
    static void _proj_each_type(const std::vector<double>& r,
                                const std::vector<std::vector<double>>& in_clkxn,
                                const std::vector<int>& l_clkxn,
                                const std::vector<ModuleBase::Vector3<double>>& qs, //< from PW_Basis::getgpluskcar
                                const double& omega,                                //< from UnitCell::omega
                                const double& tpiba,                                //< from UnitCell::tpiba
                                std::vector<std::complex<double>>& out);
    
    /**
     * @brief raw pointer version of project function. This version can support in-continuous
     * memory layout of input functions. For example if there are beta distributed in Atom
     * instances, build a vector of pointers, [it][l][zeta]. Then save l into another raw pointer
     * , indexed by [izeta]
     * 
     * @param nr number of real space grid
     * @param r C-style array of real space grid
     * @param in_clkxn C-style array of input functions, indexed by [iproj][ir]
     * @param l_clkxn C-style array of angular momentum for each input function, indexed by [iproj]
     * @param npw number of wavevectors
     * @param qs C-style array of wavevectors
     * @param omega cell volume
     * @param tpiba 2*pi/lat0
     * @param out C-style array of output values, indexed by [iproj][iq]
     */
    static void _proj_each_type(const int nr,
                                const double* r,
                                const std::vector<double*> in_clkxn,
                                const int* l_clkxn,
                                const std::vector<ModuleBase::Vector3<double>>& qs, // do I need to change to std::vector?
                                const double omega,
                                const double tpiba,
                                std::complex<double>* out);

    static void proj_all_atoms(const int nr,
                               const double* r,
                               const std::vector<std::vector<double*>>& in_clkxn,
                               const std::vector<int>& l_clkxn,
                               const std::vector<ModuleBase::Vector3<double>>& qs,
                               const double omega,
                               const double tpiba,
                               std::vector<std::complex<double>>& out);

    /**
     * @brief build a map from [it][l][zeta] to vectorized index, together with its reverse
     * 
     * @param ntype number of atom types
     * @param lmax maximal angular momentum for each kind
     * @param nzeta number of zeta for each l for each kind
     * @param map [out] map from [it][l][zeta] to vectorized index
     * @param rmap [out] reverse map from vectorized index to [it][l][zeta]
     */
    static void _indexing(const int ntype,                                  //< from UnitCell::ntype
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
    static ModulePW::PW_Basis _init_qgrid(const double& ecutwfc,
                                          const double& ecutrho,
                                          const ModuleBase::Matrix3& R); // do I need to change to std::vector?

    /**
     * @brief get the mask function for SBFFT
     * 
     * @param mask mask function
     */
    static void _mask_func(std::vector<double>& mask);

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
    static void _do_mask_on_radial(const int nr1,
                                   const double* r,
                                   const double* in,
                                   const int nr2,
                                   const double* mask,
                                   double* out);
}

#endif // RADIAL_PROJECTION_H