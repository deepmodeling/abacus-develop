#ifndef LIBXC_ABACUS_H
#define LIBXC_ABACUS_H

#ifdef USE_LIBXC

#include "source_base/matrix.h"
#include "source_base/vector3.h"

#include <xc.h>
#include <xc_funcs.h>

#include <tuple>
#include <vector>

#include <map>
#include <utility>

class Charge;

namespace XC_Functional_Libxc
{
//-------------------
//  libxc_setup.cpp
//-------------------

    // sets functional type, which allows combination of LIBXC keyword connected by "+"
    //        for example, "XC_LDA_X+XC_LDA_C_PZ"
    extern std::pair<int, std::vector<int>> set_xc_type_libxc(const std::string& xc_func_in);

    /**
     * @brief instantiate the XC functional by its ID, and set the external parameters if provided.
     * 
     * @param func_id libxc ID of functional, see https://libxc.gitlab.io/functionals/ for details
     * @param xc_polarized 0: unpolarized, 1: spin-polarized
     * @return std::vector<xc_func_type> 
     * 
     * @note the functionality of this method is extended by supporting the user-defined
     *       external parameters of xc. However, there are several functionals' external 
     *       parameters are pre-defined in the code, which herein we call those are 
     *       "in-built" parameters. If the same functional ID is found in both in-built 
     *       and external parameters, the external parameters will overwrite the in-built ones.
     *       The external parameters can be passed here by keywords xc_exch_ext and
     *       xc_corr_ext in the input file. The expected format would be an XC ID 
     *       followed by a list of parameters.
     */
    extern std::vector<xc_func_type> init_func(
        const std::vector<int> &func_id,
        const int xc_polarized,
        const double hybrid_alpha,
        const double hse_omega);

    extern void finish_func(std::vector<xc_func_type> &funcs);


//-------------------
//  libxc_pot.cpp
//-------------------

    extern std::tuple<double, double, ModuleBase::matrix> v_xc_libxc(
        const std::vector<int> &func_id,
        const int &nrxx, // number of real-space grid
        const double &omega, // volume of cell
        const double tpiba,
        const Charge* const chr, // charge density
        const int nspin,
        const bool domag,
        const bool domag_z,
        const std::map<int, double>* scaling_factor,
        const double hybrid_alpha,
        const double hse_omega);

    // for mGGA functional
    extern std::tuple<double, double, ModuleBase::matrix, ModuleBase::matrix> v_xc_meta(
        const std::vector<int> &func_id,
        const int &nrxx, // number of real-space grid
        const double &omega, // volume of cell
        const double tpiba,
        const Charge* const chr,
        const int nspin,
        const double hybrid_alpha,
        const double hse_omega);


//-------------------
//  libxc_tools.cpp
//-------------------

    // converting rho (abacus=>libxc)
    extern std::vector<double> convert_rho(
        const int nspin,
        const std::size_t nrxx,
        const Charge* const chr);

    // converting rho (abacus=>libxc)
    extern std::tuple<std::vector<double>, std::vector<double>> convert_rho_amag_nspin4(
        const int nspin,
        const std::size_t nrxx,
        const Charge* const chr);

    // calculating grho
    extern std::vector<std::vector<ModuleBase::Vector3<double>>> cal_gdr(
        const int nspin,
        const std::size_t nrxx,
        const std::vector<double> &rho,
        const double tpiba,
        const Charge* const chr);

    // converting grho (abacus=>libxc)
    extern std::vector<double> convert_sigma(
        const std::vector<std::vector<ModuleBase::Vector3<double>>> &gdr);

    /// Calculate Laplacian of density using spectral method: ∇²ρ = IFFT(-|G|²·FFT(ρ))
    /// @see cal_lapl() implementation in libxc_tools.cpp for full documentation
    extern std::vector<double> cal_lapl(
        const int nspin,
        const std::size_t nrxx,
        const std::vector<double> &rho,
        const double tpiba,
        const Charge* const chr);

    /// Calculate Laplacian of density using finite-difference kernel (bounded at high G)
    /// @see cal_lapl_fd() implementation in libxc_tools.cpp for full documentation
    extern std::vector<double> cal_lapl_fd(
        const int nspin,
        const std::size_t nrxx,
        const std::vector<double> &rho,
        const Charge* const chr);

    /// Calculate density Hessian: H_ab = ∂²ρ/∂r_a∂r_b (6 independent components per spin)
    /// @see cal_rho_hessian() implementation in libxc_tools.cpp for full documentation
    extern std::vector<std::vector<double>> cal_rho_hessian(
        const int nspin,
        const std::size_t nrxx,
        const std::vector<double> &rho,
        const Charge* const chr);

    extern std::vector<double> cal_sgn(
        const double rho_threshold,
        const double grho_threshold,
        const xc_func_type &func,
        const int nspin,
        const std::size_t nrxx,
        const std::vector<double> &rho,
        const std::vector<double> &sigma);

    // converting etxc from exc (libxc=>abacus)
    extern double convert_etxc(
        const int nspin,
        const std::size_t nrxx,
        const std::vector<double> &sgn,
        const std::vector<double> &rho,
        std::vector<double> exc);

    // converting vtxc and v from vrho and vsigma (libxc=>abacus)
    extern std::pair<double, ModuleBase::matrix> convert_vtxc_v(
        const xc_func_type &func,
        const int nspin,
        const std::size_t nrxx,
        const std::vector<double> &sgn,
        const std::vector<double> &rho,
        const std::vector<std::vector<ModuleBase::Vector3<double>>> &gdr,
        const std::vector<double> &vrho,
        const std::vector<double> &vsigma,
        const double tpiba,
        const Charge* const chr);

    // dh for gga v
    extern std::vector<std::vector<double>> cal_dh(
        const int nspin,
        const std::size_t nrxx,
        const std::vector<double> &sgn,
        const std::vector<std::vector<ModuleBase::Vector3<double>>> &gdr,
        const std::vector<double> &vsigma,
        const double tpiba,
        const Charge* const chr);

    // convert v for NSPIN=4
    extern ModuleBase::matrix convert_v_nspin4(
        const std::size_t nrxx,
        const Charge* const chr,
        const std::vector<double> &amag,
        const ModuleBase::matrix &v);


//-------------------
//  libxc_lda_wrap.cpp
//-------------------

    extern void xc_spin_libxc(
        const std::vector<int> &func_id,
        const double &rhoup,
        const double &rhodw,
        double &exc,
        double &vxcup,
        double &vxcdw,
        const double hybrid_alpha,
        const double hse_omega);


//-------------------
//  libxc_gga_wrap.cpp
//-------------------

    // the entire GGA functional, for nspin=1 case
    extern void gcxc_libxc(
        const std::vector<int> &func_id,
        const double &rho,
        const double &grho,
        double &sxc,
        double &v1xc,
        double &v2xc,
        const double hybrid_alpha,
        const double hse_omega);

    // the entire GGA functional, for nspin=2 case
    extern void gcxc_spin_libxc(
        const std::vector<int> &func_id,
        const double rhoup,
        const double rhodw,
        const ModuleBase::Vector3<double> gdr1,
        const ModuleBase::Vector3<double> gdr2,
        double &sxc,
        double &v1xcup,
        double &v1xcdw,
        double &v2xcup,
        double &v2xcdw,
        double &v2xcud,
        const double hybrid_alpha,
        const double hse_omega);


//-------------------
//  libxc_mgga_wrap.cpp
//-------------------

    /// Wrapper for meta-GGA functionals (single spin channel).
    /// Computes XC energy density and potentials including Laplacian-dependent terms.
    ///
    /// @param func_id  LibXC functional IDs (exchange + correlation)
    /// @param rho      Electron density at grid point
    /// @param grho     |∇ρ|² at grid point
    /// @param lapl_rho ∇²ρ at grid point (density Laplacian; pass 0.0 for SCAN)
    /// @param atau     Kinetic energy density τ at grid point
    /// @param sxc      [out] XC energy density
    /// @param v1xc     [out] ∂ε_xc/∂ρ
    /// @param v2xc     [out] 2·∂ε_xc/∂(|∇ρ|²)
    /// @param v3xc     [out] ∂ε_xc/∂τ
    /// @param vlapl    [out] ∂ε_xc/∂(∇²ρ) (Laplacian potential; 0.0 for SCAN)
    /// @param hybrid_alpha  Exact exchange mixing fraction (0.0 for pure GGA/meta-GGA)
    /// @param hse_omega     Range-separation parameter for HSE-type functionals
    ///
    /// @note For SCAN (non-Laplacian meta-GGA), pass lapl_rho=0.0 and vlapl=0.0.
    ///       For SCANL, pass the actual ∇²ρ value and vlapl will be populated.
    /// @note The vlapl output is used by process_vlapl_potential() to compute the
    ///       FD Laplacian kernel contribution to the XC potential.
    extern void tau_xc(
        const std::vector<int> &func_id,
        const double &rho,
        const double &grho,
        const double &lapl_rho,
        const double &atau,
        double &sxc,
         double &v1xc,
         double &v2xc,
         double &v3xc,
         double &vlapl,
         const double &hybrid_alpha,
         const double &hse_omega);

    /// Wrapper for meta-GGA functionals (spin-polarized, two channels).
    /// Same as tau_xc but handles spin-up/spin-down channels separately.
    ///
    /// @param laplup, lapldw  [in]  ∇²ρ_up, ∇²ρ_dn (Laplacian for each spin channel)
    /// @param vlaplup, vlapldw [out] ∂ε_xc/∂(∇²ρ_up), ∂ε_xc/∂(∇²ρ_dn)
    ///
    /// @note See tau_xc() for other parameter descriptions.
    extern void tau_xc_spin(
        const std::vector<int> &func_id,
        double rhoup,
        double rhodw,
        ModuleBase::Vector3<double> gdr1,
        ModuleBase::Vector3<double> gdr2,
        double laplup,
        double lapldw,
        double tauup,
        double taudw,
        double &sxc,
        double &v1xcup,
        double &v1xcdw,
        double &v2xcup,
        double &v2xcdw,
        double &v2xcud,
         double &v3xcup,
         double &v3xcdw,
         double &vlaplup,
         double &vlapldw,
         const double &hybrid_alpha,
         const double &hse_omega);

} // namespace XC_Functional_Libxc

#endif // USE_LIBXC

#endif // LIBXC_ABACUS_H
