#ifndef DFPT_XC_ANALYTIC_H
#define DFPT_XC_ANALYTIC_H

#ifdef __LIBXC

#include "dfpt_rho.h"

#include "source_base/matrix3.h"
#include "source_basis/module_pw/pw_basis.h"

#include <vector>

namespace ModuleDFPT
{

/**
 * @brief Analytic first-order XC kernel built from libxc second derivatives
 *        (C6 provider; replaces the finite-difference adapter when ABACUS
 *        is compiled with libxc).
 *
 * The response to the q-shifted density amplitude A(r) of the physical
 * first-order density drho(r) = 2 Re[A(r) e^{i q r}] is evaluated in the
 * dgradcorr form (QE LR_Modules/dgradcorr.f90): the gradient operator acts
 * on the total physical field through i G -> i (G + q), so the carrier
 * contribution i (q . f_grad) A -- exact for any q, and structurally
 * missing from grid-only finite differences of the carrier-free amplitude
 * -- is included. LDA functionals reduce to dvxc = f_rho A (their gradient
 * kernel vanishes identically, no q dependence).
 *
 * The libxc kernel fields are evaluated once at construction on the
 * ground-state density rho0 and stored in Ry units (libxc Hartree values
 * scaled by e2 = 2). Real-space points with |rho0| <= epsr or
 * |grad rho0|^2 <= epsg follow the QE dgradcorr guard: the gradient-kernel
 * fields vanish there. Construction is the only libxc entry point
 * (handles are opened and closed inside the constructor), so this header
 * does not include xc.h.
 */
class XC_First_Order_Analytic : public XC_First_Order
{
  public:
    XC_First_Order_Analytic(const std::vector<int>& func_id,
                            const int func_type,
                            const std::vector<double>& rho0_r,
                            ModulePW::PW_Basis* pw_rho);

    void set_context(const ModuleBase::Vector3<double>& q_cart) const override;

    void apply(const std::vector<std::complex<double>>& drho_r,
               std::vector<std::complex<double>>& dvxc_r) const override;

  private:
    /// grho_(ipol)[ir] = d rho0 / d r_ipol on the rho grid (real part)
    void init_gradient_(const std::vector<double>& rho0_r);

    /// kernel fields k_rr_/k_sr_/k_ss_/k_s_ (LDA fills k_rr_ only)
    void init_kernel_fields_(const std::vector<double>& rho0_r);

    int func_type_ = 1; ///< 1 = LDA, 2 = GGA (XC_Functional convention)
    std::vector<int> func_id_; ///< libxc ids, exchange first (XC_Functional order)
    ModulePW::PW_Basis* pw_rho_ = nullptr;

    mutable ModuleBase::Vector3<double> q_cart_; ///< phonon q, ucell::G units (1/lat0)
    mutable bool q_set_ = false; ///< false: apply() behaves as q = 0 (Gamma)

    std::vector<std::vector<double>> grho_; ///< [3][nrxx] ground-state gradient
    std::vector<double> k_rr_; ///< rho-rho kernel field, Ry units
    std::vector<double> k_sr_; ///< rho-sigma field (s1 = grho.gdrho coefficient)
    std::vector<double> k_ss_; ///< sigma-sigma field (s1^2 coefficient)
    std::vector<double> k_s_; ///< sigma field (gdrho coefficient)
};

} // namespace ModuleDFPT

#endif // __LIBXC
#endif // DFPT_XC_ANALYTIC_H
