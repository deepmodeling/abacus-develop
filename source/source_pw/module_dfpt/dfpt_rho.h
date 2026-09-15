#ifndef DFPT_RHO_H
#define DFPT_RHO_H

#include "dfpt_pw_data.h"
#include "source_base/matrix3.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_basis/module_pw/pw_basis_k.h"
#include "source_psi/psi.h"

#include <memory>
#include <string>
#include <vector>

namespace Base_Mixing
{
class Plain_Mixing;
}

namespace ModuleDFPT
{

class DFPT_KQ_Basis;

/**
 * @brief First-order exchange-correlation kernel contract (C6).
 *
 * Production adapters live at the esolver wiring layer (C7): the complex
 * q-shifted density amplitude drho_r is split into Re/Im parts, fed
 * through the real-space finite-difference kernel (elecstate::PotXC_FDM,
 * delta V_xc = V_xc[rho0 + drho] - V_xc[rho0]) and recombined - linear
 * superposition is exact up to O(|drho|^2). module_dfpt itself never
 * includes pot_xc_fdm.h (minimal header dependencies), mirroring the
 * DFPT_Stern::LinearOperator injection convention.
 */
class XC_First_Order
{
  public:
    virtual ~XC_First_Order() = default;

    /// dvxc_r(r) = delta V_xc[drho_r](r), complex q-shifted amplitude on
    /// the shared real-space grid. Implementations must not resize or
    /// alias drho_r; dvxc_r is resized to drho_r.size() and fully
    /// overwritten.
    virtual void apply(const std::vector<std::complex<double>>& drho_r, std::vector<std::complex<double>>& dvxc_r) const
        = 0;
};

/**
 * @brief First-order density response (C3).
 *
 * compute_drho builds the q-shifted response density
 *   drho(r) = 2 Re[ e^{i q r} A(r) ],
 *   A(r) = sum_{k,n occ} wg(k,n) u*_nk(r) du_nk(r),
 * where u/du are the periodic parts of psi_nk (k basis) and dpsi_nk (k+q
 * basis); PW_Basis_K / PW_Basis transforms are phase-free (they return the
 * periodic part), so the Bloch phases combine into the single e^{i q r}
 * factor. In reciprocal space drho_g holds the q-shifted coefficients
 * drho_Delta (coefficient of e^{i (Delta+q) r}, indexed by the rho-grid ig),
 * A_Delta = sum_{kn} wg sum_G c*_G(k,n) d_{G+Delta}(k,n); the Delta = -q
 * harmonic is dropped when -q falls on a reciprocal-lattice vector (charge
 * conservation, notably at q = Gamma).
 *
 * mix_drho mixes the q-shifted coefficients:
 *   drho_in <- drho_in + beta_g (drho_out - drho_in)
 * through Base_Mixing::Plain_Mixing (no Charge_Mixing / Charge dependency).
 * With mix_type = "plain" beta_g = beta on every shell; with mix_type =
 * "kerker" beta_g = beta * f_g, the Kerker screen
 *   f_g = |G+q|^2 / (|G+q|^2 + a^2),
 * evaluated with the same q_shifted |G+q| = |gcar + q_cart * recip| (in
 * 1/lat0^2 units) convention as v_hartree_q, so the Coulomb-stiffness
 * eigenvalues concentrated on the smallest shells (lambda ~ -2.2 on
 * {111}/{200} for the diamond smoke case, where plain mixing needs
 * beta < 2 / (1 + |lambda|)) become stabilizable at beta up to 1. The
 * screen is applied to both drho_in and drho_out before the plain mix and
 * the screened-out part is added back, i.e. the stored mixed density is
 * rin + beta_g (out - in) (physical, not screen-scaled); the |G+q| = 0
 * harmonic (f = 0) is frozen, consistent with its drop in compute_drho.
 */
class DFPT_Rho
{
  public:
    /// aggregate config for init (no defaults: every field must be set)
    struct Config
    {
        int nspin;
        int nrxx;
        ModulePW::PW_Basis* pw_rho;
        ModulePW::PW_Basis_K* pw_wfc;
        ModuleBase::Matrix3 recip_matrix;
        std::string mix_type; // "plain" or "kerker"
        double mix_beta;
        double kerker_a2; // Kerker screen a^2, 1/lat0^2 (kerker only)
    };

    DFPT_Rho();
    ~DFPT_Rho();

    void init(const Config& cfg);

    void compute_drho(const psi::Psi<std::complex<double>>& psi,
                      const ModuleBase::matrix& wg,
                      int q_idx,
                      DFPT_PW_Data& data);

    /// first-order occupation matrix (docc) for DFT+U (U0 reservation).
    void cal_docc(const psi::Psi<std::complex<double>>& psi,
                  const ModuleBase::matrix& wg,
                  int q_idx,
                  DFPT_PW_Data& data);

    void mix_drho(int q_idx, DFPT_PW_Data& data);

    /// C7: drop the mixing state of q_idx so the next perturbation at the
    /// same q restarts from a zero input density (the drho_in slot is
    /// indexed by q only, while every (atom, direction) needs its own
    /// self-consistent cycle).
    void reset_mixing(int q_idx);

    /// C6: q-shifted first-order Hartree potential in reciprocal space,
    ///   dV_H(G) = 4 pi e^2 / |G+q|^2 * drho_g,
    /// with the convention aligned with elecstate::H_Hartree_pw::v_hartree
    /// (fac = e2 * FOUR_PI / (tpiba2 * |G+q|^2)); the |G+q| = 0 component
    /// (ig = -q) is skipped. Serves both the C6 q->0 response and the C7
    /// screened potential at every q point.
    void v_hartree_q(const ModuleBase::Vector3<double>& q_cart,
                     const std::vector<std::complex<double>>& drho_g,
                     std::vector<std::complex<double>>& dv_ha_g) const;

    double get_residual(int q_idx, DFPT_PW_Data& data) const;

  private:
    /// accumulate one (ik, ib) band contribution to the real-space amplitude
    void add_band_(int ik,
                   int ib,
                   double w,
                   const std::complex<double>* c_ptr,
                   const std::vector<std::complex<double>>& dpsi,
                   const DFPT_KQ_Basis& kq,
                   std::vector<std::complex<double>>& a_r);

    /// rebuild the real-space manifest drho_r from G-space coefficients
    void make_drho_r_(const std::vector<std::complex<double>>& drho_g,
                      const ModuleBase::Vector3<double>& q_frac,
                      std::vector<double>& drho_r) const;

    /// charge conservation: zero the Delta = -q harmonic when -q is a G vector
    void zero_neg_q_(const ModuleBase::Vector3<double>& q_cart,
                     std::vector<std::complex<double>>& drho_g) const;

    /// true if q is Gamma within 1e-10
    static bool is_gamma_q_(const ModuleBase::Vector3<double>& q);

    int nspin_ = 1;
    int nrxx_ = 0;
    ModulePW::PW_Basis* pw_rho_ = nullptr;
    ModulePW::PW_Basis_K* pw_wfc_ = nullptr;
    ///< reciprocal lattice matrix in 1/lat0 (UnitCell::G convention)
    ModuleBase::Matrix3 recip_matrix_;
    double mix_beta_ = 0.7;
    ///< mixing algorithm: "plain" (beta only) or "kerker" (Kerker screen)
    std::string mix_type_;
    ///< Kerker screening parameter a^2 in 1/lat0^2 (same units as |G+q|^2)
    double kerker_a2_ = 0.0;

    std::unique_ptr<Base_Mixing::Plain_Mixing> mixer_;

    /// mixing state, q-shifted coefficients on the rho grid, [q][spin]
    std::vector<std::vector<std::vector<std::complex<double>>>> drho_in_;
    std::vector<std::vector<std::vector<std::complex<double>>>> drho_out_;
    std::vector<double> residual_;
};

} // namespace ModuleDFPT

#endif // DFPT_RHO_H
