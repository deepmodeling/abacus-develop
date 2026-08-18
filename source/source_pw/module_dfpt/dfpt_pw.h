// ============================================================
// This code is added by Mohan Chen on 2026-05-18.
// This code is currently in design phase and has not been
// put into production yet. It may change in the future.
// Please use this code with caution. Only developers who know
// what they are doing should use this code.
// ============================================================

#ifndef DFPT_PW_H
#define DFPT_PW_H

#include "source_base/matrix.h"
#include "source_cell/unitcell.h"
#include "source_psi/psi.h"

#include <string>
#include <vector>

class Plus_U;
class Structure_Factor;

namespace ModulePW {
class PW_Basis;
class PW_Basis_K;
}

namespace ModuleDFPT {

class XC_First_Order;

/**
 * @brief Density-functional perturbation theory driver (plane waves).
 *
 * C7 wiring: init receives the converged ground state (psi, wg, eig), the
 * shared-grid plane-wave bases, the real-space effective potential and the
 * first-order XC kernel contract; run() then drives, per irreducible q, the
 * per-displacement self-consistent Sternheimer cycle
 *   build_dv (bare external)
 *   -> [ dv_sc = v_hartree_q(drho_in) + xc_->apply(drho_in)
 *        -> Sternheimer solve of (H(k+q) - eps_n) P_c dpsi = -P_c (dV_ext
 *           + dV_sc)|psi_n> with the k+q occupied states as the projector
 *        -> compute_drho -> mix_drho ]*
 *   -> accumulate_electron -> assemble -> diagonalize (+ LO-TO at q = 0).
 * With null bases (design-phase skeleton) run() keeps the documented
 * first-iteration-converged fallback of the irrep bookkeeping loop.
 */
class DFPT_PW {
public:
    DFPT_PW();
    ~DFPT_PW();

    void init(UnitCell& ucell, const psi::Psi<std::complex<double>>& psi,
              ModulePW::PW_Basis* pw_rho, ModulePW::PW_Basis_K* pw_wfc,
              Structure_Factor* sf, const std::vector<double>& veff_r,
              const ModuleBase::matrix& wg, const ModuleBase::matrix& eig,
              const XC_First_Order* xc,
              double nelec, double ecutwfc, const Plus_U* dftu);

    void run();

    /// DFT+U reservation accessors (U0): with_u() reports whether a DFT+U
    /// provider is wired (dft_plus_u enabled upstream); u_active() further
    /// requires the provider to be usable (locale initialized, i.e. the LCAO
    /// orbital files are present).
    bool get_with_u() const;
    bool get_u_active() const;

    std::vector<double> get_phonon_freq(int q_idx) const;

    ModuleBase::matrix get_dielectric_tensor() const;

    ModuleBase::matrix get_born_charges(int atom_idx) const;

    /// q-point source: a q list file overrides the Monkhorst-Pack q mesh
    void set_qfile(const std::string& filename);

    void set_qmesh(int nqx, int nqy, int nqz);

    void set_conv_thr(double thr);

    void set_max_iter(int max_iter);

    void set_mix_beta(double beta);

    /// q = 0 response switches (epsilon_inf / Born charges / LO-TO)
    void set_compute_q0(bool flag);

    void set_loto(bool flag);

private:
    class Impl;
    Impl* pimpl_;
};

} // namespace ModuleDFPT

#endif // DFPT_PW_H
