// ============================================================
// This code is added by Mohan Chen on 2026-05-18.
// This code is currently in the design phase and has not been
// put into production yet. It may change in the future.
// Please use this code with caution. Only developers who know
// what they are doing should use this code.
// ============================================================

#ifndef DFPT_RHO_H
#define DFPT_RHO_H

#include "dfpt_pw_data.h"
#include "source_base/matrix3.h"
#include "source_psi/psi.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_basis/module_pw/pw_basis_k.h"
#include <string>
#include <vector>

namespace Base_Mixing
{
class Plain_Mixing;
}

namespace ModuleDFPT {

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
 * mix_drho applies plain mixing on the q-shifted coefficients:
 *   drho_in <- drho_in + beta (drho_out - drho_in)
 * through Base_Mixing::Plain_Mixing (no Charge_Mixing / Charge dependency).
 */
class DFPT_Rho {
public:
    DFPT_Rho();
    ~DFPT_Rho();
    
    void init(int nspin, int nrxx, ModulePW::PW_Basis* pw_rho,
              ModulePW::PW_Basis_K* pw_wfc,
              const ModuleBase::Matrix3& recip_matrix,
              const std::string& mix_type, double mix_beta);
    
    void compute_drho(const psi::Psi<std::complex<double>>& psi, 
                      const ModuleBase::matrix& wg, int q_idx, 
                      DFPT_PW_Data& data);
    
    /// first-order occupation matrix (docc) for DFT+U (U0 reservation).
    void cal_docc(const psi::Psi<std::complex<double>>& psi, 
                  const ModuleBase::matrix& wg, int q_idx, 
                  DFPT_PW_Data& data);
    
    void mix_drho(int q_idx, DFPT_PW_Data& data);
    
    double get_residual(int q_idx, DFPT_PW_Data& data) const;

private:
    int nspin_ = 1;
    int nrxx_ = 0;
    ModulePW::PW_Basis* pw_rho_ = nullptr;
    ModulePW::PW_Basis_K* pw_wfc_ = nullptr;
    ///< reciprocal lattice matrix in 1/lat0 (UnitCell::G convention)
    ModuleBase::Matrix3 recip_matrix_;
    double mix_beta_ = 0.7;
    
    Base_Mixing::Plain_Mixing* mixer_ = nullptr;
    
    /// mixing state, q-shifted coefficients on the rho grid, [q][spin]
    std::vector<std::vector<std::vector<std::complex<double>>>> drho_in_;
    std::vector<std::vector<std::vector<std::complex<double>>>> drho_out_;
    std::vector<double> residual_;
};

} // namespace ModuleDFPT

#endif // DFPT_RHO_H
