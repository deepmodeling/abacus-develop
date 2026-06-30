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
#include "source_psi/psi.h"
#include "source_estate/module_charge/charge_mixing.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_basis/module_pw/pw_basis_k.h"

namespace ModuleDFPT {

class DFPT_Rho {
public:
    DFPT_Rho();
    ~DFPT_Rho();
    
    void init(int nspin, int nrxx, ModulePW::PW_Basis* pw_rho, 
              ModulePW::PW_Basis_K* pw_wfc, const std::string& mix_type, 
              double mix_beta);
    
    void compute_drho(const psi::Psi<std::complex<double>>& psi, 
                      const ModuleBase::matrix& wg, int q_idx, 
                      DFPT_PW_Data& data);
    
    void mix_drho(int q_idx, DFPT_PW_Data& data);
    
    double get_residual(int q_idx, DFPT_PW_Data& data) const;

private:
    int nspin_ = 1;
    int nrxx_ = 0;
    ModulePW::PW_Basis* pw_rho_ = nullptr;
    ModulePW::PW_Basis_K* pw_wfc_ = nullptr;
    
    Charge_Mixing* mixer_ = nullptr;
    
    std::vector<std::vector<std::vector<double>>> drho_in_;
    std::vector<std::vector<std::vector<double>>> drho_out_;
};

} // namespace ModuleDFPT

#endif // DFPT_RHO_H
