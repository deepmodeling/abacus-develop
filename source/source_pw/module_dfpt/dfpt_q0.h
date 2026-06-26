// ============================================================
// This code is added by Mohan Chen on 2026-05-18.
// This code is currently in the design phase and has not been
// put into production yet. It may change in the future.
// Please use this code with caution. Only developers who know
// what they are doing should use this code.
// ============================================================

#ifndef DFPT_Q0_H
#define DFPT_Q0_H

#include "dfpt_pw_data.h"
#include "source_cell/unitcell.h"
#include "source_psi/psi.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_basis/module_pw/pw_basis_k.h"

namespace ModuleDFPT {

class DFPT_Q0 {
public:
    DFPT_Q0();
    ~DFPT_Q0();
    
    void init(UnitCell& ucell, ModulePW::PW_Basis* pw_rho, 
              ModulePW::PW_Basis_K* pw_wfc);
    
    void compute_eps(const psi::Psi<std::complex<double>>& psi, 
                     const ModuleBase::matrix& wg, DFPT_PW_Data& data);
    
    void compute_born(const psi::Psi<std::complex<double>>& psi, DFPT_PW_Data& data);
    
    void compute_q0_response(DFPT_PW_Data& data);
    
private:
    UnitCell* ucell_ = nullptr;
    ModulePW::PW_Basis* pw_rho_ = nullptr;
    ModulePW::PW_Basis_K* pw_wfc_ = nullptr;
    
    void pos_matrix(const psi::Psi<std::complex<double>>& psi,
                    std::vector<std::vector<ModuleBase::Vector3<std::complex<double>>>>& r_mat);
};

} // namespace ModuleDFPT

#endif // DFPT_Q0_H