// ============================================================
// This code is added by Mohan Chen on 2026-05-18.
// This code is currently in the design phase and has not been
// put into production yet. It may change in the future.
// Please use this code with caution. Only developers who know
// what they are doing should use this code.
// ============================================================

#ifndef DFPT_PERT_H
#define DFPT_PERT_H

#include "dfpt_pw_data.h"
#include "source_cell/unitcell.h"
#include "source_psi/psi.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_basis/module_pw/pw_basis_k.h"
#include "source_pw/module_pwdft/structure_factor.h"

namespace ModuleDFPT {

class DFPT_Pert {
public:
    DFPT_Pert();
    ~DFPT_Pert();
    
    void init(UnitCell& ucell, ModulePW::PW_Basis* pw_rho, 
              ModulePW::PW_Basis_K* pw_wfc, Structure_Factor& sf);
    
    void build_dv(int q_idx, int atom_idx, int dir, DFPT_PW_Data& data);
    
    void apply_dv(int q_idx, int k_idx, const psi::Psi<std::complex<double>>& psi, 
                  DFPT_PW_Data& data);
    
    void build_efield(const ModuleBase::Vector3<double>& field, DFPT_PW_Data& data);

private:
    UnitCell* ucell_ = nullptr;
    ModulePW::PW_Basis* pw_rho_ = nullptr;
    ModulePW::PW_Basis_K* pw_wfc_ = nullptr;
    Structure_Factor* sf_ = nullptr;
    
    void dVloc_dtau(int atom_idx, int dir, const ModuleBase::Vector3<double>& q, 
                    std::vector<std::complex<double>>& dv);
    
    void dVnl_dtau(int atom_idx, int dir, const ModuleBase::Vector3<double>& q,
                   const psi::Psi<std::complex<double>>& psi, int k_idx,
                   std::vector<std::vector<std::complex<double>>>& dv_psi);
};

} // namespace ModuleDFPT

#endif // DFPT_PERT_H
