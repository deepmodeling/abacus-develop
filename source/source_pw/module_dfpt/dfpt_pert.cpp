// ============================================================
// This code is added by Mohan Chen on 2026-05-18.
// This code is currently in the design phase and has not been
// put into production yet. It may change in the future.
// Please use this code with caution. Only developers who know
// what they are doing should use this code.
// ============================================================

#include "dfpt_pert.h"

namespace ModuleDFPT {

DFPT_Pert::DFPT_Pert() {}

DFPT_Pert::~DFPT_Pert() {}

void DFPT_Pert::init(UnitCell& ucell, ModulePW::PW_Basis* pw_rho, 
                     ModulePW::PW_Basis_K* pw_wfc, Structure_Factor& sf) {
    ucell_ = &ucell;
    pw_rho_ = pw_rho;
    pw_wfc_ = pw_wfc;
    sf_ = &sf;
}

void DFPT_Pert::build_dv(int q_idx, int atom_idx, int dir, DFPT_PW_Data& data) {
    (void)q_idx;
    (void)atom_idx;
    (void)dir;
    (void)data;
}

void DFPT_Pert::apply_dv(int q_idx, int k_idx, const psi::Psi<std::complex<double>>& psi, 
                         DFPT_PW_Data& data) {
    (void)q_idx;
    (void)k_idx;
    (void)psi;
    (void)data;
}

void DFPT_Pert::build_efield(const ModuleBase::Vector3<double>& field, DFPT_PW_Data& data) {
    (void)field;
    (void)data;
}

void DFPT_Pert::dVloc_dtau(int atom_idx, int dir, const ModuleBase::Vector3<double>& q, 
                           std::vector<std::complex<double>>& dv) {
    (void)atom_idx;
    (void)dir;
    (void)q;
    (void)dv;
}

void DFPT_Pert::dVnl_dtau(int atom_idx, int dir, const ModuleBase::Vector3<double>& q,
                          const psi::Psi<std::complex<double>>& psi, int k_idx,
                          std::vector<std::vector<std::complex<double>>>& dv_psi) {
    (void)atom_idx;
    (void)dir;
    (void)q;
    (void)psi;
    (void)k_idx;
    (void)dv_psi;
}

} // namespace ModuleDFPT