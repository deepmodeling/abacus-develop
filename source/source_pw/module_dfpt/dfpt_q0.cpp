// ============================================================
// This code is added by Mohan Chen on 2026-05-18.
// This code is currently in the design phase and has not been
// put into production yet. It may change in the future.
// Please use this code with caution. Only developers who know
// what they are doing should use this code.
// ============================================================

#include "dfpt_q0.h"

namespace ModuleDFPT {

DFPT_Q0::DFPT_Q0() {}

DFPT_Q0::~DFPT_Q0() {}

void DFPT_Q0::init(UnitCell& ucell, ModulePW::PW_Basis* pw_rho, 
                   ModulePW::PW_Basis_K* pw_wfc) {
    ucell_ = &ucell;
    pw_rho_ = pw_rho;
    pw_wfc_ = pw_wfc;
}

void DFPT_Q0::compute_eps(const psi::Psi<std::complex<double>>& psi, 
                          const ModuleBase::matrix& wg, DFPT_PW_Data& data) {
    (void)psi;
    (void)wg;
    (void)data;
}

void DFPT_Q0::compute_born(const psi::Psi<std::complex<double>>& psi, DFPT_PW_Data& data) {
    (void)psi;
    (void)data;
}

void DFPT_Q0::compute_q0_response(DFPT_PW_Data& data) {
    (void)data;
}

void DFPT_Q0::pos_matrix(const psi::Psi<std::complex<double>>& psi,
                         std::vector<std::vector<ModuleBase::Vector3<std::complex<double>>>>& r_mat) {
    (void)psi;
    (void)r_mat;
}

} // namespace ModuleDFPT