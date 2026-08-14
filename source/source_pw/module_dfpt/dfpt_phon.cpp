// ============================================================
// This code is added by Mohan Chen on 2026-05-18.
// This code is currently in the design phase and has not been
// put into production yet. It may change in the future.
// Please use this code with caution. Only developers who know
// what they are doing should use this code.
// ============================================================

#include "dfpt_phon.h"

namespace ModuleDFPT {

DFPT_Phon::DFPT_Phon() {}

DFPT_Phon::~DFPT_Phon() {}

void DFPT_Phon::init(UnitCell& ucell) {
    ucell_ = &ucell;
}

void DFPT_Phon::assemble(int q_idx, DFPT_PW_Data& data) {
    int nat = ucell_->nat;
    ModuleBase::matrix dynmat(3 * nat, 3 * nat);
    dynmat.zero_out();
    
    ModuleBase::Vector3<double> q = data.get_qvec(q_idx);
    ion_ion(q, dynmat);
    electron(q_idx, data, dynmat);
    // DFT+U dynamical-matrix term (U0 reservation, implemented in C5)
    if (data.with_u()) {
        dftu_onsite(q_idx, data);
    }
    
    data.set_dynmat(q_idx, dynmat);
}

void DFPT_Phon::diagonalize(int q_idx, DFPT_PW_Data& data) {
    ModuleBase::matrix dynmat = data.get_dynmat(q_idx);
    int nat = ucell_->nat;
    
    std::vector<double> freq(3 * nat, 0.0);
    for (int i = 0; i < 3 * nat; ++i) {
        freq[i] = static_cast<double>(i);
    }
    
    data.set_phon_freq(q_idx, freq);
}

void DFPT_Phon::add_loto(DFPT_PW_Data& data) {
    (void)data;
}

bool DFPT_Phon::check_sum_rule(int q_idx, DFPT_PW_Data& data) const {
    (void)q_idx;
    (void)data;
    return true;
}

void DFPT_Phon::ion_ion(const ModuleBase::Vector3<double>& q, ModuleBase::matrix& dyn) {
    (void)q;
    (void)dyn;
}

void DFPT_Phon::electron(int q_idx, DFPT_PW_Data& data, ModuleBase::matrix& dyn) {
    (void)q_idx;
    (void)data;
    (void)dyn;
}

void DFPT_Phon::dftu_onsite(int q_idx, DFPT_PW_Data& data) {
    // Reserved DFT+U contribution to the dynamical matrix (U0).
    // The physical implementation lands in C5 (dftu_lambda electron term):
    //   sum_nk w_nk [ <psi|dV_U|dpsi> + frozen second-order U term
    //   (~ becp * V_eff * dbecp_f contractions) ], accumulated into the
    //   dynamical matrix. dV_U itself is assembled by DFPT_Pert::build_dv_u.
    (void)q_idx;
    (void)data;
}

void DFPT_Phon::ewald_sum(const ModuleBase::Vector3<double>& q, ModuleBase::matrix& dyn) {
    (void)q;
    (void)dyn;
}

} // namespace ModuleDFPT