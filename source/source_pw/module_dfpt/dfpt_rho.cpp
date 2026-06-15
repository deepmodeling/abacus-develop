// ============================================================
// This code is added by Mohan Chen on 2026-05-18.
// This code is currently in the design phase and has not been
// put into production yet. It may change in the future.
// Please use this code with caution. Only developers who know
// what they are doing should use this code.
// ============================================================

#include "dfpt_rho.h"

namespace ModuleDFPT {

DFPT_Rho::DFPT_Rho() {}

DFPT_Rho::~DFPT_Rho() {
    if (mixer_ != nullptr) {
        delete mixer_;
        mixer_ = nullptr;
    }
}

void DFPT_Rho::init(int nspin, int nrxx, ModulePW::PW_Basis* pw_rho, 
                    ModulePW::PW_Basis_K* pw_wfc, const std::string& mix_type, 
                    double mix_beta) {
    nspin_ = nspin;
    nrxx_ = nrxx;
    pw_rho_ = pw_rho;
    pw_wfc_ = pw_wfc;
    (void)mix_type;
    (void)mix_beta;
}

void DFPT_Rho::compute_drho(const psi::Psi<std::complex<double>>& psi, 
                            const ModuleBase::matrix& wg, int q_idx, 
                            DFPT_PW_Data& data) {
    (void)psi;
    (void)wg;
    (void)q_idx;
    (void)data;
}

void DFPT_Rho::mix_drho(int q_idx, DFPT_PW_Data& data) {
    (void)q_idx;
    (void)data;
}

double DFPT_Rho::get_residual(int q_idx, DFPT_PW_Data& data) const {
    (void)q_idx;
    (void)data;
    return 0.0;
}

} // namespace ModuleDFPT