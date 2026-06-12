// ============================================================
// This code is added by Mohan Chen on 2026-05-18.
// This code is currently in the design phase and has not been
// put into production yet. It may change in the future.
// Please use this code with caution. Only developers who know
// what they are doing should use this code.
// ============================================================

#include "dfpt_metal.h"

namespace ModuleDFPT {

DFPT_Metal::DFPT_Metal() {}

DFPT_Metal::~DFPT_Metal() {}

void DFPT_Metal::init(double sigma, const std::string& smearing_type) {
    sigma_ = sigma;
    smearing_type_ = smearing_type;
}

void DFPT_Metal::dfdeps(const ModuleBase::matrix& eig, double efermi, 
                         ModuleBase::matrix& dfdeps) {
    (void)eig;
    (void)efermi;
    (void)dfdeps;
}

void DFPT_Metal::compute_dmu(int q_idx, const psi::Psi<std::complex<double>>& psi,
                             const ModuleBase::matrix& wg, const ModuleBase::matrix& dfdeps,
                             DFPT_PW_Data& data) {
    (void)q_idx;
    (void)psi;
    (void)wg;
    (void)dfdeps;
    (void)data;
}

void DFPT_Metal::compute_drho_metal(int q_idx, const psi::Psi<std::complex<double>>& psi,
                                    const ModuleBase::matrix& wg, const ModuleBase::matrix& dfdeps,
                                    double dmu, DFPT_PW_Data& data) {
    (void)q_idx;
    (void)psi;
    (void)wg;
    (void)dfdeps;
    (void)dmu;
    (void)data;
}

double DFPT_Metal::fd_dfdeps(double e, double efermi) {
    (void)e;
    (void)efermi;
    return 0.0;
}

double DFPT_Metal::gauss_dfdeps(double e, double efermi) {
    (void)e;
    (void)efermi;
    return 0.0;
}

} // namespace ModuleDFPT