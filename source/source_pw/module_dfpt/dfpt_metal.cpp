// ============================================================
// This code is added by Mohan Chen on 2026-05-18.
// This code is currently in the design phase and has not been
// put into production yet. It may change in the future.
// Please use this code with caution. Only developers who know
// what they are doing should use this code.
// ============================================================

#include "dfpt_metal.h"
#include "source_base/tool_quit.h"

namespace ModuleDFPT {

DFPT_Metal::DFPT_Metal() {}

DFPT_Metal::~DFPT_Metal() {}

void DFPT_Metal::init(double sigma, const std::string& smearing_type) {
    sigma_ = sigma;
    smearing_type_ = smearing_type;
}

void DFPT_Metal::dfdeps(const ModuleBase::matrix& eig, double efermi, 
                         ModuleBase::matrix& dfdeps) {
    // C4 interface reservation: the metallic DFPT branch (smearing
    // derivatives, Fermi-level shift dmu and the occupation-response part of
    // the first-order density) is intentionally NOT implemented in this
    // design-phase iteration; insulating systems only. The data members
    // (sigma_, smearing_type_) and the is_metal_/dmu_ slots of DFPT_PW_Data
    // are already in place for the future implementation.
    (void)eig;
    (void)efermi;
    (void)dfdeps;
    ModuleBase::WARNING_QUIT("DFPT_Metal",
                             "metallic DFPT (dfdeps) is not supported in the design phase");
}

void DFPT_Metal::compute_dmu(int q_idx, const psi::Psi<std::complex<double>>& psi,
                             const ModuleBase::matrix& wg, const ModuleBase::matrix& dfdeps,
                             DFPT_PW_Data& data) {
    (void)q_idx;
    (void)psi;
    (void)wg;
    (void)dfdeps;
    (void)data;
    ModuleBase::WARNING_QUIT("DFPT_Metal",
                             "metallic DFPT (compute_dmu) is not supported in the design phase");
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
    ModuleBase::WARNING_QUIT("DFPT_Metal",
                             "metallic DFPT (compute_drho_metal) is not supported in the design phase");
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
