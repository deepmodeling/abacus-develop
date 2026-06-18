// ============================================================
// This code is added by Mohan Chen on 2026-05-18.
// This code is currently in the design phase and has not been
// put into production yet. It may change in the future.
// Please use this code with caution. Only developers who know
// what they are doing should use this code.
// ============================================================

#include "dfpt_stern.h"

namespace ModuleDFPT {

DFPT_Stern::DFPT_Stern() {}

DFPT_Stern::~DFPT_Stern() {}

void DFPT_Stern::init(int nk, int nbands, int npw_max, const ModuleBase::matrix& eig, 
                      const ModuleBase::matrix& wg, double alpha) {
    nk_ = nk;
    nbands_ = nbands;
    npw_max_ = npw_max;
    eig_ = eig;
    wg_ = wg;
    alpha_ = alpha;
}

void DFPT_Stern::solve(const psi::Psi<std::complex<double>>& psi, 
                       const std::vector<std::complex<double>>& dv_psi,
                       int q_idx, int k_idx, int band_idx, double omega, 
                       DFPT_PW_Data& data) {
    (void)psi;
    (void)dv_psi;
    (void)q_idx;
    (void)k_idx;
    (void)band_idx;
    (void)omega;
    (void)data;
}

void DFPT_Stern::apply_pv(const std::vector<std::complex<double>>& x, 
                          std::vector<std::complex<double>>& px) {
    (void)x;
    (void)px;
}

void DFPT_Stern::apply_op(const std::vector<std::complex<double>>& x,
                          std::vector<std::complex<double>>& y,
                          int k_idx, int band_idx) {
    (void)x;
    (void)y;
    (void)k_idx;
    (void)band_idx;
}

void DFPT_Stern::cg_solve(const std::vector<std::complex<double>>& b,
                          std::vector<std::complex<double>>& x,
                          int k_idx, int band_idx, double& residual) {
    (void)b;
    (void)x;
    (void)k_idx;
    (void)band_idx;
    (void)residual;
}

} // namespace ModuleDFPT