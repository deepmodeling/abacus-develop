/**
 * @file dfpt_irrep_data.cpp
 * @brief Implementation of the irrep-indexed DFPT data adapter (Phase 4).
 * @author Mohan Chen (added on 2026-05-18)
 * @note Phase 4 (DFPT wiring) interim layer: simulates the option-2.
 *       signatures on top of the current per-q storage. The irrep index is
 *       forwarded to the underlying per-q accessors (an empty, fully
 *       symmetric placeholder irrep), matching the one-irrep-per-q behavior
 *       currently produced by ModuleCell::QList::get_irreps().
 */
#include "dfpt_irrep_data.h"

namespace ModuleDFPT {

DFPT_IrrepData::DFPT_IrrepData(DFPT_PW_Data& data) : data_(data) {}

int DFPT_IrrepData::get_nq() const
{
    return data_.get_nq();
}

int DFPT_IrrepData::get_nirr(int q_idx) const
{
    return data_.get_nirr(q_idx);
}

std::vector<int> DFPT_IrrepData::get_irrep_modes(int q_idx, int irrep) const
{
    return data_.get_irrep_modes(q_idx, irrep);
}

void DFPT_IrrepData::set_dpsi(int q_idx, int irrep, int k_idx, int band_idx,
                              const std::vector<std::complex<double>>& psi)
{
    (void)irrep;
    data_.set_dpsi(q_idx, k_idx, band_idx, psi);
}

std::vector<std::complex<double>> DFPT_IrrepData::get_dpsi(int q_idx, int irrep, int k_idx,
                                                           int band_idx) const
{
    (void)irrep;
    return data_.get_dpsi(q_idx, k_idx, band_idx);
}

void DFPT_IrrepData::set_drho_r(int q_idx, int irrep, int spin, const std::vector<double>& rho)
{
    (void)irrep;
    data_.set_drho_r(q_idx, spin, rho);
}

std::vector<double> DFPT_IrrepData::get_drho_r(int q_idx, int irrep, int spin) const
{
    (void)irrep;
    return data_.get_drho_r(q_idx, spin);
}

void DFPT_IrrepData::set_drho_g(int q_idx, int irrep, int spin,
                                const std::vector<std::complex<double>>& rho)
{
    (void)irrep;
    data_.set_drho_g(q_idx, spin, rho);
}

std::vector<std::complex<double>> DFPT_IrrepData::get_drho_g(int q_idx, int irrep, int spin) const
{
    (void)irrep;
    return data_.get_drho_g(q_idx, spin);
}

void DFPT_IrrepData::set_dv_r(int q_idx, int irrep, int spin, const std::vector<double>& v)
{
    (void)irrep;
    data_.set_dv_r(q_idx, spin, v);
}

std::vector<double> DFPT_IrrepData::get_dv_r(int q_idx, int irrep, int spin) const
{
    (void)irrep;
    return data_.get_dv_r(q_idx, spin);
}

void DFPT_IrrepData::set_converged(int q_idx, int irrep, bool flag)
{
    converged_[std::make_pair(q_idx, irrep)] = flag;
}

bool DFPT_IrrepData::get_converged(int q_idx, int irrep) const
{
    std::map<std::pair<int, int>, bool>::const_iterator it
        = converged_.find(std::make_pair(q_idx, irrep));
    return it != converged_.end() ? it->second : false;
}

void DFPT_IrrepData::add_residual(int q_idx, int irrep, double r)
{
    residuals_[std::make_pair(q_idx, irrep)].push_back(r);
}

std::vector<double> DFPT_IrrepData::get_residuals(int q_idx, int irrep) const
{
    std::map<std::pair<int, int>, std::vector<double>>::const_iterator it
        = residuals_.find(std::make_pair(q_idx, irrep));
    return it != residuals_.end() ? it->second : std::vector<double>();
}

void DFPT_IrrepData::set_current_iter(int q_idx, int irrep, int iter)
{
    current_iter_[std::make_pair(q_idx, irrep)] = iter;
}

int DFPT_IrrepData::get_current_iter(int q_idx, int irrep) const
{
    std::map<std::pair<int, int>, int>::const_iterator it
        = current_iter_.find(std::make_pair(q_idx, irrep));
    return it != current_iter_.end() ? it->second : 0;
}

} // namespace ModuleDFPT