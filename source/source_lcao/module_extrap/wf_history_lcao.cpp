#include "source_lcao/module_extrap/wf_history_lcao.h"

#include "source_lcao/module_extrap/wf_orthonormalize_lcao.h"
#include "source_base/timer.h"

#include <algorithm>
#include <complex>

namespace ModuleExtrap
{


template <typename TK>
WfHistoryLCAO<TK>::WfHistoryLCAO(const WfcExtrapMethod method, const std::size_t max_depth)
    : method_(method), max_depth_(std::max<std::size_t>(max_depth, 1))
{
}

template <typename TK>
WfcExtrapMethod WfHistoryLCAO<TK>::method() const noexcept
{
    return this->method_;
}

template <typename TK>
bool WfHistoryLCAO<TK>::enabled() const noexcept
{
    return this->method_ != WfcExtrapMethod::None;
}

template <typename TK>
std::size_t WfHistoryLCAO<TK>::size() const noexcept
{
    return this->snapshots_.size();
}

template <typename TK>
std::size_t WfHistoryLCAO<TK>::max_depth() const noexcept
{
    return this->max_depth_;
}

template <typename TK>
bool WfHistoryLCAO<TK>::empty() const noexcept
{
    return this->snapshots_.empty();
}

template <typename TK>
void WfHistoryLCAO<TK>::set_method(const WfcExtrapMethod method) noexcept
{
    this->method_ = method;
}

template <typename TK>
void WfHistoryLCAO<TK>::set_max_depth(const std::size_t max_depth)
{
    this->max_depth_ = std::max<std::size_t>(max_depth, 1);
    this->prune_history();
}

template <typename TK>
void WfHistoryLCAO<TK>::clear() noexcept
{
    this->snapshots_.clear();
}

template <typename TK>
void WfHistoryLCAO<TK>::update_after_scf(const int istep, const psi::Psi<TK>& psi, const ModuleBase::matrix& wg)
{
    if (!this->enabled())
    {
        return;
    }

    this->snapshots_.emplace_back(istep, psi, wg);
    this->prune_history();
}

template <typename TK>
bool WfHistoryLCAO<TK>::latest_snapshot(WfSnapshotLCAO<TK>& snapshot) const
{
    if (this->snapshots_.empty())
    {
        return false;
    }
    snapshot = this->snapshots_.back();
    return true;
}

template <typename TK>
WfExtrapApplyResult WfHistoryLCAO<TK>::try_use_prev_wf_gamma(const TK*,
                                                             const Parallel_Orbitals&,
                                                             psi::Psi<TK>&,
                                                             const ModuleBase::matrix&,
                                                             const double,
                                                             const double)
{
    WfExtrapApplyResult result;
    result.status = WfcExtrapStatus::Unsupported;
    return result;
}

template <>
WfExtrapApplyResult WfHistoryLCAO<double>::try_use_prev_wf_gamma(const double* current_overlap,
                                                                 const Parallel_Orbitals& pv,
                                                                 psi::Psi<double>& psi,
                                                                 const ModuleBase::matrix& wg_now,
                                                                 const double pivot_threshold,
                                                                 const double check_tolerance)
{
    WfExtrapApplyResult result;

    if (this->method_ != WfcExtrapMethod::UsePrevWf)
    {
        result.status = WfcExtrapStatus::Disabled;
        return result;
    }

    if (this->snapshots_.empty())
    {
        result.status = WfcExtrapStatus::EmptyHistory;
        return result;
    }

    if (current_overlap == nullptr || !(pivot_threshold >= 0.0) || !(check_tolerance > 0.0))
    {
        result.status = WfcExtrapStatus::InvalidInput;
        return result;
    }

    // The first PR only supports the real Gamma-only path. In this path the
    // first Psi dimension may label spin channels, but not physical k-points.
    if (!psi.get_k_first())
    {
        result.status = WfcExtrapStatus::Unsupported;
        return result;
    }

    const WfSnapshotLCAO<double>& snapshot = this->snapshots_.back();
    result.snapshot_istep = snapshot.istep;

    if (!snapshot.compatible_with(psi, wg_now))
    {
        result.status = WfcExtrapStatus::DimensionMismatch;
        return result;
    }

    // Work on an owned trial Psi first. The caller's Psi is not overwritten
    // unless the stored WFN can be loaded and reorthonormalized successfully.
    psi::Psi<double> psi_trial(psi);
    ModuleBase::timer::start("WFN_Extrap", "restore_snapshot");
    const bool snapshot_loaded = snapshot.load_to(psi_trial);
    ModuleBase::timer::end("WFN_Extrap", "restore_snapshot");
    if (!snapshot_loaded)
    {
        result.status = WfcExtrapStatus::DimensionMismatch;
        return result;
    }

    ModuleBase::timer::start("WFN_Extrap", "orthonormalize");
    const WfOrthonormalizeResult orth_result = reorthonormalize_gamma_lcao(current_overlap,
                                                                           pv,
                                                                           psi_trial,
                                                                           snapshot.wg,
                                                                           1.0e-12,
                                                                           pivot_threshold,
                                                                           check_tolerance);
    ModuleBase::timer::end("WFN_Extrap", "orthonormalize");
    if (!orth_result.ok())
    {
        result.status = orth_result.status;
        result.failed_state = orth_result.failed_state;
        result.failed_pivot_index = orth_result.failed_pivot_index;
        result.failed_pivot = orth_result.failed_pivot;
        result.min_metric_diag = orth_result.min_metric_diag;
        result.max_metric_diag = orth_result.max_metric_diag;
        result.max_metric_abs = orth_result.max_metric_abs;
        result.max_metric_asymmetry = orth_result.max_metric_asymmetry;
        result.nstate = orth_result.nstate;
        result.nbands = orth_result.nbands;
        result.nbasis = orth_result.nbasis;
        result.nactive_bands = orth_result.nactive_bands;
        result.max_orthonormality_deviation = orth_result.max_deviation;
        return result;
    }

    psi = psi_trial;
    result.status = WfcExtrapStatus::Success;
    result.failed_state = -1;
    result.nstate = orth_result.nstate;
    result.nbands = orth_result.nbands;
    result.nbasis = orth_result.nbasis;
    result.nactive_bands = orth_result.nactive_bands;
    result.max_orthonormality_deviation = orth_result.max_deviation;
    return result;
}


template <typename TK>
void WfHistoryLCAO<TK>::prune_history()
{
    while (this->snapshots_.size() > this->max_depth_)
    {
        this->snapshots_.pop_front();
    }
}

template class WfHistoryLCAO<double>;
template class WfHistoryLCAO<std::complex<double>>;

} // namespace ModuleExtrap
