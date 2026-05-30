#include "source_lcao/module_extrap/wf_history_lcao.h"

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
