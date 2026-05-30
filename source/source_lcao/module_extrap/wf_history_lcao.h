#ifndef SOURCE_LCAO_MODULE_EXTRAP_WF_HISTORY_LCAO_H
#define SOURCE_LCAO_MODULE_EXTRAP_WF_HISTORY_LCAO_H

#include "source_lcao/module_extrap/wf_extrap_method.h"
#include "source_lcao/module_extrap/wf_snapshot_lcao.h"

#include <cstddef>
#include <deque>

namespace ModuleExtrap
{

template <typename TK>
class WfHistoryLCAO
{
  public:
    explicit WfHistoryLCAO(WfcExtrapMethod method = WfcExtrapMethod::None, std::size_t max_depth = 1);

    WfcExtrapMethod method() const noexcept;
    bool enabled() const noexcept;

    std::size_t size() const noexcept;
    std::size_t max_depth() const noexcept;
    bool empty() const noexcept;

    void set_method(WfcExtrapMethod method) noexcept;
    void set_max_depth(std::size_t max_depth);
    void clear() noexcept;

    void update_after_scf(const int istep, const psi::Psi<TK>& psi, const ModuleBase::matrix& wg);

    // Return a copy of the most recent snapshot.  This avoids exposing an internal
    // pointer/reference whose lifetime would be invalidated by the next history update.
    bool latest_snapshot(WfSnapshotLCAO<TK>& snapshot) const;

  private:
    void prune_history();

    WfcExtrapMethod method_ = WfcExtrapMethod::None;
    std::size_t max_depth_ = 1;
    std::deque<WfSnapshotLCAO<TK>> snapshots_;
};

} // namespace ModuleExtrap

#endif // SOURCE_LCAO_MODULE_EXTRAP_WF_HISTORY_LCAO_H
