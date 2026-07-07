#ifndef SOURCE_LCAO_MODULE_EXTRAP_WF_HISTORY_LCAO_H
#define SOURCE_LCAO_MODULE_EXTRAP_WF_HISTORY_LCAO_H

#include "source_lcao/module_extrap/wf_extrap_method.h"
#include "source_lcao/module_extrap/wf_snapshot_lcao.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include <cstddef>
#include <deque>
#include <string>

namespace hamilt
{
template <typename TK, typename TR>
class HamiltLCAO;
}

namespace elecstate
{
template <typename TK, typename TR>
class DensityMatrix;
}

class Charge;

namespace ModuleExtrap
{

struct WfExtrapApplyResult
{
    WfcExtrapStatus status = WfcExtrapStatus::Disabled;
    double max_orthonormality_deviation = 0.0;

    // Propagated diagnostics from WfOrthonormalizeResult.
    double failed_pivot = 0.0;
    double min_metric_diag = 0.0;
    double max_metric_diag = 0.0;
    double max_metric_abs = 0.0;
    double max_metric_asymmetry = 0.0;
    int failed_state = -1;
    int failed_pivot_index = -1;
    int nstate = 0;
    int nbands = 0;
    int nbasis = 0;
    int nactive_bands = 0;
    int snapshot_istep = -1;

    bool ok() const noexcept
    {
        return this->status == WfcExtrapStatus::Success;
    }
};

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

    /**
     * Restore the latest Gamma-only NAO wavefunction and reorthonormalize it with
     * the current overlap matrix.
     *
     * This method only updates the Psi object after the full restore and
     * reorthonormalization path succeeds.  DMK/DMR/rho rebuilding is deliberately
     * left to the ESolver integration layer, where the existing density-matrix code
     * can be called with the same control flow as the current initialization path.
     *
     * The first PR only supports the real Gamma-only path.  Multi-k support requires
     * complex WFN snapshots, per-k S(k), and phase/convention handling, and should be
     * added in a separate path instead of silently reusing this helper.
     */
    WfExtrapApplyResult try_use_prev_wf_gamma(const TK* current_overlap,
                                              const Parallel_Orbitals& pv,
                                              psi::Psi<TK>& psi,
                                              const ModuleBase::matrix& wg_now,
                                              double pivot_threshold = 1.0e-14,
                                              double check_tolerance = 1.0e-8);

  private:
    void prune_history();

    WfcExtrapMethod method_ = WfcExtrapMethod::None;
    std::size_t max_depth_ = 1;
    std::deque<WfSnapshotLCAO<TK>> snapshots_;
};

/**
 * Restore a Gamma-only WFN snapshot and rebuild the corresponding density.
 *
 * This runtime integration is kept outside WfHistoryLCAO so the history and
 * orthonormalization core remains independently unit-testable.
 */
template <typename TK, typename TR>
bool initialize_gamma_density_from_history(WfHistoryLCAO<TK>& history,
                                           hamilt::HamiltLCAO<TK, TR>& hamiltonian,
                                           const Parallel_Orbitals& pv,
                                           psi::Psi<TK>& psi,
                                           const ModuleBase::matrix& wg_now,
                                           elecstate::DensityMatrix<TK, double>& dmat,
                                           Charge& charge,
                                           int nspin,
                                           const std::string& ks_solver);

} // namespace ModuleExtrap

#endif // SOURCE_LCAO_MODULE_EXTRAP_WF_HISTORY_LCAO_H
