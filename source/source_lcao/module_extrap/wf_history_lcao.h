#ifndef SOURCE_LCAO_MODULE_EXTRAP_WF_HISTORY_LCAO_H
#define SOURCE_LCAO_MODULE_EXTRAP_WF_HISTORY_LCAO_H

#include "source_lcao/module_extrap/wf_extrap_status.h"
#include "source_lcao/module_extrap/wf_snapshot_lcao.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_hamilt/hamilt.h"
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
    WfcExtrapStatus status = WfcExtrapStatus::InvalidInput;
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
    WfHistoryLCAO() = default;

    void update_after_scf(const int istep, const psi::Psi<TK>& psi, const ModuleBase::matrix& wg);

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
                                              double pivot_threshold,
                                              double check_tolerance);

    void initialize_gamma_density(hamilt::Hamilt<TK>& hamiltonian,
                                  const Parallel_Orbitals& pv,
                                  psi::Psi<TK>& psi,
                                  const ModuleBase::matrix& wg_now,
                                  elecstate::DensityMatrix<TK, double>& dmat,
                                  Charge& charge);

  private:
    bool has_snapshot_ = false;
    WfSnapshotLCAO<TK> snapshot_;
};

} // namespace ModuleExtrap

#endif // SOURCE_LCAO_MODULE_EXTRAP_WF_HISTORY_LCAO_H
