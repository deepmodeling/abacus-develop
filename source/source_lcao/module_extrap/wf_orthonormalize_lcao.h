#ifndef SOURCE_LCAO_MODULE_EXTRAP_WF_ORTHONORMALIZE_LCAO_H
#define SOURCE_LCAO_MODULE_EXTRAP_WF_ORTHONORMALIZE_LCAO_H

#include "source_lcao/module_extrap/wf_extrap_method.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_base/matrix.h"
#include "source_psi/psi.h"

namespace ModuleExtrap
{

struct WfOrthonormalizeResult
{
    WfcExtrapStatus status = WfcExtrapStatus::Success;
    double max_deviation = 0.0;

    // Diagnostics for orthogonalization failures.  These values are kept in the
    // result object so that the ESolver layer can print them without exposing
    // implementation details of the orthonormalization helper.
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

    bool ok() const noexcept
    {
        return this->status == WfcExtrapStatus::Success;
    }
};

/**
 * Reorthonormalize Gamma-only NAO wavefunctions with the current overlap matrix.
 *
 * The input wavefunction coefficients are assumed to be stored as band-major rows,
 * i.e. C(iband, ibasis) for the current Gamma/spin channel.  For each channel this
 * routine computes
 *
 *     O = C S C^T
 *
 * and applies a Cholesky-based transformation C_new = L^{-1} C, where O = L L^T.
 * Thus the updated coefficients satisfy C_new S C_new^T = I up to numerical noise.
 *
 * In Gamma-only LCAO, Psi::get_nk() may label spin channels rather than physical
 * k-points.  This routine therefore treats the first Psi dimension as a generic
 * state/channel index.
 *
 * Only bands with non-negligible occupation are reorthonormalized.  For the
 * initial-density use case of use_prev_wf, completely unoccupied bands do not
 * contribute to DMK/DMR/rho and orthonormalizing them can make the metric nearly
 * singular when all basis-sized bands are present.
 *
 * The implementation is MPI-aware for the Gamma-only block-cyclic NAO layout.
 * It assembles dense global work matrices from distributed local blocks, performs
 * the small Cholesky transformation redundantly on all ranks, and scatters the
 * transformed coefficients back to the local Psi layout.  This favors correctness
 * and a simple reviewable backend over peak performance; a future optimized path
 * can replace the dense work arrays with fully distributed PBLAS operations.
 */
WfOrthonormalizeResult reorthonormalize_gamma_lcao(const double* overlap,
                                                   const Parallel_Orbitals& pv,
                                                   psi::Psi<double>& wfc,
                                                   const ModuleBase::matrix& occupations,
                                                   double occupation_threshold = 1.0e-12,
                                                   double pivot_threshold = 1.0e-14,
                                                   double check_tolerance = 1.0e-8);

} // namespace ModuleExtrap

#endif // SOURCE_LCAO_MODULE_EXTRAP_WF_ORTHONORMALIZE_LCAO_H
