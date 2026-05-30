#ifndef SOURCE_LCAO_MODULE_EXTRAP_WF_ORTHONORMALIZE_LCAO_H
#define SOURCE_LCAO_MODULE_EXTRAP_WF_ORTHONORMALIZE_LCAO_H

#include "source_lcao/module_extrap/wf_extrap_method.h"
#include "source_psi/psi.h"

namespace ModuleExtrap
{

struct WfOrthonormalizeResult
{
    WfcExtrapStatus status = WfcExtrapStatus::Success;
    double max_deviation = 0.0;
    int failed_state = -1;

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
 * This is a dense serial utility for the first use_prev_wf implementation.  The
 * distributed/k-point path should be implemented separately instead of silently
 * reusing this routine for block-cyclic matrices.
 */
WfOrthonormalizeResult reorthonormalize_gamma_lcao(const double* overlap,
                                                   psi::Psi<double>& wfc,
                                                   double pivot_threshold = 1.0e-14,
                                                   double check_tolerance = 1.0e-8);

} // namespace ModuleExtrap

#endif // SOURCE_LCAO_MODULE_EXTRAP_WF_ORTHONORMALIZE_LCAO_H
