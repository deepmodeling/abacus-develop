#ifndef HSOLVER_DIAG_SK_CUT_H
#define HSOLVER_DIAG_SK_CUT_H

#include "diag_lcao.h"

#include <iosfwd>

namespace hsolver
{
/**
 * Solve H C = S C e in the retained overlap subspace and return C in the AO layout.
 * Use the original generalized problem when all modes are retained. Preserve the
 * unfactored overlap across successive calls.
 */
template <typename T>
void diag_sk_cut(ModuleBase::MatrixBlock<T>& h,
                const ModuleBase::MatrixBlock<T>& s,
                psi::Psi<T, base_device::DEVICE_CPU>& psi,
                double* eigenvalues,
                const LcaoDiagOptions& options,
                int kpoint,
                double cutoff,
                std::ostream* log);
} // namespace hsolver

#endif
