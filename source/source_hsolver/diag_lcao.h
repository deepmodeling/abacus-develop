#ifndef HSOLVER_DIAG_LCAO_H
#define HSOLVER_DIAG_LCAO_H

#include <string>

class Parallel_Blkcyclic;

namespace ModuleBase
{
template <typename T> struct MatrixBlock;
}

namespace psi
{
template <typename T, typename Device> class Psi;
}

namespace base_device
{
struct DEVICE_CPU;
}

namespace hsolver
{
/** Global problem dimensions and the layout used by the selected backend. */
struct LcaoDiagOptions
{
    const std::string& method; // Borrowed backend name.
    int order;
    int bands;
    bool use_gpu;
    const Parallel_Blkcyclic* pool; // Non-null for a temporary k-point pool.
    // True uses the backend's persistent ELPA factorization state for the
    // caller's S buffer. The caller resets that state when S is rebuilt.
    // False starts a fresh factorization of a scratch S and keeps its state
    // local to this solve, as needed by the Löwdin overlap cutoff.
    bool reuse_overlap_factor;
};

/** Solve H C = S C e using the configured LCAO backend and its optional k-point pool. */
template <typename T>
void diag_lcao(ModuleBase::MatrixBlock<T>& h,
               ModuleBase::MatrixBlock<T>& s,
               psi::Psi<T, base_device::DEVICE_CPU>& psi,
               double* eigenvalues,
               const LcaoDiagOptions& options);
} // namespace hsolver

#endif
