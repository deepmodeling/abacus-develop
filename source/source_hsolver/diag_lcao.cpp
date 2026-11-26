#include "diag_lcao.h"

#include "diago_lapack.h"
#include "source_base/parallel_blkcyc.h"
#include "source_base/tool_quit.h"

#ifdef __MPI
#include "diago_scalapack.h"
#endif
#ifdef __ELPA
#include "diago_elpa.h"
#include "diago_elpa_native.h"
#endif
#ifdef __CUDA
#include "diago_cusolver.h"
#ifdef __CUSOLVERMP
#include "diago_cusolvermp.h"
#endif
#endif

#include <complex>

namespace hsolver
{
namespace
{
template <typename Solver, typename T>
inline void run_diagonalizer(Solver& solver,
                             ModuleBase::MatrixBlock<T>& h,
                             ModuleBase::MatrixBlock<T>& s,
                             psi::Psi<T>& psi,
                             double* eigenvalues,
                             const Parallel_Blkcyclic* pool)
{
#ifdef __MPI
    if (pool != nullptr)
    {
        MPI_Comm pool_comm = pool->comm();
        solver.diag_pool(h, s, psi, eigenvalues, pool_comm);
        return;
    }
#endif
    solver.diag(h, s, psi, eigenvalues);
}
} // namespace

template <typename T>
void diag_lcao(ModuleBase::MatrixBlock<T>& h,
               ModuleBase::MatrixBlock<T>& s,
               psi::Psi<T>& psi,
               double* eigenvalues,
               const LcaoDiagOptions& options)
{
    const auto& method = options.method;
    const int order = options.order;
    const int bands = options.bands;
#ifdef __MPI
    if (method == "scalapack_gvx")
    {
        DiagoScalapack<T> solver(order, bands);
        run_diagonalizer(solver, h, s, psi, eigenvalues, options.pool);
    }
    else
#endif
#ifdef __ELPA
    // ELPA updates its factorization flag in place. For a scratch overlap,
    // start at zero and keep the update local so it cannot mark the caller's
    // overlap as factored or make a later scratch solve reuse a stale factor.
    if (method == "genelpa")
    {
        int workspace_factor = 0;
        auto& overlap_factor = options.reuse_overlap_factor ? DiagoElpa<T>::DecomposedState : workspace_factor;
        DiagoElpa<T> solver(order, bands, overlap_factor);
        run_diagonalizer(solver, h, s, psi, eigenvalues, options.pool);
    }
    else if (method == "elpa")
    {
        int workspace_factor = 0;
        auto& overlap_factor = options.reuse_overlap_factor ? DiagoElpaNative<T>::DecomposedState : workspace_factor;
        DiagoElpaNative<T> solver(order, bands, options.use_gpu, overlap_factor);
        run_diagonalizer(solver, h, s, psi, eigenvalues, options.pool);
    }
    else
#endif
#ifdef __CUDA
    if (method == "cusolver")
    {
        DiagoCusolver<T> solver(order, bands);
        solver.diag(h, s, psi, eigenvalues);
    }
    else
#ifdef __CUSOLVERMP
    if (method == "cusolvermp")
    {
        DiagoCusolverMP<T> solver(order, bands);
        solver.diag(h, s, psi, eigenvalues);
    }
    else
#endif
#endif
    if (method == "lapack")
    {
        DiagoLapack<T> solver(order, bands);
        run_diagonalizer(solver, h, s, psi, eigenvalues, options.pool);
    }
    else
    {
        const std::string message = "Unsupported LCAO eigensolver: " + method;
        ModuleBase::WARNING_QUIT("diag_lcao", message);
    }
}

template void diag_lcao<double>(ModuleBase::MatrixBlock<double>&,
                                ModuleBase::MatrixBlock<double>&,
                                psi::Psi<double>&,
                                double*,
                                const LcaoDiagOptions&);
template void diag_lcao<std::complex<double>>(ModuleBase::MatrixBlock<std::complex<double>>&,
                                              ModuleBase::MatrixBlock<std::complex<double>>&,
                                              psi::Psi<std::complex<double>>&,
                                              double*,
                                              const LcaoDiagOptions&);
} // namespace hsolver
