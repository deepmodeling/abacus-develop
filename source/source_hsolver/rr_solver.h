#ifndef HSOLVER_RR_SOLVER_H
#define HSOLVER_RR_SOLVER_H

#include "source_base/module_container/ATen/core/tensor.h"
#include <complex>
#include <cstddef>
#include <memory>
#include <vector>
#ifdef __MPI
#include <mpi.h>
#endif

namespace hsolver
{
namespace rr
{
// Generalized Hermitian Rayleigh–Ritz backends, independent of the outer solver.
// Callers supply 1 <= wanted <= dimension and valid column-major buffers:
// dimension^2 input entries, dimension values and dimension*wanted vectors.
// Dimensions and storage sizes must fit the backend integer and allocation types.
// Root-local HEGVX scratch; reuse one instance across consecutive RR solves.
struct SerialGvWorkspace
{
    using Wide = std::complex<double>;
    // Inputs are preserved. Return LAPACK INFO, or a negative workspace/output status.
    int solve(const Wide* h, const Wide* g, int dimension, int wanted, double* values, Wide* vectors);
    // Local owned capacity; an unused workspace contributes zero on nonroot ranks.
    std::size_t storage_bytes() const;

  private:
    int capacity_dimension_ = 0;
    int work_size_ = 0;
    ct::Tensor h_{ct::DataType::DT_COMPLEX_DOUBLE, ct::DeviceType::CpuDevice, {0}};
    ct::Tensor g_{ct::DataType::DT_COMPLEX_DOUBLE, ct::DeviceType::CpuDevice, {0}};
    ct::Tensor work_{ct::DataType::DT_COMPLEX_DOUBLE, ct::DeviceType::CpuDevice, {0}};
    ct::Tensor real_work_{ct::DataType::DT_DOUBLE, ct::DeviceType::CpuDevice, {0}};
    std::vector<int> int_work_;
    std::vector<int> failures_;
};

/** Explicit CPU RR resources. The communicator is borrowed and must stay alive
 * until this object is destroyed or rebound. Calls are collective. Matrix data
 * are read from rank zero and vectors returned there, with no implicit broadcast.
 * The caller supplies a valid communicator and block_size >= 0 (0 selects automatically).
 * One current layout/query plus high-water buffers bounds memory independently
 * of the number of different tail dimensions seen during SCF.
 */
class ParallelRREigensolver
{
  public:
#ifdef __MPI
    ParallelRREigensolver();
    ~ParallelRREigensolver();
    ParallelRREigensolver(const ParallelRREigensolver&) = delete;
    ParallelRREigensolver& operator=(const ParallelRREigensolver&) = delete;

    // Returns wanted eigenvalues and n*wanted column-major coefficients on rank zero.
    void solve_generalized(const std::complex<double>* h,
        const std::complex<double>* s,
        const int n,
        const int bands,
        double* values,
        std::complex<double>* vectors,
        const MPI_Comm& comm,
        const int block_size);

    // Local owned workspace capacity; zero before the first solve.
    std::size_t storage_bytes() const;

  private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
#endif
};

/** CUDA generalized Hermitian RR with caller-owned, capacity-reused resources.
 * Inputs and outputs are host, column-major complex<double>. Only the first
 * wanted eigenpairs are returned, with vector leading dimension n. The input
 * is preserved; the upper triangle defines the Hermitian matrix. Resources
 * are created lazily, so CPU/ROCm callers can own an unused instance safely.
 * The caller keeps the same CUDA device current throughout the instance's use
 * and destruction, and must not use it concurrently.
 */
class GpuRREigensolver
{
  public:
    GpuRREigensolver();
    ~GpuRREigensolver();
    GpuRREigensolver(GpuRREigensolver&&) noexcept;
    GpuRREigensolver(const GpuRREigensolver&) = delete;
    GpuRREigensolver& operator=(const GpuRREigensolver&) = delete;

    static bool available();

    // Type-1 generalized Hermitian RR, with the same host layout and ownership.
    // Use HEGVDX for a partial spectrum, HEGVD for the full spectrum.
    // False means the overlap is not positive definite; other numerical or
    // API failures throw. Inputs are preserved.
    bool solve_generalized(const std::complex<double>* matrix,
        const std::complex<double>* metric,
        const int n,
        const int wanted,
        double* values,
        std::complex<double>* vectors);

    // Local owned workspace capacity; zero before the first solve.
    std::size_t storage_bytes() const;

  private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace rr
} // namespace hsolver
#endif
