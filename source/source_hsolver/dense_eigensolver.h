#ifndef HSOLVER_DENSE_EIGENSOLVER_H
#define HSOLVER_DENSE_EIGENSOLVER_H

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
namespace dense
{
// Dense generalized Hermitian eigensolvers, with 1 <= wanted <= dimension.
// Column-major buffers hold dimension^2 input entries, dimension values
// and dimension*wanted eigenvector coefficients.
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

/** ScaLAPACK generalized eigensolver with reusable block-cyclic storage and workspace.
 * Calls are collective; rank zero supplies the matrices and receives eigenvectors.
 * The borrowed communicator stays alive until this object is destroyed or rebound.
 * A positive block_size sets the block width; zero selects it automatically.
 * Buffers retain their largest allocated capacity across solves.
 */
class ParallelDenseEigensolver
{
  public:
#ifdef __MPI
    ParallelDenseEigensolver();
    ~ParallelDenseEigensolver();
    ParallelDenseEigensolver(const ParallelDenseEigensolver&) = delete;
    ParallelDenseEigensolver& operator=(const ParallelDenseEigensolver&) = delete;

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

/** CUDA generalized Hermitian eigensolver with reusable cuSOLVER resources.
 * Inputs and outputs are host, column-major complex<double>. The first wanted
 * eigenpairs are returned, with vector leading dimension n. The upper triangle
 * defines each input matrix, whose contents are preserved.
 * Resources are allocated on the first solve and reused across subsequent calls.
 * Calls on an instance are serialized on one CUDA device, which remains current
 * through destruction.
 */
class GpuDenseEigensolver
{
  public:
    GpuDenseEigensolver();
    ~GpuDenseEigensolver();
    GpuDenseEigensolver(const GpuDenseEigensolver&) = delete;
    GpuDenseEigensolver& operator=(const GpuDenseEigensolver&) = delete;

    static bool available();

    // Solve H C = S C Lambda (type 1).
    // Use HEGVDX for a partial spectrum, HEGVD for the full spectrum.
    // False means the overlap is not positive definite; other numerical or
    // API failures throw.
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

} // namespace dense
} // namespace hsolver
#endif
