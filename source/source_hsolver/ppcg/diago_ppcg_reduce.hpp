#include "source_base/parallel_reduce.h"

namespace hsolver {
namespace {

template <typename Value>
void reduce_pool_if_mpi_ready(Value& value)
{
#ifdef __MPI
    int initialized = 0;
    int finalized = 0;
    MPI_Initialized(&initialized);
    MPI_Finalized(&finalized);
    if (initialized && !finalized)
        Parallel_Reduce::reduce_pool(value);
#endif
}

template <typename Value>
void reduce_pool_if_mpi_ready(Value* value, const int n)
{
#ifdef __MPI
    int initialized = 0;
    int finalized = 0;
    MPI_Initialized(&initialized);
    MPI_Finalized(&finalized);
    if (initialized && !finalized)
        Parallel_Reduce::reduce_pool(value, n);
#endif
}

template <typename T, typename Real>
Real max_generalized_residual(
    const T* hpsi,
    const T* spsi,
    const Real* eigenvalue,
    int ld,
    int n_dim,
    int ncol)
{
    Real max_res = 0;
    std::vector<double> nrm2_all(ncol, 0.0);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (n_dim * ncol > 4096)
#endif
    for (int j = 0; j < ncol; ++j)
    {
        double nrm2 = 0.0;
        for (int ig = 0; ig < n_dim; ++ig)
        {
            const T r = hpsi[ig + j * ld] - T(eigenvalue[j]) * spsi[ig + j * ld];
            nrm2 += static_cast<double>(std::norm(r));
        }
        nrm2_all[j] = nrm2;
    }
    reduce_pool_if_mpi_ready(nrm2_all.data(), ncol);
    for (int j = 0; j < ncol; ++j)
    {
        max_res = std::max(max_res, std::sqrt(static_cast<Real>(nrm2_all[j])));
    }
    return max_res;
}

template <typename T>
inline void set_zero(std::vector<T>& x)
{
    std::fill(x.begin(), x.end(), T(0));
}

} // anonymous namespace
} // namespace hsolver
