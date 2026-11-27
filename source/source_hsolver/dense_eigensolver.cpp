#include "dense_eigensolver.h"
#include "source_base/memory_recorder.h"
#include "source_base/module_container/base/third_party/lapack.h"
#include "source_base/timer.h"
#ifdef __MPI
#include "source_base/module_external/blacs_connector.h"
#include "source_base/module_external/scalapack_connector.h"
#include "source_base/parallel_blkcyc.h"
#endif
#ifdef __CUDA
#include <cuda_runtime_api.h>
#include <cusolverDn.h>
#endif
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>

namespace hsolver
{
namespace dense
{
std::size_t SerialGvWorkspace::storage_bytes() const
{
    return h_.buffer().GetAllocatedBytes() + g_.buffer().GetAllocatedBytes() + work_.buffer().GetAllocatedBytes() +
           real_work_.buffer().GetAllocatedBytes() + int_work_.capacity() * sizeof(int) +
           failures_.capacity() * sizeof(int);
}

int SerialGvWorkspace::solve(
    const Wide* h, const Wide* g, const int dimension, const int wanted, double* values, Wide* vectors)
{
    ModuleBase::timer::Scope timer("SerialGvWorkspace", "solve");
    const int64_t square = static_cast<int64_t>(dimension) * dimension;
    const int type = 1;
    const char job = 'V';
    const char range = 'I';
    const char triangle = 'U';
    const double bound = 0;
    const double tolerance = 0;
    const int first = 1;
    int found = 0;
    int info = 0;
    if (dimension > capacity_dimension_)
    {
        h_.resize({square});
        g_.resize({square});
        const int64_t real_work_size = static_cast<int64_t>(7) * dimension;
        real_work_.resize({real_work_size});
        int_work_.resize(static_cast<std::size_t>(5) * dimension);
        failures_.resize(dimension);
        const int query_size = -1;
        Wide query(0);
        ct::lapackConnector::hegvx(type,
            job,
            range,
            triangle,
            dimension,
            h_.data<Wide>(),
            dimension,
            g_.data<Wide>(),
            dimension,
            bound,
            bound,
            first,
            wanted,
            tolerance,
            &found,
            values,
            vectors,
            dimension,
            &query,
            query_size,
            real_work_.data<double>(),
            int_work_.data(),
            failures_.data(),
            info);
        if (info != 0)
        {
            return info;
        }
        const double requested = query.real();
        if (!std::isfinite(requested) || requested > std::numeric_limits<int>::max() ||
            dimension > std::numeric_limits<int>::max() / 2)
        {
            return -1001;
        }
        work_size_ = std::max(2 * dimension, static_cast<int>(requested));
        work_.resize({work_size_});
        capacity_dimension_ = dimension;
    }
    const std::size_t bytes = storage_bytes();
    ModuleBase::Memory::record("RR::serial_gv", bytes, false);
    std::copy_n(h, square, h_.data<Wide>());
    std::copy_n(g, square, g_.data<Wide>());
    ct::lapackConnector::hegvx(type,
        job,
        range,
        triangle,
        dimension,
        h_.data<Wide>(),
        dimension,
        g_.data<Wide>(),
        dimension,
        bound,
        bound,
        first,
        wanted,
        tolerance,
        &found,
        values,
        vectors,
        dimension,
        work_.data<Wide>(),
        work_size_,
        real_work_.data<double>(),
        int_work_.data(),
        failures_.data(),
        info);
    if (info == 0 && found != wanted)
    {
        return -1002;
    }
    return info;
}

#ifdef __MPI
namespace
{
template <typename T>
void grow_workspace(std::vector<T>& buffer, const std::size_t count)
{
    if (count > buffer.size())
    {
        buffer.resize(count);
    }
}

int checked_workspace_size(const double size)
{
    if (!std::isfinite(size) || size < 1 || size > std::numeric_limits<int>::max())
    {
        throw std::runtime_error("Invalid ScaLAPACK workspace size.");
    }
    return static_cast<int>(std::ceil(size));
}
class GvxWorkspace
{
    using T = std::complex<double>;
    using Real = double;

  public:
    GvxWorkspace() = default;
    GvxWorkspace(const GvxWorkspace&) = delete;
    GvxWorkspace& operator=(const GvxWorkspace&) = delete;

    // A BLACS context integer may be recycled after rebinding a communicator.
    void invalidate_layout() { bands_ = 0; }

    std::size_t storage_bytes() const;

    void solve(const int* desc,
        const int nrow,
        const int ncol,
        const int nbands,
        const T* h,
        const T* s,
        double* values,
        T* vectors);

  private:
    std::array<int, 9> descriptor_ = {};
    int bands_ = 0;
    int lwork_ = 0;
    int lrwork_base_ = 0;
    int liwork_ = 0;
    std::vector<T> h_;
    std::vector<T> s_;
    std::vector<T> work_;
    std::vector<Real> rwork_;
    std::vector<Real> gap_;
    std::vector<int> iwork_;
    std::vector<int> ifail_;
    std::vector<int> iclustr_;
};

std::size_t GvxWorkspace::storage_bytes() const
{
    return sizeof(T) * (h_.capacity() + s_.capacity() + work_.capacity()) +
           sizeof(Real) * (rwork_.capacity() + gap_.capacity()) +
           sizeof(int) * (iwork_.capacity() + ifail_.capacity() + iclustr_.capacity());
}

void GvxWorkspace::solve(
    const int* desc, const int nrow, const int ncol, const int nbands, const T* h, const T* s, Real* values, T* vectors)
{
    ModuleBase::timer::Scope timer("GvxWorkspace", "solve");
    ModuleBase::Memory::record("RR::gvx_work", 0, false);
    int nprow = 1;
    int npcol = 1;
    int myprow = 0;
    int mypcol = 0;
    Cblacs_gridinfo(desc[1], &nprow, &npcol, &myprow, &mypcol);
    const int n = desc[2];
    const int dsize = nprow * npcol;
    const std::size_t local_size = static_cast<std::size_t>(desc[8]) * ncol;
    const std::size_t storage_size = std::max<std::size_t>(1, local_size);
    grow_workspace(h_, storage_size);
    grow_workspace(s_, storage_size);
    const int cluster_entries = 2 * dsize;
    ifail_.resize(n);
    iclustr_.resize(cluster_entries);
    gap_.resize(dsize);
    const char jobz = 'V';
    const char range = 'I';
    const char uplo = 'U';
    const int itype = 1;
    const int one = 1;
    const Real abstol = SCALAPACK_ABSTOL;
    const Real orfac = SCALAPACK_ORFAC;
    const Real zero = 0;
    int m = 0;
    int nz = 0;
    int info = 0;
    const bool new_query = bands_ != nbands || !std::equal(descriptor_.begin(), descriptor_.end(), desc);
    if (new_query)
    {
        T work_query[3] = {};
        Real real_query[3] = {};
        int int_query[3] = {};
        int lwork = -1;
        int lrwork = -1;
        int liwork = -1;
        pzhegvx_(&itype,
            &jobz,
            &range,
            &uplo,
            &n,
            h_.data(),
            &one,
            &one,
            desc,
            s_.data(),
            &one,
            &one,
            desc,
            &zero,
            &zero,
            &one,
            &nbands,
            &abstol,
            &m,
            &nz,
            values,
            &orfac,
            vectors,
            &one,
            &one,
            desc,
            work_query,
            &lwork,
            real_query,
            &lrwork,
            int_query,
            &liwork,
            ifail_.data(),
            iclustr_.data(),
            gap_.data(),
            &info);

        if (info != 0)
        {
            throw std::runtime_error("ScaLAPACK GVX query: info=" + std::to_string(info));
        }
        const double work_size = std::real(work_query[0]);
        lwork_ = checked_workspace_size(work_size);
        const double int_size = int_query[0];
        liwork_ = checked_workspace_size(int_size);
        const double real_size = real_query[0];
        lrwork_base_ = checked_workspace_size(real_size);
        const int work_storage = std::max(3, lwork_);
        grow_workspace(work_, work_storage);
        grow_workspace(iwork_, liwork_);
        std::copy(desc, desc + 9, descriptor_.begin());
        bands_ = nbands;
    }

    // PZHEGVX requires query RWORK plus (CLUSTERSIZE-1)*N for cluster
    // reorthogonalization. Reserve the full N-eigenvalue cluster once.
    // Restore this call's original inputs: GVX overwrites both matrices.
    for (int j = 0; j < ncol; ++j)
    {
        for (int i = 0; i < nrow; ++i)
        {
            const std::size_t index = i + static_cast<std::size_t>(j) * desc[8];
            h_[index] = h[index];
            s_[index] = s[index];
        }
    }
    const double real_size = lrwork_base_ + static_cast<double>(n - 1) * n;
    int lrwork = checked_workspace_size(real_size);
    const int real_storage = std::max(3, lrwork);
    grow_workspace(rwork_, real_storage);
    const std::size_t bytes = storage_bytes();
    ModuleBase::Memory::record("RR::gvx_work", bytes, false);
    pzhegvx_(&itype,
        &jobz,
        &range,
        &uplo,
        &n,
        h_.data(),
        &one,
        &one,
        desc,
        s_.data(),
        &one,
        &one,
        desc,
        &zero,
        &zero,
        &one,
        &nbands,
        &abstol,
        &m,
        &nz,
        values,
        &orfac,
        vectors,
        &one,
        &one,
        desc,
        work_.data(),
        &lwork_,
        rwork_.data(),
        &lrwork,
        iwork_.data(),
        &liwork_,
        ifail_.data(),
        iclustr_.data(),
        gap_.data(),
        &info);
    if (info != 0)
    {
        throw std::runtime_error("ScaLAPACK GVX failed: info=" + std::to_string(info));
    }
    if (m != nbands || nz != nbands)
    {
        throw std::runtime_error("ScaLAPACK GVX returned an incomplete eigensystem.");
    }
}

} // namespace

struct ParallelDenseEigensolver::Impl
{
    using Z = std::complex<double>;
    MPI_Comm communicator = MPI_COMM_NULL;
    Parallel_Blkcyclic local;
    Parallel_Blkcyclic root;
    int dimension = 0;
    int block = 0;
    std::vector<Z> h;
    std::vector<Z> s;
    std::vector<Z> vectors;
    GvxWorkspace gvx;

    void prepare(const int n, const int bands, const MPI_Comm& comm, const int block_size)
    {
        int nb = block_size;
        if (nb == 0)
        {
            nb = bands > 500 ? 32 : 16;
        }
        // Each communicator handle gets its own BLACS grid.
        if (comm != communicator || n != dimension || nb != block)
        {
            if (comm != communicator)
            {
                local.init(n, n, nb, comm, false);
                communicator = comm;
                gvx.invalidate_layout();
            }
            else
            {
                local.set(n, n, nb, local.ictxt());
            }
            root.set(n, n, n, local.ictxt());
            dimension = n;
            block = nb;
        }
        const int columns = local.ncloc();
        const std::size_t elements = static_cast<std::size_t>(local.desc()[8]) * columns;
        const std::size_t storage = std::max<std::size_t>(1, elements);
        grow_workspace(h, storage);
        grow_workspace(s, storage);
        grow_workspace(vectors, storage);
    }

    void scatter(const Z* input, Z* output)
    {
        // GEMR2D only reads its source; the connector accepts a mutable pointer.
        Z* source = const_cast<Z*>(input);
        Cpxgemr2d(dimension, dimension, source, 1, 1, root.desc(), output, 1, 1, local.desc(), local.ictxt());
    }

    void gather(Z* output, const int columns)
    {
        Cpxgemr2d(dimension, columns, vectors.data(), 1, 1, local.desc(), output, 1, 1, root.desc(), local.ictxt());
    }
};

ParallelDenseEigensolver::ParallelDenseEigensolver() = default;
ParallelDenseEigensolver::~ParallelDenseEigensolver() = default;

std::size_t ParallelDenseEigensolver::storage_bytes() const
{
    if (!impl_)
    {
        return 0;
    }
    const Impl& p = *impl_;
    return sizeof(Impl::Z) * (p.h.capacity() + p.s.capacity() + p.vectors.capacity()) + p.gvx.storage_bytes();
}

void ParallelDenseEigensolver::solve_generalized(const std::complex<double>* h,
    const std::complex<double>* s,
    const int n,
    const int bands,
    double* values,
    std::complex<double>* vectors,
    const MPI_Comm& comm,
    const int block_size)
{
    ModuleBase::timer::Scope timer("ParallelDenseEigensolver", "solve_generalized");
    ModuleBase::Memory::record("RR::parallel_stage", 0, false);
    if (!impl_)
    {
        impl_.reset(new Impl);
    }
    Impl& p = *impl_;
    p.prepare(n, bands, comm, block_size);
    p.scatter(h, p.h.data());
    p.scatter(s, p.s.data());
    const int rows = p.local.nrloc();
    const int columns = p.local.ncloc();
    p.gvx.solve(p.local.desc(), rows, columns, bands, p.h.data(), p.s.data(), values, p.vectors.data());
    p.gather(vectors, bands);
    const std::size_t stage_bytes =
        (p.h.capacity() + p.s.capacity() + p.vectors.capacity()) * sizeof(std::complex<double>);
    ModuleBase::Memory::record("RR::parallel_stage", stage_bytes, false);
}

#endif // __MPI

namespace
{
#ifdef __CUDA
void check_cuda(const cudaError_t status, const char* operation)
{
    if (status != cudaSuccess)
    {
        throw std::runtime_error(std::string(operation) + ": " + cudaGetErrorString(status));
    }
}

void check_solver(const cusolverStatus_t status, const char* operation)
{
    if (status != CUSOLVER_STATUS_SUCCESS)
    {
        throw std::runtime_error(std::string(operation) + ": cuSOLVER status " + std::to_string(status));
    }
}

#endif
} // namespace

struct GpuDenseEigensolver::Impl
{
#ifdef __CUDA
    cusolverDnHandle_t handle = nullptr;
    int query_n = 0;
    int query_wanted = 0;
    int lwork = 0;
    ct::Tensor matrix{ct::DataType::DT_COMPLEX_DOUBLE, ct::DeviceType::GpuDevice, {0}};
    ct::Tensor metric{ct::DataType::DT_COMPLEX_DOUBLE, ct::DeviceType::GpuDevice, {0}};
    ct::Tensor eigenvalues{ct::DataType::DT_DOUBLE, ct::DeviceType::GpuDevice, {0}};
    ct::Tensor info{ct::DataType::DT_INT, ct::DeviceType::GpuDevice, {0}};
    ct::Tensor work{ct::DataType::DT_COMPLEX_DOUBLE, ct::DeviceType::GpuDevice, {0}};

    std::size_t storage_bytes() const
    {
        return matrix.buffer().GetAllocatedBytes() + metric.buffer().GetAllocatedBytes() +
               eigenvalues.buffer().GetAllocatedBytes() + info.buffer().GetAllocatedBytes() +
               work.buffer().GetAllocatedBytes();
    }

    void reserve(ct::Tensor& buffer, const int64_t count)
    {
        const std::size_t previous = buffer.buffer().GetAllocatedBytes();
        buffer.resize({count});
        // The container allocator returns null on allocation failure.
        if (count > 0 && buffer.data() == nullptr)
        {
            throw std::runtime_error("GPU RR allocation failed.");
        }
        const std::size_t capacity = buffer.buffer().GetAllocatedBytes();
        if (capacity > previous)
        {
            ModuleBase::Memory::record("RR::gpu_work", 0, false);
            const std::size_t bytes = storage_bytes();
            ModuleBase::Memory::record_gpu("RR::gpu_work", bytes, false);
        }
    }

    Impl()
    {
        const cusolverStatus_t status = cusolverDnCreate(&handle);
        check_solver(status, "GPU RR handle creation");
    }

    ~Impl() { cusolverDnDestroy(handle); }
#endif
};

GpuDenseEigensolver::GpuDenseEigensolver() = default;
GpuDenseEigensolver::~GpuDenseEigensolver() = default;

bool GpuDenseEigensolver::available()
{
#ifdef __CUDA
    return true;
#else
    return false;
#endif
}

std::size_t GpuDenseEigensolver::storage_bytes() const
{
#ifdef __CUDA
    if (impl_)
    {
        return impl_->storage_bytes();
    }
#endif
    return 0;
}

bool GpuDenseEigensolver::solve_generalized(const std::complex<double>* matrix,
    const std::complex<double>* metric,
    const int n,
    const int wanted,
    double* values,
    std::complex<double>* vectors)
{
    ModuleBase::timer::Scope timer("GpuDenseEigensolver", "solve_generalized");
#ifdef __CUDA
    if (!impl_)
    {
        impl_.reset(new Impl);
    }
    Impl& state = *impl_;
    const std::size_t elements = static_cast<std::size_t>(n) * n;
    const std::size_t matrix_bytes = elements * sizeof(std::complex<double>);
    const std::size_t info_bytes = sizeof(int);
    const int64_t matrix_size = elements;
    state.reserve(state.matrix, matrix_size);
    state.reserve(state.metric, matrix_size);
    state.reserve(state.eigenvalues, n);
    state.reserve(state.info, 1);
    cuDoubleComplex* a = static_cast<cuDoubleComplex*>(state.matrix.data());
    cuDoubleComplex* b = static_cast<cuDoubleComplex*>(state.metric.data());
    double* eigenvalues = state.eigenvalues.data<double>();
    int* device_info = state.info.data<int>();
    const cudaError_t copy_h = cudaMemcpy(a, matrix, matrix_bytes, cudaMemcpyHostToDevice);
    check_cuda(copy_h, "GPU GV matrix upload");
    const cudaError_t copy_g = cudaMemcpy(b, metric, matrix_bytes, cudaMemcpyHostToDevice);
    check_cuda(copy_g, "GPU GV metric upload");
    const bool partial = wanted < n;
    const cusolverEigRange_t range = partial ? CUSOLVER_EIG_RANGE_I : CUSOLVER_EIG_RANGE_ALL;
    int found = 0;
    if (state.query_n != n || state.query_wanted != wanted)
    {
        int lwork = 0;
        cusolverStatus_t status;
        if (!partial)
        {
            status = cusolverDnZhegvd_bufferSize(state.handle,
                CUSOLVER_EIG_TYPE_1,
                CUSOLVER_EIG_MODE_VECTOR,
                CUBLAS_FILL_MODE_UPPER,
                n,
                a,
                n,
                b,
                n,
                eigenvalues,
                &lwork);
        }
        else
        {
            status = cusolverDnZhegvdx_bufferSize(state.handle,
                CUSOLVER_EIG_TYPE_1,
                CUSOLVER_EIG_MODE_VECTOR,
                range,
                CUBLAS_FILL_MODE_UPPER,
                n,
                a,
                n,
                b,
                n,
                0.0,
                0.0,
                1,
                wanted,
                &found,
                eigenvalues,
                &lwork);
        }
        check_solver(status, "GPU GV workspace query");
        if (lwork < 1)
        {
            throw std::runtime_error("Invalid GPU GV workspace query result.");
        }
        state.reserve(state.work, lwork);
        state.lwork = lwork;
        state.query_n = n;
        state.query_wanted = wanted;
    }
    cuDoubleComplex* work = static_cast<cuDoubleComplex*>(state.work.data());
    cusolverStatus_t status;
    if (!partial)
    {
        status = cusolverDnZhegvd(state.handle,
            CUSOLVER_EIG_TYPE_1,
            CUSOLVER_EIG_MODE_VECTOR,
            CUBLAS_FILL_MODE_UPPER,
            n,
            a,
            n,
            b,
            n,
            eigenvalues,
            work,
            state.lwork,
            device_info);
    }
    else
    {
        status = cusolverDnZhegvdx(state.handle,
            CUSOLVER_EIG_TYPE_1,
            CUSOLVER_EIG_MODE_VECTOR,
            range,
            CUBLAS_FILL_MODE_UPPER,
            n,
            a,
            n,
            b,
            n,
            0.0,
            0.0,
            1,
            wanted,
            &found,
            eigenvalues,
            work,
            state.lwork,
            device_info);
    }
    check_solver(status, "GPU GV eigensolve");

    int info = 0;
    const cudaError_t copy_info = cudaMemcpy(&info, device_info, info_bytes, cudaMemcpyDeviceToHost);
    check_cuda(copy_info, "GPU GV status download");
    if (info > n && info <= 2 * static_cast<long long>(n))
    {

        return false;
    }
    if (info != 0 || (partial && found != wanted))
    {

        throw std::runtime_error("GPU GV failed: INFO=" + std::to_string(info) + ", found=" + std::to_string(found));
    }
    const std::size_t returned = static_cast<std::size_t>(n) * wanted;
    const std::size_t returned_value_bytes = static_cast<std::size_t>(wanted) * sizeof(double);
    const std::size_t returned_vector_bytes = returned * sizeof(std::complex<double>);
    const cudaError_t copy_values = cudaMemcpy(values, eigenvalues, returned_value_bytes, cudaMemcpyDeviceToHost);
    check_cuda(copy_values, "GPU GV eigenvalue download");
    const cudaError_t copy_vectors = cudaMemcpy(vectors, a, returned_vector_bytes, cudaMemcpyDeviceToHost);
    check_cuda(copy_vectors, "GPU GV coefficients download");
    return true;
#else
    throw std::runtime_error("GPU generalized RR requires a CUDA build.");
#endif
}

} // namespace dense
} // namespace hsolver
