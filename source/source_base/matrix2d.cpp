#include "source_base/matrix2d.h"
#include "source_base/matrix2d_collective.h"
#include "source_base/matrix2d_reduction.h"
#include "source_base/matrix2d_scalapack.h"
#include "source_base/matrix2d_view.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <memory>
#include <stdexcept>

namespace
{
namespace elementwise
{
template <typename T>
T conj(const T& value)
{
    return value;
}

inline std::complex<float> conj(const std::complex<float>& value)
{
    return std::conj(value);
}

inline std::complex<double> conj(const std::complex<double>& value)
{
    return std::conj(value);
}

template <typename T>
double real(const T& value)
{
    return static_cast<double>(value);
}

inline double real(const std::complex<float>& value)
{
    return static_cast<double>(value.real());
}

inline double real(const std::complex<double>& value)
{
    return value.real();
}

template <typename T>
double imag(const T&)
{
    return 0.0;
}

inline double imag(const std::complex<float>& value)
{
    return static_cast<double>(value.imag());
}

inline double imag(const std::complex<double>& value)
{
    return value.imag();
}
} // namespace elementwise

} // namespace

namespace matrix2d_impl
{
namespace scalapack
{
template <typename T>
void transpose(Matrix2dView<T> result, Matrix2dView<const T> source)
{
    static_cast<void>(result);
    static_cast<void>(source);
    throw std::invalid_argument(
        "ScaLAPACK transpose supports only double and complex<double>");
}

template <typename T>
void conjugate_transpose(Matrix2dView<T> result,
                         Matrix2dView<const T> source)
{
    static_cast<void>(result);
    static_cast<void>(source);
    throw std::invalid_argument(
        "ScaLAPACK conjugate transpose supports only double and complex<double>");
}

template <typename T>
void matmul(Matrix2dView<T> result,
            Matrix2dView<const T> left,
            Matrix2dView<const T> right)
{
    static_cast<void>(result);
    static_cast<void>(left);
    static_cast<void>(right);
    throw std::invalid_argument(
        "ScaLAPACK matmul supports only double and complex<double>");
}
} // namespace scalapack
} // namespace matrix2d_impl

template <typename T>
void NDArray<T, 2>::assert_compatible_layout(
    const NDArray<T, 2>& other) const
{
    if (nrglb() != other.nrglb()
        || ncglb() != other.ncglb())
    {
        throw std::invalid_argument("matrix global shapes are incompatible");
    }

    const bool distributed = has_distribution();
    const bool distributed_ = other.has_distribution();
    if (distributed != distributed_)
    {
        throw std::invalid_argument("host and distributed matrices are incompatible");
    }

    if (!distributed)
    {
        return;
    }

    if (distrib_->nb() != other.distribution()->nb())
    {
        throw std::invalid_argument("matrix block sizes are incompatible");
    }
    if (distrib_->nrproc() != other.distribution()->nrproc()
        || distrib_->ncproc() != other.distribution()->ncproc())
    {
        throw std::invalid_argument("matrix process grids are incompatible");
    }
    if (distrib_ != other.distribution())
    {
        throw std::invalid_argument("matrix distribution handles are incompatible");
    }
}

template <typename T>
Matrix2dView<T> NDArray<T, 2>::view()
{
    if (distrib_ != nullptr)
    {
        return Matrix2dView<T>(data(), distrib_);
    }
    return Matrix2dView<T>(data(), nrglb(), ncglb());
}

template <typename T>
Matrix2dView<const T> NDArray<T, 2>::view() const
{
    if (distrib_ != nullptr)
    {
        return Matrix2dView<const T>(data(), distrib_);
    }
    return Matrix2dView<const T>(data(), nrglb(), ncglb());
}

template <typename T>
void NDArray<T, 2>::fill(const T& value)
{
    if (distrib_ == nullptr)
    {
        std::fill(storage_.begin(), storage_.end(), value);
        return;
    }
    const int nrloc = this->nrloc();
    const int ncloc = this->ncloc();
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
    for (int irloc = 0; irloc < nrloc; ++irloc)
    {
        for (int icloc = 0; icloc < ncloc; ++icloc)
        {
            atloc(irloc, icloc) = value;
        }
    }
}

template <typename T>
void NDArray<T, 2>::scale(const T& value)
{
    if (distrib_ == nullptr)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (std::size_t i = 0; i < storage_.size(); ++i)
        {
            storage_[i] *= value;
        }
        return;
    }
    const int nrloc = this->nrloc();
    const int ncloc = this->ncloc();
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
    for (int irloc = 0; irloc < nrloc; ++irloc)
    {
        for (int icloc = 0; icloc < ncloc; ++icloc)
        {
            atloc(irloc, icloc) *= value;
        }
    }
}

template <typename T>
NDArray<T, 2>& NDArray<T, 2>::operator+=(const NDArray<T, 2>& other)
{
    assert_compatible_layout(other);
    const int nrloc = this->nrloc();
    const int ncloc = this->ncloc();
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
    for (int irloc = 0; irloc < nrloc; ++irloc)
    {
        for (int icloc = 0; icloc < ncloc; ++icloc)
        {
            atloc(irloc, icloc) += other.atloc(irloc, icloc);
        }
    }
    return *this;
}

template <typename T>
NDArray<T, 2>& NDArray<T, 2>::operator-=(const NDArray<T, 2>& other)
{
    assert_compatible_layout(other);
    const int nrloc = this->nrloc();
    const int ncloc = this->ncloc();
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
    for (int irloc = 0; irloc < nrloc; ++irloc)
    {
        for (int icloc = 0; icloc < ncloc; ++icloc)
        {
            atloc(irloc, icloc) -= other.atloc(irloc, icloc);
        }
    }
    return *this;
}

template <typename T>
NDArray<T, 2>& NDArray<T, 2>::operator*=(const NDArray<T, 2>& other)
{
    assert_compatible_layout(other);
    const int nrloc = this->nrloc();
    const int ncloc = this->ncloc();
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
    for (int irloc = 0; irloc < nrloc; ++irloc)
    {
        for (int icloc = 0; icloc < ncloc; ++icloc)
        {
            atloc(irloc, icloc) *= other.atloc(irloc, icloc);
        }
    }
    return *this;
}

template <typename T>
NDArray<T, 2>& NDArray<T, 2>::conjugate()
{
    const int nrloc = this->nrloc();
    const int ncloc = this->ncloc();
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
    for (int i = 0; i < nrloc; ++i)
    {
        for (int j = 0; j < ncloc; ++j)
        {
            atloc(i, j) = elementwise::conj(atloc(i, j));
        }
    }
    return *this;
}

template <typename T>
NDArray<T, 2> NDArray<T, 2>::transpose() const
{
    if (!is_active())
    {
        // Serial degenerate path: a serial frozen distribution owns the full
        // matrix on one rank with column-major local layout; treat it like the
        // host path so distributed operators degrade to the host side.
        Matrix2d<T> result(ncloc_, nrloc_);
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
        for (int i = 0; i < nrloc_; ++i)
        {
            for (int j = 0; j < ncloc_; ++j)
            {
                result.atloc(j, i) = atloc(i, j);
            }
        }
        return result;
    }

#ifdef __MPI
    const int nr = ncglb();
    const int nc = nrglb();
    const int nb = distrib_->nb();
    const std::shared_ptr<const Parallel_Blkcyclic> result_distrib
        = make_borrowed_distribution(nr, nc, nb, distrib_);
    Matrix2d<T> result(result_distrib);
    matrix2d_impl::scalapack::transpose(result.view(), view());
    return result;
#else
    throw std::invalid_argument("distributed transpose requires MPI support");
#endif
}

template <typename T>
NDArray<T, 2> NDArray<T, 2>::dagger() const
{
    if (!is_active())
    {
        Matrix2d<T> result(ncloc_, nrloc_);
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
        for (int i = 0; i < nrloc_; ++i)
        {
            for (int j = 0; j < ncloc_; ++j)
            {
                result.atloc(j, i) = elementwise::conj(atloc(i, j));
            }
        }
        return result;
    }

#ifdef __MPI
    const int nr = ncglb();
    const int nc = nrglb();
    const int nb = distrib_->nb();
    const std::shared_ptr<const Parallel_Blkcyclic> result_distrib
        = make_borrowed_distribution(nr, nc, nb, distrib_);
    Matrix2d<T> result(result_distrib);
    matrix2d_impl::scalapack::conjugate_transpose(result.view(), view());
    return result;
#else
    throw std::invalid_argument(
        "distributed conjugate transpose requires MPI support");
#endif
}

template <typename T>
NDArray<T, 2> NDArray<T, 2>::matmul(const NDArray<T, 2>& other) const
{
    const int m = nrglb();
    const int k = ncglb();
    const int n = other.ncglb();
    {
        const int k_ = other.nrglb();
        if (k != k_)
        {
            throw std::invalid_argument("matrix inner dimensions are incompatible");
        }
    }

    const bool serial = !is_active();
    const bool serial_ = !other.is_active();
    if (serial != serial_)
    {
        throw std::invalid_argument("host and distributed matrices are incompatible");
    }

    if (serial)
    {
        // Serial degenerate path: both sides are host or serial frozen
        // distributions, so matmul degrades to the host triple loop.
        Matrix2d<T> result(m, n);
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
        for (int i = 0; i < m; ++i)
        {
            for (int j = 0; j < n; ++j)
            {
                T value = T();
                for (int k_ = 0; k_ < k; ++k_)
                {
                    value += atloc(i, k_) * other.atloc(k_, j);
                }
                result.atloc(i, j) = value;
            }
        }
        return result;
    }

#ifdef __MPI
    const std::shared_ptr<const Parallel_Blkcyclic> other_distrib
        = other.distribution();
    if (distrib_->nb() != other_distrib->nb())
    {
        throw std::invalid_argument("matrix block sizes are incompatible");
    }
    if (distrib_->nrproc() != other_distrib->nrproc()
        || distrib_->ncproc() != other_distrib->ncproc())
    {
        throw std::invalid_argument("matrix process grids are incompatible");
    }
    if (distrib_->ictxt() != other_distrib->ictxt())
    {
        throw std::invalid_argument("matrix BLACS contexts are incompatible");
    }

    const int nb = distrib_->nb();
    const std::shared_ptr<const Parallel_Blkcyclic> result_distrib
        = make_borrowed_distribution(m, n, nb, distrib_);
    Matrix2d<T> result(result_distrib);
    matrix2d_impl::scalapack::matmul(result.view(), view(), other.view());
    return result;
#else
    throw std::invalid_argument("distributed matmul requires MPI support");
#endif
}

template <typename T>
NDArray<double, 2> NDArray<T, 2>::real() const
{
    Matrix2d<double> result = has_distribution()
        ? Matrix2d<double>(distribution())
        : Matrix2d<double>(nrglb(), ncglb());
    const int nrloc = this->nrloc();
    const int ncloc = this->ncloc();
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
    for (int i = 0; i < nrloc; ++i)
    {
        for (int j = 0; j < ncloc; ++j)
        {
            result.atloc(i, j) = elementwise::real(atloc(i, j));
        }
    }
    return result;
}

template <typename T>
NDArray<double, 2> NDArray<T, 2>::imag() const
{
    Matrix2d<double> result = has_distribution()
        ? Matrix2d<double>(distribution())
        : Matrix2d<double>(nrglb(), ncglb());
    const int nrloc = this->nrloc();
    const int ncloc = this->ncloc();
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
    for (int i = 0; i < nrloc; ++i)
    {
        for (int j = 0; j < ncloc; ++j)
        {
            result.atloc(i, j) = elementwise::imag(atloc(i, j));
        }
    }
    return result;
}

template <typename T>
T NDArray<T, 2>::sum() const
{
    return matrix2d_impl::sum(*this);
}

template <typename T>
T NDArray<T, 2>::dot(const NDArray<T, 2>& other) const
{
    assert_compatible_layout(other);
    return matrix2d_impl::dot(*this, other);
}

template <typename T>
double NDArray<T, 2>::norm() const
{
    return matrix2d_impl::norm(*this);
}

template <typename T>
double NDArray<T, 2>::max_abs() const
{
    return matrix2d_impl::max_abs(*this);
}

template <typename T>
std::vector<T> NDArray<T, 2>::diagonal() const
{
    return matrix2d_impl::diagonal(*this);
}

template <typename T>
NDArray<T, 2> NDArray<T, 2>::gather() const
{
    if (distrib_ == nullptr)
    {
        return *this;
    }
    if (distrib_->is_serial())
    {
        Matrix2d<T> result(nrglb(), ncglb());
        for (int irloc = 0; irloc < nrloc(); ++irloc)
        {
            const int irglb = distrib_->irglb(irloc);
            for (int icloc = 0; icloc < ncloc(); ++icloc)
            {
                const int icglb = distrib_->icglb(icloc);
                result.atloc(irglb, icglb) = atloc(irloc, icloc);
            }
        }
        return result;
    }
#ifdef __MPI
    const MPI_Comm comm = distrib_->comm();
    if (comm == MPI_COMM_NULL)
    {
        throw std::invalid_argument("gather requires a valid MPI comm");
    }
    return matrix2d_impl::allgather(*this, comm);
#else
    return *this;
#endif
}

template <typename T>
NDArray<T, 2> NDArray<T, 2>::gather_on_root(const int iproot) const
{
    if (distrib_ == nullptr)
    {
        if (iproot == 0)
        {
            return *this;
        }
        return NDArray(EmptyMatrixTag());
    }
    if (distrib_->is_serial())
    {
        if (iproot != 0)
        {
            return NDArray(EmptyMatrixTag());
        }
        Matrix2d<T> result(nrglb(), ncglb());
        for (int irloc = 0; irloc < nrloc(); ++irloc)
        {
            const int irglb = distrib_->irglb(irloc);
            for (int icloc = 0; icloc < ncloc(); ++icloc)
            {
                const int icglb = distrib_->icglb(icloc);
                result.atloc(irglb, icglb) = atloc(irloc, icloc);
            }
        }
        return result;
    }
#ifdef __MPI
    const MPI_Comm comm = distrib_->comm();
    if (comm == MPI_COMM_NULL)
    {
        throw std::invalid_argument("gather_on_root requires a valid MPI comm");
    }
    std::unique_ptr<Matrix2d<T>> gathered
        = matrix2d_impl::gather_on_root(*this, iproot, comm);
    if (gathered)
    {
        return *gathered;
    }
    return NDArray(EmptyMatrixTag());
#else
    static_cast<void>(iproot);
    return NDArray(EmptyMatrixTag());
#endif
}

template <typename T>
void NDArray<T, 2>::scatter(const NDArray<T, 2>& matglb)
{
    if (distrib_ == nullptr)
    {
        if (matglb.has_distribution())
        {
            throw std::invalid_argument("host scatter requires a host global matrix");
        }
        if (nrloc_ != matglb.nrglb()
            || ncloc_ != matglb.ncglb())
        {
            throw std::invalid_argument("scatter global matrix shape does not match matrix");
        }
        *this = matglb;
        return;
    }

    if (distrib_->is_serial())
    {
        if (matglb.has_distribution())
        {
            throw std::invalid_argument("scatter requires a host global matrix");
        }
        if (nrglb() != matglb.nrglb()
            || ncglb() != matglb.ncglb())
        {
            throw std::invalid_argument("scatter global matrix shape does not match distribution");
        }
        for (int irloc = 0; irloc < nrloc(); ++irloc)
        {
            for (int icloc = 0; icloc < ncloc(); ++icloc)
            {
                atloc(irloc, icloc)
                    = matglb.atloc(distrib_->irglb(irloc),
                                   distrib_->icglb(icloc));
            }
        }
        return;
    }

#ifdef __MPI
    matrix2d_impl::scatter_from_host(matglb, *this);
#else
    throw std::invalid_argument("active distributed scatter requires MPI support");
#endif
}

template class NDArray<int, 2>;
template class NDArray<float, 2>;
template class NDArray<double, 2>;
template class NDArray<std::complex<float>, 2>;
template class NDArray<std::complex<double>, 2>;
