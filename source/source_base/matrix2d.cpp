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
template <typename T>
T conjugate_element(const T& value)
{
    return value;
}

template <>
std::complex<float> conjugate_element<std::complex<float>>(
    const std::complex<float>& value)
{
    return std::conj(value);
}

template <>
std::complex<double> conjugate_element<std::complex<double>>(
    const std::complex<double>& value)
{
    return std::conj(value);
}

template <typename T>
double real_element(const T& value)
{
    return static_cast<double>(value);
}

template <>
double real_element<std::complex<float>>(const std::complex<float>& value)
{
    return static_cast<double>(value.real());
}

template <>
double real_element<std::complex<double>>(const std::complex<double>& value)
{
    return value.real();
}

template <typename T>
double imaginary_element(const T&)
{
    return 0.0;
}

template <>
double imaginary_element<std::complex<float>>(const std::complex<float>& value)
{
    return static_cast<double>(value.imag());
}

template <>
double imaginary_element<std::complex<double>>(const std::complex<double>& value)
{
    return value.imag();
}

template <typename T>
Matrix2d<double> extract_real_or_imaginary(const Matrix2d<T>& matrix,
                                           const bool extract_real)
{
    Matrix2d<double> result = matrix.is_distributed()
        ? Matrix2d<double>(matrix.distribution())
        : Matrix2d<double>(matrix.nrglb(),
                           matrix.ncglb());

    const int nrloc = matrix.nrloc();
    const int ncloc = matrix.ncloc();
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
    for (int i = 0; i < nrloc; ++i)
    {
        for (int j = 0; j < ncloc; ++j)
        {
            const T& value = matrix.at_local(i, j);
            result.at_local(i, j)
                = extract_real ? real_element(value) : imaginary_element(value);
        }
    }
    return result;
}

template <typename T>
void write_real_or_imaginary(const Matrix2d<T>& matrix,
                             Matrix2d<double>& output,
                             const bool extract_real)
{
    const int nrloc = matrix.nrloc();
    const int ncloc = matrix.ncloc();
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
    for (int i = 0; i < nrloc; ++i)
    {
        for (int j = 0; j < ncloc; ++j)
        {
            const T& value = matrix.at_local(i, j);
            output.at_local(i, j)
                = extract_real ? real_element(value) : imaginary_element(value);
        }
    }
}

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
Matrix2dView<T> NDArray<T, 2>::view()
{
#ifdef __MPI
    if (dist_ != nullptr)
    {
        return Matrix2dView<T>(data(), dist_);
    }
#endif
    return Matrix2dView<T>(data(), nrglb(), ncglb());
}

template <typename T>
Matrix2dView<const T> NDArray<T, 2>::view() const
{
#ifdef __MPI
    if (dist_ != nullptr)
    {
        return Matrix2dView<const T>(data(), dist_);
    }
#endif
    return Matrix2dView<const T>(data(), nrglb(), ncglb());
}

template <typename T>
void NDArray<T, 2>::fill(const T& value)
{
    if (dist_ == nullptr)
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
            at_local(irloc, icloc) = value;
        }
    }
}

template <typename T>
void NDArray<T, 2>::scale(const T& value)
{
    if (dist_ == nullptr)
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
            at_local(irloc, icloc) *= value;
        }
    }
}

template <typename T>
NDArray<T, 2>& NDArray<T, 2>::operator+=(const NDArray<T, 2>& other)
{
    validate_compatible_layout(other);
    const int nrloc = this->nrloc();
    const int ncloc = this->ncloc();
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
    for (int irloc = 0; irloc < nrloc; ++irloc)
    {
        for (int icloc = 0; icloc < ncloc; ++icloc)
        {
            at_local(irloc, icloc) += other.at_local(irloc, icloc);
        }
    }
    return *this;
}

template <typename T>
NDArray<T, 2>& NDArray<T, 2>::operator-=(const NDArray<T, 2>& other)
{
    validate_compatible_layout(other);
    const int nrloc = this->nrloc();
    const int ncloc = this->ncloc();
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
    for (int irloc = 0; irloc < nrloc; ++irloc)
    {
        for (int icloc = 0; icloc < ncloc; ++icloc)
        {
            at_local(irloc, icloc) -= other.at_local(irloc, icloc);
        }
    }
    return *this;
}

template <typename T>
NDArray<T, 2>& NDArray<T, 2>::operator*=(const NDArray<T, 2>& other)
{
    validate_compatible_layout(other);
    const int nrloc = this->nrloc();
    const int ncloc = this->ncloc();
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
    for (int irloc = 0; irloc < nrloc; ++irloc)
    {
        for (int icloc = 0; icloc < ncloc; ++icloc)
        {
            at_local(irloc, icloc) *= other.at_local(irloc, icloc);
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
            at_local(i, j) = conjugate_element(at_local(i, j));
        }
    }
    return *this;
}

template <typename T>
NDArray<T, 2> NDArray<T, 2>::transpose() const
{
    if (dist_ == nullptr || dist_->is_serial())
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
                result.at_local(j, i) = at_local(i, j);
            }
        }
        return result;
    }

#ifdef __MPI
    const int result_rows = ncglb();
    const int result_cols = nrglb();
    const int nb = dist_->nb();
    const std::shared_ptr<const Parallel_Blkcyclic> result_dist
        = make_borrowed_distribution(result_rows, result_cols, nb, dist_);
    Matrix2d<T> result(result_dist);
    matrix2d_impl::scalapack::transpose(result.view(), view());
    return result;
#else
    throw std::invalid_argument("distributed transpose requires MPI support");
#endif
}

template <typename T>
NDArray<T, 2> NDArray<T, 2>::dagger() const
{
    if (dist_ == nullptr || dist_->is_serial())
    {
        Matrix2d<T> result(ncloc_, nrloc_);
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
        for (int i = 0; i < nrloc_; ++i)
        {
            for (int j = 0; j < ncloc_; ++j)
            {
                result.at_local(j, i) = conjugate_element(at_local(i, j));
            }
        }
        return result;
    }

#ifdef __MPI
    const int result_rows = ncglb();
    const int result_cols = nrglb();
    const int nb = dist_->nb();
    const std::shared_ptr<const Parallel_Blkcyclic> result_dist
        = make_borrowed_distribution(result_rows, result_cols, nb, dist_);
    Matrix2d<T> result(result_dist);
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
    const int left_rows = nrglb();
    const int shared_dimension = ncglb();
    const int other_shared_dimension = other.nrglb();
    const int result_cols = other.ncglb();
    if (shared_dimension != other_shared_dimension)
    {
        throw std::invalid_argument("matrix inner dimensions are incompatible");
    }

    const bool this_is_serial
        = dist_ == nullptr || dist_->is_serial();
    const bool other_is_serial
        = other.distribution() == nullptr || other.distribution()->is_serial();
    if (this_is_serial != other_is_serial)
    {
        throw std::invalid_argument("host and distributed matrices are incompatible");
    }

    if (this_is_serial)
    {
        // Serial degenerate path: both sides are host or serial frozen
        // distributions, so matmul degrades to the host triple loop.
        Matrix2d<T> result(left_rows, result_cols);
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
        for (int i = 0; i < left_rows; ++i)
        {
            for (int j = 0; j < result_cols; ++j)
            {
                T value = T();
                for (int k = 0; k < shared_dimension; ++k)
                {
                    value += at_local(i, k) * other.at_local(k, j);
                }
                result.at_local(i, j) = value;
            }
        }
        return result;
    }

#ifdef __MPI
    const std::shared_ptr<const Parallel_Blkcyclic> other_dist
        = other.distribution();
    if (dist_->nb() != other_dist->nb())
    {
        throw std::invalid_argument("matrix block sizes are incompatible");
    }
    if (dist_->nrproc() != other_dist->nrproc()
        || dist_->ncproc() != other_dist->ncproc())
    {
        throw std::invalid_argument("matrix process grids are incompatible");
    }
    if (dist_->ictxt() != other_dist->ictxt())
    {
        throw std::invalid_argument("matrix BLACS contexts are incompatible");
    }

    const int nb = dist_->nb();
    const std::shared_ptr<const Parallel_Blkcyclic> result_dist
        = make_borrowed_distribution(left_rows, result_cols, nb, dist_);
    Matrix2d<T> result(result_dist);
    matrix2d_impl::scalapack::matmul(result.view(), view(), other.view());
    return result;
#else
    throw std::invalid_argument("distributed matmul requires MPI support");
#endif
}

template <typename T>
NDArray<double, 2> NDArray<T, 2>::real() const
{
    return extract_real_or_imaginary(*this, true);
}

template <typename T>
NDArray<double, 2> NDArray<T, 2>::imag() const
{
    return extract_real_or_imaginary(*this, false);
}

template <typename T>
void NDArray<T, 2>::real(NDArray<double, 2>& output) const
{
    validate_compatible_layout(output);
    write_real_or_imaginary(*this, output, true);
}

template <typename T>
void NDArray<T, 2>::imag(NDArray<double, 2>& output) const
{
    validate_compatible_layout(output);
    write_real_or_imaginary(*this, output, false);
}

template <typename T>
T NDArray<T, 2>::sum() const
{
    return matrix2d_impl::sum(*this);
}

template <typename T>
double NDArray<T, 2>::sum_squares() const
{
    return matrix2d_impl::sum_squares(*this);
}

template <typename T>
T NDArray<T, 2>::dot(const NDArray<T, 2>& other) const
{
    validate_compatible_layout(other);
    return matrix2d_impl::dot(*this, other);
}

template <typename T>
double NDArray<T, 2>::norm2() const
{
    return matrix2d_impl::norm2(*this);
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
    if (dist_ == nullptr)
    {
        return *this;
    }
    if (dist_->is_serial())
    {
        Matrix2d<T> result(nrglb(), ncglb());
        for (int irloc = 0; irloc < nrloc(); ++irloc)
        {
            const int irglb = dist_->irglb(irloc);
            for (int icloc = 0; icloc < ncloc(); ++icloc)
            {
                const int icglb = dist_->icglb(icloc);
                result.at_local(irglb, icglb) = at_local(irloc, icloc);
            }
        }
        return result;
    }
#ifdef __MPI
    const MPI_Comm comm = dist_->comm();
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
    if (dist_ == nullptr)
    {
        if (iproot == 0)
        {
            return *this;
        }
        return NDArray(EmptyMatrixTag());
    }
    if (dist_->is_serial())
    {
        if (iproot != 0)
        {
            return NDArray(EmptyMatrixTag());
        }
        Matrix2d<T> result(nrglb(), ncglb());
        for (int irloc = 0; irloc < nrloc(); ++irloc)
        {
            const int irglb = dist_->irglb(irloc);
            for (int icloc = 0; icloc < ncloc(); ++icloc)
            {
                const int icglb = dist_->icglb(icloc);
                result.at_local(irglb, icglb) = at_local(irloc, icloc);
            }
        }
        return result;
    }
#ifdef __MPI
    const MPI_Comm comm = dist_->comm();
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
    if (dist_ == nullptr)
    {
        if (matglb.is_distributed())
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

    if (dist_->is_serial())
    {
        if (matglb.is_distributed())
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
                at_local(irloc, icloc)
                    = matglb.at_local(dist_->irglb(irloc),
                                             dist_->icglb(icloc));
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
