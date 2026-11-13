#ifndef ABACUS_SOURCE_SOURCE_BASE_MATRIX2D_H
#define ABACUS_SOURCE_SOURCE_BASE_MATRIX2D_H

#include "source_base/parallel_blkcyc.h"
#include "source_base/ndarray.h"

#include <algorithm>
#include <cassert>
#include <complex>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

template <typename T>
class Matrix2dView;

template <typename T>
class NDArray<T, 2> final
{
  public:
    NDArray() = delete;

    NDArray(const int nrloc, const int ncloc)
        : nrloc_(nrloc),
          ncloc_(ncloc)
    {
        if (nrloc_ <= 0)
        {
            throw std::invalid_argument("nrloc must be positive");
        }
        if (ncloc_ <= 0)
        {
            throw std::invalid_argument("ncloc must be positive");
        }
        const std::size_t n
            = static_cast<std::size_t>(nrloc_) * static_cast<std::size_t>(ncloc_);
        storage_.resize(n);
    }

#ifdef __MPI
    NDArray(const int nrglb,
            const int ncglb,
            const MPI_Comm comm)
    {
        const std::shared_ptr<const Parallel_Blkcyclic> dist
            = make_frozen_distribution(nrglb, ncglb, 1, comm);
        initialize_distributed_storage(dist);
    }

    NDArray(const int nrglb,
            const int ncglb,
            const int nb,
            const MPI_Comm comm)
    {
        const std::shared_ptr<const Parallel_Blkcyclic> dist
            = make_frozen_distribution(nrglb, ncglb, nb, comm);
        initialize_distributed_storage(dist);
    }

    NDArray(const int nrglb,
            const int ncglb,
            const int nrproc,
            const int ncproc,
            const int nb,
            const MPI_Comm comm)
    {
        const std::shared_ptr<const Parallel_Blkcyclic> dist
            = make_frozen_distribution(nrglb, ncglb, nrproc, ncproc, nb, comm);
        initialize_distributed_storage(dist);
    }
#endif

    explicit NDArray(std::shared_ptr<const Parallel_Blkcyclic> dist)
        : dist_(std::move(dist))
    {
        if (dist_ == nullptr)
        {
            throw std::invalid_argument("distribution must not be null");
        }
        initialize_distributed_storage();
    }

    NDArray(const NDArray&) = default;
    NDArray(NDArray&&) noexcept = default;
    ~NDArray() = default;

    NDArray& operator=(const NDArray&) = default;
    NDArray& operator=(NDArray&&) noexcept = default;

    void fill(const T& value);
    void scale(const T& value);
    NDArray& operator+=(const NDArray& other);
    NDArray& operator-=(const NDArray& other);
    NDArray& operator*=(const NDArray& other);

    NDArray& conjugate();
    NDArray transpose() const;
    NDArray dagger() const;
    NDArray matmul(const NDArray& other) const;

    NDArray<double, 2> real() const;
    NDArray<double, 2> imag() const;
    void real(NDArray<double, 2>& output) const;
    void imag(NDArray<double, 2>& output) const;

    T sum() const;
    double sum_squares() const;
    T dot(const NDArray& other) const;
    double norm2() const;
    double max_abs() const;
    std::vector<T> diagonal() const;

    NDArray gather() const;
    NDArray gather_on_root(const int iproot) const;
    void scatter(const NDArray& matglb);

    T& at_global(const int i, const int j)
    {
        assert(i >= 0);
        assert(i < nrglb());
        assert(j >= 0);
        assert(j < ncglb());

        if (dist_ == nullptr)
        {
            return at_local(i, j);
        }

        const std::size_t iglb = static_cast<std::size_t>(i);
        const std::size_t jglb = static_cast<std::size_t>(j);
        const bool owns = dist_->in_this_processor(iglb, jglb);
        assert(owns);

        const int irloc = dist_->irloc(i);
        const int icloc = dist_->icloc(j);
        assert(irloc != -1);
        assert(icloc != -1);
        return at_local(irloc, icloc);
    }

    const T& at_global(const int i, const int j) const
    {
        assert(i >= 0);
        assert(i < nrglb());
        assert(j >= 0);
        assert(j < ncglb());

        if (dist_ == nullptr)
        {
            return at_local(i, j);
        }

        const std::size_t iglb = static_cast<std::size_t>(i);
        const std::size_t jglb = static_cast<std::size_t>(j);
        const bool owns = dist_->in_this_processor(iglb, jglb);
        assert(owns);

        const int irloc = dist_->irloc(i);
        const int icloc = dist_->icloc(j);
        assert(irloc != -1);
        assert(icloc != -1);
        return at_local(irloc, icloc);
    }

    T& at_local(const int i, const int j)
    {
        assert(i >= 0);
        assert(i < nrloc());
        assert(j >= 0);
        assert(j < ncloc());

        std::size_t ioff = 0;
        if (dist_ == nullptr)
        {
            ioff = static_cast<std::size_t>(i) * static_cast<std::size_t>(ncloc_)
                + static_cast<std::size_t>(j);
        }
        else
        {
            const int lld = dist_->lld();
            ioff = static_cast<std::size_t>(j) * static_cast<std::size_t>(lld)
                + static_cast<std::size_t>(i);
        }
        assert(ioff < storage_.size());
        return storage_[ioff];
    }

    const T& at_local(const int i, const int j) const
    {
        assert(i >= 0);
        assert(i < nrloc());
        assert(j >= 0);
        assert(j < ncloc());

        std::size_t ioff = 0;
        if (dist_ == nullptr)
        {
            ioff = static_cast<std::size_t>(i) * static_cast<std::size_t>(ncloc_)
                + static_cast<std::size_t>(j);
        }
        else
        {
            const int lld = dist_->lld();
            ioff = static_cast<std::size_t>(j) * static_cast<std::size_t>(lld)
                + static_cast<std::size_t>(i);
        }
        assert(ioff < storage_.size());
        return storage_[ioff];
    }

    bool owns_global_index(const int i, const int j) const
    {
        assert(i >= 0);
        assert(i < nrglb());
        assert(j >= 0);
        assert(j < ncglb());
        if (dist_ == nullptr)
        {
            return true;
        }
        const std::size_t iglb = static_cast<std::size_t>(i);
        const std::size_t jglb = static_cast<std::size_t>(j);
        return dist_->in_this_processor(iglb, jglb);
    }

    T* data() { return storage_.data(); }
    const T* data() const { return storage_.data(); }

    int nrglb() const { return dist_ == nullptr ? nrloc_ : dist_->nrglb(); }
    int ncglb() const { return dist_ == nullptr ? ncloc_ : dist_->ncglb(); }
    int nrloc() const { return dist_ == nullptr ? nrloc_ : dist_->nrloc(); }
    int ncloc() const { return dist_ == nullptr ? ncloc_ : dist_->ncloc(); }
    std::size_t size() const
    {
        return dist_ == nullptr
            ? static_cast<std::size_t>(nrloc_) * static_cast<std::size_t>(ncloc_)
            : dist_->size();
    }
    int nrproc() const { return dist_ == nullptr ? 1 : dist_->nrproc(); }
    int ncproc() const { return dist_ == nullptr ? 1 : dist_->ncproc(); }
    int irproc() const { return dist_ == nullptr ? 0 : dist_->irproc(); }
    int icproc() const { return dist_ == nullptr ? 0 : dist_->icproc(); }
    int lld() const { return dist_ == nullptr ? ncloc_ : dist_->lld(); }

    const int* desc() const
    {
#ifdef __MPI
        return dist_ != nullptr ? dist_->desc() : nullptr;
#else
        return nullptr;
#endif
    }

    int ictxt() const
    {
#ifdef __MPI
        return dist_ != nullptr ? dist_->ictxt() : -1;
#else
        return -1;
#endif
    }

    const std::shared_ptr<const Parallel_Blkcyclic>& distribution() const { return dist_; }
    bool is_distributed() const { return dist_ != nullptr; }

    Matrix2dView<T> view();
    Matrix2dView<const T> view() const;

  private:
    struct EmptyMatrixTag
    {
    };

    explicit NDArray(EmptyMatrixTag)
    {
    }

    template <typename U>
    void validate_compatible_layout(const NDArray<U, 2>& other) const
    {
        if (nrglb() != other.nrglb()
            || ncglb() != other.ncglb())
        {
            throw std::invalid_argument("matrix global shapes are incompatible");
        }

        const bool this_is_distributed = dist_ != nullptr;
        const bool other_is_distributed = other.distribution() != nullptr;
        if (this_is_distributed != other_is_distributed)
        {
            throw std::invalid_argument("host and distributed matrices are incompatible");
        }

        if (!this_is_distributed)
        {
            return;
        }

        if (dist_->nb() != other.distribution()->nb())
        {
            throw std::invalid_argument("matrix block sizes are incompatible");
        }
        if (dist_->nrproc() != other.distribution()->nrproc()
            || dist_->ncproc() != other.distribution()->ncproc())
        {
            throw std::invalid_argument("matrix process grids are incompatible");
        }
        if (dist_ != other.distribution())
        {
            throw std::invalid_argument("matrix distribution handles are incompatible");
        }
    }

    void initialize_distributed_storage(
        const std::shared_ptr<const Parallel_Blkcyclic>& dist)
    {
        if (dist == nullptr)
        {
            throw std::invalid_argument("distribution must not be null");
        }
        dist_ = dist;
        initialize_distributed_storage();
    }

    void initialize_distributed_storage()
    {
        nrloc_ = dist_->nrloc();
        ncloc_ = dist_->ncloc();
        assert(nrloc_ >= 0);
        assert(ncloc_ >= 0);

        const int lld = std::max(nrloc_, 1);
        const int ncmax = std::max(ncloc_, lld);
        const std::size_t ncm = static_cast<std::size_t>(lld) * static_cast<std::size_t>(ncloc_);
        const std::size_t nrm = static_cast<std::size_t>(nrloc_) * static_cast<std::size_t>(ncmax);
        const std::size_t n = std::max(ncm, nrm);
        storage_.resize(n);
    }

    int nrloc_ = 0;
    int ncloc_ = 0;
    std::vector<T> storage_;
    std::shared_ptr<const Parallel_Blkcyclic> dist_;
};

template <typename T>
using Matrix2d = NDArray<T, 2>;

static_assert(std::is_same<Matrix2d<int>, NDArray<int, 2>>::value,
              "Matrix2d<T> and NDArray<T, 2> must be the same type.");
static_assert(!std::is_same<NDArray<int>, NDArray<int, 2>>::value,
              "NDArray<T> and NDArray<T, 2> must be distinct types.");
static_assert(!std::is_base_of<NDArray<int>, NDArray<int, 2>>::value,
              "The rank-two specialization must compose storage rather than inherit.");

#endif
