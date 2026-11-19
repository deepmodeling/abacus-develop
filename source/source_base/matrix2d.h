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
        const std::shared_ptr<const Parallel_Blkcyclic> distrib // distribution, not distance
            = make_frozen_distribution(nrglb, ncglb, 1, comm);
        initialize_distributed_storage(distrib);
    }

    NDArray(const int nrglb,
            const int ncglb,
            const int nb,
            const MPI_Comm comm)
    {
        const std::shared_ptr<const Parallel_Blkcyclic> distrib
            = make_frozen_distribution(nrglb, ncglb, nb, comm);
        initialize_distributed_storage(distrib);
    }

    NDArray(const int nrglb,
            const int ncglb,
            const int nrproc,
            const int ncproc,
            const int nb,
            const MPI_Comm comm)
    {
        const std::shared_ptr<const Parallel_Blkcyclic> distrib
            = make_frozen_distribution(nrglb, ncglb, nrproc, ncproc, nb, comm);
        initialize_distributed_storage(distrib);
    }
#endif

    explicit NDArray(std::shared_ptr<const Parallel_Blkcyclic> distrib)
        : distrib_(std::move(distrib))
    {
        if (distrib_ == nullptr)
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

    /** Return a newly allocated matrix containing the real parts. */
    NDArray<double, 2> real() const;
    /** Return a newly allocated matrix containing the imaginary parts. */
    NDArray<double, 2> imag() const;

    T sum() const;
    T dot(const NDArray& other) const;
    /** Return the Frobenius norm sqrt(sum_ij |a_ij|^2). */
    double norm() const;
    double max_abs() const;
    std::vector<T> diagonal() const;

    NDArray gather() const;

    /**
     * @brief Gather the distributed matrix on one rank.
     *
     * MPI collectives require every rank to call this function. When
     * iproot == 0, ranks other than rank zero return a 0 x 0 host placeholder.
     * Such placeholder results, or a root-only gathered matrix, must not be fed
     * back to scatter() because scatter() requires a full host matrix on every
     * rank.
     */
    NDArray gather_on_root(const int iproot) const;

    /**
     * @brief Scatter a full host matrix into the frozen distribution.
     *
     * The MPI implementation is collective: every rank must call scatter()
     * with a complete host matrix of the global shape. Only the matrix supplied
     * by communicator rank zero is read; matrices supplied by other ranks are
     * ignored.
     */
    void scatter(const NDArray& matglb);

    /** Access an owned element by global row and column indices. */
    T& atglb(const int irglb, const int icglb)
    {
        assert(irglb >= 0);
        assert(irglb < nrglb());
        assert(icglb >= 0);
        assert(icglb < ncglb());

        if (distrib_ == nullptr)
        {
            return atloc(irglb, icglb);
        }

        const bool owns = distrib_->in_this_processor(irglb, icglb);
        assert(owns);

        const int irloc = distrib_->irloc(irglb);
        const int icloc = distrib_->icloc(icglb);
        assert(irloc != -1);
        assert(icloc != -1);
        return atloc(irloc, icloc);
    }

    /** Access an owned element by global row and column indices. */
    const T& atglb(const int irglb, const int icglb) const
    {
        assert(irglb >= 0);
        assert(irglb < nrglb());
        assert(icglb >= 0);
        assert(icglb < ncglb());

        if (distrib_ == nullptr)
        {
            return atloc(irglb, icglb);
        }

        const bool owns = distrib_->in_this_processor(irglb, icglb);
        assert(owns);

        const int irloc = distrib_->irloc(irglb);
        const int icloc = distrib_->icloc(icglb);
        assert(irloc != -1);
        assert(icloc != -1);
        return atloc(irloc, icloc);
    }

    /** Access an element by local row and column indices. */
    T& atloc(const int irloc, const int icloc)
    {
        assert(irloc >= 0);
        assert(irloc < nrloc());
        assert(icloc >= 0);
        assert(icloc < ncloc());

        std::size_t offset = 0;
        if (distrib_ == nullptr)
        {
            offset = static_cast<std::size_t>(irloc)
                    * static_cast<std::size_t>(ncloc_)
                + static_cast<std::size_t>(icloc);
        }
        else
        {
            const int lld = distrib_->lld();
            offset = static_cast<std::size_t>(icloc)
                    * static_cast<std::size_t>(lld)
                + static_cast<std::size_t>(irloc);
        }
        assert(offset < storage_.size());
        return storage_[offset];
    }

    /** Access an element by local row and column indices. */
    const T& atloc(const int irloc, const int icloc) const
    {
        assert(irloc >= 0);
        assert(irloc < nrloc());
        assert(icloc >= 0);
        assert(icloc < ncloc());

        std::size_t offset = 0;
        if (distrib_ == nullptr)
        {
            offset = static_cast<std::size_t>(irloc)
                    * static_cast<std::size_t>(ncloc_)
                + static_cast<std::size_t>(icloc);
        }
        else
        {
            const int lld = distrib_->lld();
            offset = static_cast<std::size_t>(icloc)
                    * static_cast<std::size_t>(lld)
                + static_cast<std::size_t>(irloc);
        }
        assert(offset < storage_.size());
        return storage_[offset];
    }

    /** Return whether this processor owns the given global indices. */
    bool in_this_processor(const int irglb, const int icglb) const
    {
        assert(irglb >= 0);
        assert(irglb < nrglb());
        assert(icglb >= 0);
        assert(icglb < ncglb());
        if (distrib_ == nullptr)
        {
            return true;
        }
        return distrib_->in_this_processor(irglb, icglb);
    }

    T* data() { return storage_.data(); }
    const T* data() const { return storage_.data(); }

    int nrglb() const { return distrib_ == nullptr ? nrloc_ : distrib_->nrglb(); }
    int ncglb() const { return distrib_ == nullptr ? ncloc_ : distrib_->ncglb(); }
    int nrloc() const { return distrib_ == nullptr ? nrloc_ : distrib_->nrloc(); }
    int ncloc() const { return distrib_ == nullptr ? ncloc_ : distrib_->ncloc(); }
    std::size_t size() const
    {
        return distrib_ == nullptr
            ? static_cast<std::size_t>(nrloc_) * static_cast<std::size_t>(ncloc_)
            : distrib_->size();
    }

    /**
     * @brief Number of T elements in the physical local storage buffer.
     *
     * This includes column-major padding and the minimum element reserved on
     * ranks that own no local rows or columns.
     */
    std::size_t local_storage_size() const { return storage_.size(); }
    int nrproc() const { return distrib_ == nullptr ? 1 : distrib_->nrproc(); }
    int ncproc() const { return distrib_ == nullptr ? 1 : distrib_->ncproc(); }
    int irproc() const { return distrib_ == nullptr ? 0 : distrib_->irproc(); }
    int icproc() const { return distrib_ == nullptr ? 0 : distrib_->icproc(); }
    int lld() const { return distrib_ == nullptr ? ncloc_ : distrib_->lld(); }

    const int* desc() const
    {
#ifdef __MPI
        return distrib_ != nullptr ? distrib_->desc() : nullptr;
#else
        return nullptr;
#endif
    }

    int ictxt() const
    {
#ifdef __MPI
        return distrib_ != nullptr ? distrib_->ictxt() : -1;
#else
        return -1;
#endif
    }

    const std::shared_ptr<const Parallel_Blkcyclic>& distribution() const { return distrib_; }
    bool is_distributed() const { return distrib_ != nullptr; }

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
    void assert_compatible_layout(const NDArray<U, 2>& other) const
    {
        if (nrglb() != other.nrglb()
            || ncglb() != other.ncglb())
        {
            throw std::invalid_argument("matrix global shapes are incompatible");
        }

        const bool distributed = distrib_ != nullptr;
        const bool distributed_ = other.distribution() != nullptr;
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

    void initialize_distributed_storage(
        const std::shared_ptr<const Parallel_Blkcyclic>& distrib)
    {
        if (distrib == nullptr)
        {
            throw std::invalid_argument("distribution must not be null");
        }
        distrib_ = distrib;
        initialize_distributed_storage();
    }

    void initialize_distributed_storage()
    {
        nrloc_ = distrib_->nrloc();
        ncloc_ = distrib_->ncloc();
        assert(nrloc_ >= 0);
        assert(ncloc_ >= 0);

        const int lld_distrib = distrib_->lld();
        const int lld = std::max(lld_distrib, 1);
        const std::size_t ld = static_cast<std::size_t>(lld);
        // The actual BLACS/LAPACK buffer is column-major. The row-major extent
        // is kept in the maximum as a guard for local logical accesses.
        const std::size_t column_major_elements
            = ld * static_cast<std::size_t>(ncloc_);
        const std::size_t row_major_elements
            = static_cast<std::size_t>(nrloc_)
                * static_cast<std::size_t>(ncloc_);
        const std::size_t logical_extent
            = std::max(column_major_elements, row_major_elements);
        const std::size_t empty_rank_minimum = 1;
        const std::size_t storage_elements
            = std::max(logical_extent, empty_rank_minimum);
        storage_.resize(storage_elements);
    }

    int nrloc_ = 0;
    int ncloc_ = 0;
    std::vector<T> storage_;
    std::shared_ptr<const Parallel_Blkcyclic> distrib_;
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
