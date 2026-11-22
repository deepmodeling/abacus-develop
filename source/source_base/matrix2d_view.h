#ifndef ABACUS_SOURCE_SOURCE_BASE_MATRIX2D_VIEW_H
#define ABACUS_SOURCE_SOURCE_BASE_MATRIX2D_VIEW_H

/**
 * @file matrix2d_view.h
 * @brief Non-owning two-dimensional views over host or distributed storage.
 *
 * A `Matrix2dView<T>` is a lightweight descriptor for an existing buffer. It
 * can wrap caller-owned host data or the local buffer of a frozen
 * `Parallel_Blkcyclic` distribution. Host views use logical row-major indexing.
 * Distributed views use the local column-major layout reported by the
 * distribution, including serial one-by-one process grids.
 *
 * The view never owns `data()`. The caller must keep the data buffer alive and
 * unchanged in shape for the entire view lifetime. A distributed view shares
 * ownership of its distribution handle so the descriptor metadata remains
 * valid. For distributed storage, `lld()` is the column-major leading
 * dimension; `atloc(i, j)` addresses column `j` at offset `j * lld() + i`.
 * Local bounds violations are `assert()` failures rather than exceptions.
 */

#include "source_base/parallel_blkcyc.h"

#include <cassert>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <utility>

/** @brief Non-owning two-dimensional matrix view. */
template <typename T>
class Matrix2dView
{
  public:
    /**
     * @brief Convert a mutable view to a read-only view of the same buffer.
     *
     * The source view and its owner must outlive the converted view. No data
     * elements are copied.
     */
    template <
        typename U,
        typename std::enable_if<std::is_same<T, const U>::value, int>::type = 0>
    Matrix2dView(const Matrix2dView<U>& other)
        : data_(other.data()),
          nrglb_(other.nrglb()),
          ncglb_(other.ncglb()),
          nrloc_(other.nrloc()),
          ncloc_(other.ncloc()),
          lld_(other.lld()),
          desc_(other.desc()),
          distrib_(other.distribution())
    {
    }

    /**
     * @brief Construct a host view over an existing row-major buffer.
     *
     * @param data Caller-owned buffer with at least `nrglb * ncglb` elements.
     * @param nrglb Number of global rows.
     * @param ncglb Number of global columns.
     */
    Matrix2dView(T* data, const int nrglb, const int ncglb)
        : data_(data),
          nrglb_(nrglb),
          ncglb_(ncglb),
          nrloc_(nrglb),
          ncloc_(ncglb),
          lld_(ncglb),
          desc_(nullptr),
          distrib_(nullptr)
    {
        assert(nrglb_ >= 0);
        assert(ncglb_ >= 0);
    }

    /**
     * @brief Construct a view over a frozen block-cyclic local buffer.
     *
     * @param data Caller-owned column-major local buffer.
     * @param distrib Distribution describing the local shape and leading
     * dimension. The view shares ownership of this handle.
     */
    Matrix2dView(T* data,
                 std::shared_ptr<const Parallel_Blkcyclic> distrib)
        : data_(data),
          nrglb_(0),
          ncglb_(0),
          nrloc_(0),
          ncloc_(0),
          lld_(0),
          desc_(nullptr),
          distrib_(std::move(distrib))
    {
        if (distrib_ == nullptr)
        {
            throw std::invalid_argument(
                "matrix view distribution must not be null");
        }
        const Parallel_Blkcyclic& distribution = *distrib_;
        nrglb_ = distribution.nrglb();
        ncglb_ = distribution.ncglb();
        nrloc_ = distribution.nrloc();
        ncloc_ = distribution.ncloc();
        lld_ = distribution.lld();
#ifdef __MPI
        desc_ = distribution.desc();
#endif
    }

    /**
     * @brief Return the caller-owned data pointer.
     *
     * The pointer is writable only when `T` is a non-const type.
     */
    T* data() const { return data_; }
    /** @brief Return the global row count. */
    int nrglb() const { return nrglb_; }
    /** @brief Return the global column count. */
    int ncglb() const { return ncglb_; }
    /** @brief Return the local row count. */
    int nrloc() const { return nrloc_; }
    /** @brief Return the local column count. */
    int ncloc() const { return ncloc_; }
    /**
     * @brief Return the column-major leading dimension.
     *
     * Host views use the column count. Distributed views use the value from
     * their frozen distribution and may include padding between columns.
     */
    int lld() const { return lld_; }
    /** @brief Return the number of logical local entries. */
    std::size_t size() const
    {
        return static_cast<std::size_t>(nrloc_)
            * static_cast<std::size_t>(ncloc_);
    }

    /**
     * @brief Return whether a distribution handle is attached.
     *
     * This includes serial frozen distributions; it does not imply that more
     * than one MPI rank participates.
     */
    bool is_distributed() const { return distrib_ != nullptr; }
    /**
     * @brief Return the attached distribution, or null for a host view.
     *
     * The returned shared pointer keeps the distribution metadata alive while
     * it remains in scope. It does not own `data()`.
     */
    const std::shared_ptr<const Parallel_Blkcyclic>&
    distribution() const
    {
        return distrib_;
    }

    /**
     * @brief Return the ScaLAPACK descriptor for a distributed MPI view.
     *
     * Returns null for host views and for all views in non-MPI builds.
     */
    const int* desc() const
    {
        return desc_;
    }

    /**
     * @brief Access a local element without copying.
     *
     * Host views use row-major addressing. Distributed views use
     * `j * lld() + i` column-major addressing. Index violations are checked
     * with `assert()`, so they are not exceptions.
     */
    T& atloc(const int i, const int j) const
    {
        assert(i >= 0);
        assert(i < nrloc_);
        assert(j >= 0);
        assert(j < ncloc_);
        const std::size_t offset = distrib_ == nullptr
            ? static_cast<std::size_t>(i) * static_cast<std::size_t>(ncloc_)
                + static_cast<std::size_t>(j)
            : static_cast<std::size_t>(j) * static_cast<std::size_t>(lld_)
                + static_cast<std::size_t>(i);
        return data_[offset];
    }

  private:
    T* data_;
    int nrglb_;
    int ncglb_;
    int nrloc_;
    int ncloc_;
    int lld_;
    const int* desc_;
    std::shared_ptr<const Parallel_Blkcyclic> distrib_;
};

#endif
