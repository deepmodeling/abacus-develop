#ifndef ABACUS_SOURCE_SOURCE_BASE_MATRIX2D_VIEW_H
#define ABACUS_SOURCE_SOURCE_BASE_MATRIX2D_VIEW_H

#include "source_base/parallel_blkcyc.h"

#include <cassert>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <utility>

/**
 * @brief Non-owning two-dimensional matrix view.
 *
 * Host views use row-major storage. Distributed views borrow a frozen
 * block-cyclic distribution and use its local column-major layout.
 */
template <typename T>
class Matrix2dView
{
  public:
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
          distributed_(other.is_distributed()),
          desc_(other.desc()),
          dist_(other.distribution())
    {
    }

    Matrix2dView(T* data, const int nrglb, const int ncglb)
        : data_(data),
          nrglb_(nrglb),
          ncglb_(ncglb),
          nrloc_(nrglb),
          ncloc_(ncglb),
          lld_(ncglb),
          distributed_(false),
          desc_(nullptr)
    {
        assert(nrglb_ >= 0);
        assert(ncglb_ >= 0);
    }

#ifdef __MPI
    Matrix2dView(T* data,
                 std::shared_ptr<const Parallel_Blkcyclic> dist)
        : data_(data),
          nrglb_(checked_distribution(dist).nrglb()),
          ncglb_(checked_distribution(dist).ncglb()),
          nrloc_(checked_distribution(dist).nrloc()),
          ncloc_(checked_distribution(dist).ncloc()),
          lld_(checked_distribution(dist).lld()),
          distributed_(true),
          desc_(dist->desc()),
          dist_(std::move(dist))
    {
    }

    Matrix2dView(T* data,
                 const int nrloc,
                 const int ncloc,
                 const int* desc)
        : data_(data),
          nrglb_(checked_descriptor(desc)[2]),
          ncglb_(checked_descriptor(desc)[3]),
          nrloc_(nrloc),
          ncloc_(ncloc),
          lld_(checked_descriptor(desc)[8]),
          distributed_(true),
          desc_(desc)
    {
    }
#endif

    T* data() const { return data_; }
    int nrglb() const { return nrglb_; }
    int ncglb() const { return ncglb_; }
    int nrloc() const { return nrloc_; }
    int ncloc() const { return ncloc_; }
    int lld() const { return lld_; }
    std::size_t size() const
    {
        return static_cast<std::size_t>(nrloc_)
            * static_cast<std::size_t>(ncloc_);
    }

    bool is_distributed() const { return distributed_; }
    const std::shared_ptr<const Parallel_Blkcyclic>&
    distribution() const
    {
        return dist_;
    }

    const int* desc() const
    {
        return desc_;
    }

    T& at_local(const int i, const int j) const
    {
        assert(i >= 0);
        assert(i < nrloc_);
        assert(j >= 0);
        assert(j < ncloc_);
        const std::size_t offset = dist_ == nullptr
            ? static_cast<std::size_t>(i) * static_cast<std::size_t>(ncloc_)
                + static_cast<std::size_t>(j)
            : static_cast<std::size_t>(j) * static_cast<std::size_t>(lld_)
                + static_cast<std::size_t>(i);
        return data_[offset];
    }

  private:
    static const Parallel_Blkcyclic& checked_distribution(
        const std::shared_ptr<const Parallel_Blkcyclic>& dist)
    {
        if (dist == nullptr)
        {
            throw std::invalid_argument("matrix view distribution must not be null");
        }
        return *dist;
    }

    static const int* checked_descriptor(const int* desc)
    {
        if (desc == nullptr)
        {
            throw std::invalid_argument("matrix view descriptor must not be null");
        }
        return desc;
    }

    T* data_;
    int nrglb_;
    int ncglb_;
    int nrloc_;
    int ncloc_;
    int lld_;
    bool distributed_;
    const int* desc_;
    std::shared_ptr<const Parallel_Blkcyclic> dist_;
};

#endif
