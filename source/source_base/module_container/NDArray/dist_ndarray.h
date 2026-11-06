#ifndef ABACUS_SOURCE_SOURCE_BASE_MODULE_CONTAINER_NDARRAY_DIST_NDARRAY_H
#define ABACUS_SOURCE_SOURCE_BASE_MODULE_CONTAINER_NDARRAY_DIST_NDARRAY_H

#include "source_base/block_cyclic_distribution.h"
#include "source_base/ndarray.h"

#include <cassert>
#include <cstddef>
#include <memory>
#include <utility>

enum class StorageMajorOrder
{
    RowMajor,
    ColumnMajor
};

/// A read-only binding to a block-cyclic distribution. The distribution is
/// shared because multiple matrices commonly use the same layout. The unbound
/// host-side container is NDArray<T>; this type always requires a distribution.
template <typename T, typename Tdata = double>
class DistNDArray
{
  public:
    DistNDArray(std::shared_ptr<const BlockCyclicDistribution> distribution,
                StorageMajorOrder major_order,
                int ld)
        : distribution_(distribution),
          local_block_(1, 1),
          major_order_(major_order),
          ld_(ld)
    {
        assert(distribution_ != nullptr);

        const int nrloc = distribution_->nrloc();
        const int ncloc = distribution_->ncloc();
        assert(nrloc >= 0);
        assert(ncloc >= 0);

        const int nloc
            = major_order_ == StorageMajorOrder::ColumnMajor ? nrloc : ncloc;
        const int min_ld = nloc > 1 ? nloc : 1;
        assert(ld_ == min_ld);

#ifdef __MPI
        if (!distribution_->is_serial())
        {
            const int* const descriptor = distribution_->desc();
            assert(descriptor != nullptr);
            const int desc_ld = descriptor[DESC_LLD];
            const int lld = nrloc > 1 ? nrloc : 1;
            assert(desc_ld == lld);
        }
#endif

        int nrstore = ld_;
        int ncstore = ncloc;
        if (major_order_ == StorageMajorOrder::RowMajor)
        {
            nrstore = nrloc;
            ncstore = ld_;
        }

        NDArray<T> local_block(nrstore, ncstore);
        const std::size_t nelem = static_cast<std::size_t>(nrstore) * static_cast<std::size_t>(ncstore);
        assert(nelem == local_block.size());
        local_block_ = std::move(local_block);
    }

    ~DistNDArray() = default;

    DistNDArray(const DistNDArray&) = delete;
    DistNDArray& operator=(const DistNDArray&) = delete;

    DistNDArray(DistNDArray&& other) noexcept
        : distribution_(std::move(other.distribution_)),
          local_block_(std::move(other.local_block_)),
          major_order_(other.major_order_),
          ld_(other.ld_)
    {
    }

    DistNDArray& operator=(DistNDArray&& other) noexcept
    {
        if (this != &other)
        {
            distribution_ = std::move(other.distribution_);
            local_block_ = std::move(other.local_block_);
            major_order_ = other.major_order_;
            ld_ = other.ld_;
        }
        return *this;
    }

    T& at_local(const int i, const int j)
    {
        assert(i >= 0);
        assert(i < distribution_->nrloc());
        assert(j >= 0);
        assert(j < distribution_->ncloc());

        std::ptrdiff_t iflat = 0;
        if (major_order_ == StorageMajorOrder::ColumnMajor)
        {
            iflat = static_cast<std::ptrdiff_t>(j) * ld_ + i;
        }
        else
        {
            iflat = static_cast<std::ptrdiff_t>(i) * ld_ + j;
        }

        assert(iflat >= 0);
        const std::size_t istore = static_cast<std::size_t>(iflat);
        assert(istore < local_block_.size());
        return local_block_.data()[istore];
    }

    const T& at_local(const int i, const int j) const
    {
        assert(i >= 0);
        assert(i < distribution_->nrloc());
        assert(j >= 0);
        assert(j < distribution_->ncloc());

        std::ptrdiff_t iflat = 0;
        if (major_order_ == StorageMajorOrder::ColumnMajor)
        {
            iflat = static_cast<std::ptrdiff_t>(j) * ld_ + i;
        }
        else
        {
            iflat = static_cast<std::ptrdiff_t>(i) * ld_ + j;
        }

        assert(iflat >= 0);
        const std::size_t istore = static_cast<std::size_t>(iflat);
        assert(istore < local_block_.size());
        return local_block_.data()[istore];
    }

    T& at_global(const int i, const int j)
    {
        assert(i >= 0);
        assert(i < distribution_->nrglb());
        assert(j >= 0);
        assert(j < distribution_->ncglb());

        const bool owned = distribution_->in_this_processor(i, j);
        assert(owned);

        const int il = distribution_->irloc(i);
        const int jl = distribution_->icloc(j);
        assert(il != -1);
        assert(jl != -1);
        return at_local(il, jl);
    }

    const T& at_global(const int i, const int j) const
    {
        assert(i >= 0);
        assert(i < distribution_->nrglb());
        assert(j >= 0);
        assert(j < distribution_->ncglb());

        const bool owned = distribution_->in_this_processor(i, j);
        assert(owned);

        const int il = distribution_->irloc(i);
        const int jl = distribution_->icloc(j);
        assert(il != -1);
        assert(jl != -1);
        return at_local(il, jl);
    }

    bool owns_global_index(const int i, const int j) const
    {
        assert(i >= 0);
        assert(i < distribution_->nrglb());
        assert(j >= 0);
        assert(j < distribution_->ncglb());

        return distribution_->in_this_processor(i, j);
    }

    /// Returns the beginning of contiguous storage. Use major_order() and
    /// ld() for raw indexing, or at_local() for logical indexing.
    T* data()
    {
        return local_block_.data();
    }

    const T* data() const
    {
        return local_block_.data();
    }

    const std::shared_ptr<const BlockCyclicDistribution>& distribution() const
    {
        return distribution_;
    }

    StorageMajorOrder major_order() const
    {
        return major_order_;
    }

    int ld() const
    {
        return ld_;
    }

    int nrloc() const
    {
        return distribution_->nrloc();
    }

    int ncloc() const
    {
        return distribution_->ncloc();
    }

  private:
    std::shared_ptr<const BlockCyclicDistribution> distribution_;
    NDArray<T> local_block_;
    StorageMajorOrder major_order_;
    int ld_ = 0;
    static constexpr int DESC_LLD = 8;
};

#endif // ABACUS_SOURCE_SOURCE_BASE_MODULE_CONTAINER_NDARRAY_DIST_NDARRAY_H
