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

        const int local_extent
            = major_order_ == StorageMajorOrder::ColumnMajor ? nrloc : ncloc;
        const int minimum_ld = local_extent > 1 ? local_extent : 1;
        assert(ld_ == minimum_ld);

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

    T& at_local(const int local_row, const int local_col)
    {
        assert(local_row >= 0);
        assert(local_row < distribution_->nrloc());
        assert(local_col >= 0);
        assert(local_col < distribution_->ncloc());

        std::ptrdiff_t flat_index = 0;
        if (major_order_ == StorageMajorOrder::ColumnMajor)
        {
            flat_index = static_cast<std::ptrdiff_t>(local_col) * ld_ + local_row;
        }
        else
        {
            flat_index = static_cast<std::ptrdiff_t>(local_row) * ld_ + local_col;
        }

        assert(flat_index >= 0);
        const std::size_t storage_index = static_cast<std::size_t>(flat_index);
        assert(storage_index < local_block_.size());
        return local_block_.data()[storage_index];
    }

    const T& at_local(const int local_row, const int local_col) const
    {
        assert(local_row >= 0);
        assert(local_row < distribution_->nrloc());
        assert(local_col >= 0);
        assert(local_col < distribution_->ncloc());

        std::ptrdiff_t flat_index = 0;
        if (major_order_ == StorageMajorOrder::ColumnMajor)
        {
            flat_index = static_cast<std::ptrdiff_t>(local_col) * ld_ + local_row;
        }
        else
        {
            flat_index = static_cast<std::ptrdiff_t>(local_row) * ld_ + local_col;
        }

        assert(flat_index >= 0);
        const std::size_t storage_index = static_cast<std::size_t>(flat_index);
        assert(storage_index < local_block_.size());
        return local_block_.data()[storage_index];
    }

    T& at_global(const int global_row, const int global_col)
    {
        assert(global_row >= 0);
        assert(global_row < distribution_->nrglb());
        assert(global_col >= 0);
        assert(global_col < distribution_->ncglb());

        const bool owned = distribution_->in_this_processor(global_row, global_col);
        assert(owned);

        const int local_row = distribution_->irloc(global_row);
        const int local_col = distribution_->icloc(global_col);
        assert(local_row != -1);
        assert(local_col != -1);
        return at_local(local_row, local_col);
    }

    const T& at_global(const int global_row, const int global_col) const
    {
        assert(global_row >= 0);
        assert(global_row < distribution_->nrglb());
        assert(global_col >= 0);
        assert(global_col < distribution_->ncglb());

        const bool owned = distribution_->in_this_processor(global_row, global_col);
        assert(owned);

        const int local_row = distribution_->irloc(global_row);
        const int local_col = distribution_->icloc(global_col);
        assert(local_row != -1);
        assert(local_col != -1);
        return at_local(local_row, local_col);
    }

    bool owns_global_index(const int global_row, const int global_col) const
    {
        assert(global_row >= 0);
        assert(global_row < distribution_->nrglb());
        assert(global_col >= 0);
        assert(global_col < distribution_->ncglb());

        return distribution_->in_this_processor(global_row, global_col);
    }

    /// Returns the local block beginning. Its leading dimension is ld().
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

    NDArray<T>& local_block()
    {
        return local_block_;
    }

    const NDArray<T>& local_block() const
    {
        return local_block_;
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
