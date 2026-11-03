#ifndef ABACUS_SOURCE_SOURCE_BASE_MODULE_CONTAINER_NDARRAY_DISTRIBUTED_NDARRAY_H
#define ABACUS_SOURCE_SOURCE_BASE_MODULE_CONTAINER_NDARRAY_DISTRIBUTED_NDARRAY_H

#include "source_base/block_cyclic_distribution_information.h"
#include "source_base/ndarray.h"

#include <cassert>
#include <cstddef>
#include <memory>
#include <utility>

namespace distributed_ndarray_detail
{
constexpr int DESCRIPTOR_BLACS_CONTEXT_INDEX = 1;
constexpr int DESCRIPTOR_LEADING_DIM_INDEX = 8;
constexpr int DESCRIPTOR_ENTRY_COUNT = 9;
} // namespace distributed_ndarray_detail

enum class StorageMajorOrder
{
    RowMajor,
    ColumnMajor
};

/// A read-only binding to a block-cyclic distribution. The distribution is
/// shared because multiple matrices commonly use the same layout. The unbound
/// host-side container is NDArray<T>; this type always requires a distribution.
template <typename T, typename Tdata = double>
class DistributedNDArray
{
  public:
    DistributedNDArray(std::shared_ptr<const BlockCyclicDistributionInformation> distribution,
                       StorageMajorOrder major_order,
                       int leading_dimension)
        : distribution_(distribution),
          local_block_(1, 1),
          major_order_(major_order),
          leading_dimension_(leading_dimension)
    {
        assert(distribution_ != nullptr);

        const int nrloc = distribution_->nrloc();
        const int ncloc = distribution_->ncloc();
        assert(nrloc >= 0);
        assert(ncloc >= 0);

        const int minimum_leading_dimension = nrloc > 1 ? nrloc : 1;
        assert(leading_dimension_ == minimum_leading_dimension);

#ifdef __MPI
        const int* const descriptor = distribution_->desc();
        assert(descriptor != nullptr);
        assert(leading_dimension_
               == descriptor[distributed_ndarray_detail::DESCRIPTOR_LEADING_DIM_INDEX]);
#endif

        int nrstore = leading_dimension_;
        int ncstore = ncloc;
        if (major_order_ == StorageMajorOrder::RowMajor)
        {
            nrstore = nrloc;
            if (ncstore < leading_dimension_)
            {
                ncstore = leading_dimension_;
            }
        }

        NDArray<T> local_block(nrstore, ncstore);
        const std::size_t nelem = static_cast<std::size_t>(nrstore) * static_cast<std::size_t>(ncstore);
        assert(nelem == local_block.size());
        local_block_ = std::move(local_block);
    }

    ~DistributedNDArray() = default;

    DistributedNDArray(const DistributedNDArray&) = delete;
    DistributedNDArray& operator=(const DistributedNDArray&) = delete;

    DistributedNDArray(DistributedNDArray&& other) noexcept
        : distribution_(std::move(other.distribution_)),
          local_block_(std::move(other.local_block_)),
          major_order_(other.major_order_),
          leading_dimension_(other.leading_dimension_)
    {
    }

    DistributedNDArray& operator=(DistributedNDArray&& other) noexcept
    {
        if (this != &other)
        {
            distribution_ = std::move(other.distribution_);
            local_block_ = std::move(other.local_block_);
            major_order_ = other.major_order_;
            leading_dimension_ = other.leading_dimension_;
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
            flat_index = static_cast<std::ptrdiff_t>(local_col) * leading_dimension_ + local_row;
        }
        else
        {
            flat_index = static_cast<std::ptrdiff_t>(local_row) * leading_dimension_ + local_col;
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
            flat_index = static_cast<std::ptrdiff_t>(local_col) * leading_dimension_ + local_row;
        }
        else
        {
            flat_index = static_cast<std::ptrdiff_t>(local_row) * leading_dimension_ + local_col;
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

    /// Returns the local block beginning. Its leading dimension is
    /// leading_dimension().
    T* data()
    {
        return local_block_.data();
    }

    const T* data() const
    {
        return local_block_.data();
    }

    const std::shared_ptr<const BlockCyclicDistributionInformation>& distribution() const
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

    int leading_dimension() const
    {
        return leading_dimension_;
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
    std::shared_ptr<const BlockCyclicDistributionInformation> distribution_;
    NDArray<T> local_block_;
    StorageMajorOrder major_order_;
    int leading_dimension_ = 0;
};

#endif // ABACUS_SOURCE_SOURCE_BASE_MODULE_CONTAINER_NDARRAY_DISTRIBUTED_NDARRAY_H
