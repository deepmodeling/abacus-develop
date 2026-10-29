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
constexpr int kDescriptorBlacsContextIndex = 1;
constexpr int kDescriptorLeadingDimIndex = 8;
constexpr int kDescriptorEntryCount = 9;
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

        const int number_of_local_rows = distribution_->get_number_of_local_rows();
        const int number_of_local_cols = distribution_->get_number_of_local_cols();
        assert(number_of_local_rows >= 0);
        assert(number_of_local_cols >= 0);

        const int minimum_leading_dimension = number_of_local_rows > 1 ? number_of_local_rows : 1;
        assert(leading_dimension_ == minimum_leading_dimension);

#ifdef __MPI
        const int* const descriptor = distribution_->get_descriptor();
        assert(descriptor != nullptr);
        assert(leading_dimension_
               == descriptor[distributed_ndarray_detail::kDescriptorLeadingDimIndex]);
#endif

        int number_of_storage_rows = leading_dimension_;
        int number_of_storage_cols = number_of_local_cols;
        if (major_order_ == StorageMajorOrder::RowMajor)
        {
            number_of_storage_rows = number_of_local_rows;
            if (number_of_storage_cols < leading_dimension_)
            {
                number_of_storage_cols = leading_dimension_;
            }
        }

        NDArray<T> local_block(number_of_storage_rows, number_of_storage_cols);
        const std::size_t required_number_of_elements
            = static_cast<std::size_t>(number_of_storage_rows)
            * static_cast<std::size_t>(number_of_storage_cols);
        assert(required_number_of_elements == local_block.size());
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
        assert(local_row < distribution_->get_number_of_local_rows());
        assert(local_col >= 0);
        assert(local_col < distribution_->get_number_of_local_cols());

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
        assert(local_row < distribution_->get_number_of_local_rows());
        assert(local_col >= 0);
        assert(local_col < distribution_->get_number_of_local_cols());

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
        assert(global_row < distribution_->get_number_of_global_rows());
        assert(global_col >= 0);
        assert(global_col < distribution_->get_number_of_global_cols());

        const std::size_t global_row_index = static_cast<std::size_t>(global_row);
        const std::size_t global_col_index = static_cast<std::size_t>(global_col);
        const bool owns_index = distribution_->owns_index(global_row_index, global_col_index);
        assert(owns_index);

        const int local_row = distribution_->global_to_local_row(global_row);
        const int local_col = distribution_->global_to_local_col(global_col);
        assert(local_row != -1);
        assert(local_col != -1);
        return at_local(local_row, local_col);
    }

    const T& at_global(const int global_row, const int global_col) const
    {
        assert(global_row >= 0);
        assert(global_row < distribution_->get_number_of_global_rows());
        assert(global_col >= 0);
        assert(global_col < distribution_->get_number_of_global_cols());

        const std::size_t global_row_index = static_cast<std::size_t>(global_row);
        const std::size_t global_col_index = static_cast<std::size_t>(global_col);
        const bool owns_index = distribution_->owns_index(global_row_index, global_col_index);
        assert(owns_index);

        const int local_row = distribution_->global_to_local_row(global_row);
        const int local_col = distribution_->global_to_local_col(global_col);
        assert(local_row != -1);
        assert(local_col != -1);
        return at_local(local_row, local_col);
    }

    bool owns_global_index(const int global_row, const int global_col) const
    {
        assert(global_row >= 0);
        assert(global_row < distribution_->get_number_of_global_rows());
        assert(global_col >= 0);
        assert(global_col < distribution_->get_number_of_global_cols());

        const std::size_t global_row_index = static_cast<std::size_t>(global_row);
        const std::size_t global_col_index = static_cast<std::size_t>(global_col);
        return distribution_->owns_index(global_row_index, global_col_index);
    }

    /// Returns the local block beginning. Its leading dimension is
    /// get_leading_dimension().
    T* data()
    {
        return local_block_.data();
    }

    const T* data() const
    {
        return local_block_.data();
    }

    const std::shared_ptr<const BlockCyclicDistributionInformation>& get_distribution() const
    {
        return distribution_;
    }

    NDArray<T>& get_local_block()
    {
        return local_block_;
    }

    const NDArray<T>& get_local_block() const
    {
        return local_block_;
    }

    StorageMajorOrder get_major_order() const
    {
        return major_order_;
    }

    int get_leading_dimension() const
    {
        return leading_dimension_;
    }

    int get_number_of_local_rows() const
    {
        return distribution_->get_number_of_local_rows();
    }

    int get_number_of_local_cols() const
    {
        return distribution_->get_number_of_local_cols();
    }

  private:
    std::shared_ptr<const BlockCyclicDistributionInformation> distribution_;
    NDArray<T> local_block_;
    StorageMajorOrder major_order_;
    int leading_dimension_ = 0;
};

#endif // ABACUS_SOURCE_SOURCE_BASE_MODULE_CONTAINER_NDARRAY_DISTRIBUTED_NDARRAY_H
