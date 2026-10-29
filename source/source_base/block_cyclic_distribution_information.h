#ifndef ABACUS_SOURCE_SOURCE_BASE_BLOCK_CYCLIC_DISTRIBUTION_INFORMATION_H
#define ABACUS_SOURCE_SOURCE_BASE_BLOCK_CYCLIC_DISTRIBUTION_INFORMATION_H

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

#ifdef __MPI
#include <mpi.h>
#endif

// Responsibilities:
// 1. Describe a block-cyclic layout: process grid, block size, local shape,
//    global/local index maps, and owner queries.
// 2. Own or borrow the BLACS context and maintain the ScaLAPACK descriptor.
// It does not own physical matrices, atom ownership, or trace information.
class BlockCyclicDistributionInformation
{
  public:
    BlockCyclicDistributionInformation() = default;
    ~BlockCyclicDistributionInformation();

    BlockCyclicDistributionInformation(const BlockCyclicDistributionInformation&) = delete;
    BlockCyclicDistributionInformation& operator=(const BlockCyclicDistributionInformation&) = delete;
    BlockCyclicDistributionInformation(BlockCyclicDistributionInformation&& rhs) noexcept;
    BlockCyclicDistributionInformation& operator=(BlockCyclicDistributionInformation&& rhs) noexcept;

    int get_number_of_local_rows() const
    {
        return number_of_local_rows_;
    }

    int get_number_of_local_cols() const
    {
        return number_of_local_cols_;
    }

    int get_number_of_global_rows() const;
    int get_number_of_global_cols() const;

    int64_t get_number_of_local_elements() const
    {
        return number_of_local_elements_;
    }

    int get_block_size() const
    {
        return block_size_;
    }

    int get_grid_rows() const
    {
        return grid_rows_;
    }

    int get_grid_cols() const
    {
        return grid_cols_;
    }

    int get_coord_row() const
    {
        return coord_[0];
    }

    int get_coord_col() const
    {
        return coord_[1];
    }

    int global_to_local_row(const int global_row) const
    {
        return global_to_local_row_[global_row];
    }

    int global_to_local_col(const int global_col) const
    {
        return global_to_local_col_[global_col];
    }

    int local_to_global_row(const int local_row) const
    {
        return local_to_global_row_[local_row];
    }

    int local_to_global_col(const int local_col) const
    {
        return local_to_global_col_[local_col];
    }

    bool owns_index(const std::size_t iw1_all, const std::size_t iw2_all) const;
    int owner_of_index(const std::size_t iw1_all, const std::size_t iw2_all) const;

    bool owns_blacs_index(const int iprow, const int ipcol) const
    {
        return iprow == coord_[0] && ipcol == coord_[1];
    }

    void set_coord(const int coord_row, const int coord_col)
    {
        coord_[0] = coord_row;
        coord_[1] = coord_col;
    }

    bool is_serial() const
    {
        return is_serial_;
    }

    void initialize_serial(const int global_rows, const int global_cols);

#ifdef __MPI
    int initialize(const int grid_rows, const int grid_cols, const int block_size, const MPI_Comm comm);
    int initialize(const int grid_rows,
                   const int grid_cols,
                   const int block_size,
                   const MPI_Comm comm,
                   const bool swap_grid);
    int bind_descriptor(const int grid_rows,
                        const int grid_cols,
                        const int block_size,
                        const int blacs_context);

    const int* get_descriptor() const
    {
        return descriptor_.data();
    }

    int get_blacs_context() const
    {
        return blacs_context_;
    }

    MPI_Comm get_mpi_communicator() const;

    // DEPRECATED-ALIAS: get_mutable_descriptor() is only for tests/specializations.
    // DELETION-COMMIT: remove after descriptor writes use initialize/bind_descriptor; trigger: grep count == 0.
    int* get_mutable_descriptor()
    {
        return descriptor_.data();
    }
#endif

    int get_row_size() const
    {
        return get_number_of_local_rows();
    }

    int get_col_size() const
    {
        return get_number_of_local_cols();
    }

    int get_global_row_size() const
    {
        return get_number_of_global_rows();
    }

    int get_global_col_size() const
    {
        return get_number_of_global_cols();
    }

    int64_t get_local_size() const
    {
        return get_number_of_local_elements();
    }

    int get_dim0() const
    {
        return get_grid_rows();
    }

    int get_dim1() const
    {
        return get_grid_cols();
    }

    int global2local_row(const std::size_t igr) const
    {
        return global_to_local_row(static_cast<int>(igr));
    }

    int global2local_col(const std::size_t igc) const
    {
        return global_to_local_col(static_cast<int>(igc));
    }

    int local2global_row(const int ilr) const
    {
        return local_to_global_row(ilr);
    }

    int local2global_col(const int ilc) const
    {
        return local_to_global_col(ilc);
    }

    bool in_this_processor(const std::size_t iw1_all, const std::size_t iw2_all) const
    {
        return owns_index(iw1_all, iw2_all);
    }

    int owner_processor(const std::size_t iw1_all, const std::size_t iw2_all) const
    {
        return owner_of_index(iw1_all, iw2_all);
    }

    bool blacs_in_this_processor(const int iprow, const int ipcol) const
    {
        return owns_blacs_index(iprow, ipcol);
    }

    void set_serial(const int global_rows, const int global_cols)
    {
        initialize_serial(global_rows, global_cols);
    }

#ifdef __MPI
    int init(const int grid_rows,
             const int grid_cols,
             const int block_size,
             const MPI_Comm comm,
             const bool swap_grid = false);
    int set(const int grid_rows,
            const int grid_cols,
            const int block_size,
            const int blacs_context);

    const int* get_desc() const
    {
        return get_descriptor();
    }

    int get_blacs_ctxt() const
    {
        return get_blacs_context();
    }

    MPI_Comm comm() const
    {
        return get_mpi_communicator();
    }
#endif

  private:
    static constexpr int kDescriptorEntryCount = 9;
    static constexpr int kDescriptorBlacsContextIndex = 1;
    static constexpr int kDescriptorGlobalRowsIndex = 2;
    static constexpr int kDescriptorGlobalColsIndex = 3;
    static constexpr int kDescriptorLeadingDimIndex = 8;
    static constexpr int kUnownedBlacsContext = -1;
    static constexpr int kCoordinateDimension = 2;
    static constexpr int kBlacsWhatSystemContext = 10;

    int number_of_local_rows_ = 0;
    int number_of_local_cols_ = 0;
    int64_t number_of_local_elements_ = 0;
    int block_size_ = 1;
    int grid_rows_ = 0;
    int grid_cols_ = 0;

#ifdef __MPI
    int blacs_context_ = kUnownedBlacsContext;
    std::array<int, kDescriptorEntryCount> descriptor_ = {};
#endif

    std::vector<int> global_to_local_row_;
    std::vector<int> global_to_local_col_;
    std::vector<int> local_to_global_row_;
    std::vector<int> local_to_global_col_;
    std::array<int, kCoordinateDimension> coord_ = {{kUnownedBlacsContext, kUnownedBlacsContext}};
    bool is_serial_ = false;

#ifdef __MPI
    bool owns_blacs_context_ = false;

    void initialize_blacs_grid(const MPI_Comm comm, const bool swap_grid);
    void set_distribution_info(const int global_rows, const int global_cols, const int block_size);
    void release_blacs_grid();
#endif
};

#endif
