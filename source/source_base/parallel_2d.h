#ifndef _PARALLEL_2D_H_
#define _PARALLEL_2D_H_

#include <cstdint>
#include <vector>

#include "source_base/parallel_comm.h"

/// @brief  This class packs the basic information of
/// 2D-block-cyclic parallel distribution of an arbitrary matrix.
class Parallel_2D
{
  public:
    Parallel_2D() = default;
    ~Parallel_2D() = default;

    Parallel_2D& operator=(Parallel_2D&& rhs) = default;
    Parallel_2D(Parallel_2D&& rhs) = default;

    /// number of local rows
    int get_row_size() const
    {
        return nrow;
    };

    /// number of local columns
    int get_col_size() const
    {
        return ncol;
    };

    /// number of global rows
    int get_global_row_size() const;

    /// number of global columns
    int get_global_col_size() const;

    /// number of local matrix elements
    int64_t get_local_size() const
    {
        return nloc;
    };

    /// get the local index of a global index (row)
    int global2local_row(const int igr) const
    {
        return global2local_row_[igr];
    }

    /// get the local index of a global index (col)
    int global2local_col(const int igc) const
    {
        return global2local_col_[igc];
    }

    /// get the global index of a local index (row)
    int local2global_row(const int ilr) const
    {
        return local2global_row_[ilr];
    }

    /// get the global index of a local index (col)
    int local2global_col(const int ilc) const
    {
        return local2global_col_[ilc];
    }

    /// check whether a global index is in this process
    bool in_this_processor(const int iw1_all, const int iw2_all) const;

    /// side length of 2d square block
    int get_block_size() const
    {
        return nb;
    };

    /// number of processes in row dimension of the MPI Cartesian grid
    int get_dim0() const
    {
        return dim0;
    };

    /// number of processes in column dimension of the MPI Cartesian grid
    int get_dim1() const
    {
        return dim1;
    };

    /// row coordinate in the BLACS grid
    int get_coord_row() const
    {
        return coord[0];
    };

    /// column coordinate in the BLACS grid
    int get_coord_col() const
    {
        return coord[1];
    };

    /// check whether a given BLACS grid coordinate is this process
    bool blacs_in_this_processor(const int iprow, const int ipcol) const
    {
        return iprow == coord[0] && ipcol == coord[1];
    };

    /// set process coordinate in the BLACS grid (intended for testing)
    void set_coord(const int row, const int col)
    {
        coord[0] = row;
        coord[1] = col;
    };

#ifdef __MPI

    /// ScaLAPACK descriptor
    const int* get_desc() const
    {
        return desc;
    };

    /**
     * @brief Initialize a BLACS grid with the given MPI communicator
     * and set up the info of a block-cyclic distribution.
     *
     */
    int init(const int mg,
             const int ng,
             const int nb, // square block is assumed
             const MPI_Comm comm,
             bool mode = false);

    /**
     * @brief Set up the info of a block-cyclic distribution using given
     * BLACS context.
     *
     */
    int set(const int mg,
            const int ng,
            const int nb, // square block is assumed
            const int blacs_ctxt);

    /// BLACS context
    int blacs_ctxt = -1;

    /// ScaLAPACK descriptor
    int desc[9] = {};

    MPI_Comm comm() const;
#endif

    void set_serial(const int mg, const int ng);

    // FIXME the following variables should be private, but they are
    // widely used in the code currently. Public visibility is kept
    // for now, but might be changed in the future.

    /// local size (nloc = nrow * ncol)
    int nrow = 0;
    int ncol = 0;
    int64_t nloc = 0;
    // NOTE: ScaLAPACK descriptors use int type for the number of rows and columns of
    // both the global and local matrices, so nrow & ncol have to be int type. Their
    // product, however, can exceed the range of int type.

    /// block size
    int nb = 1;

    /// number of processes in each dimension of the MPI Cartesian grid
    int dim0 = 0;
    int dim1 = 0;

  private:
    /// process coordinate in the BLACS grid
    int coord[2] = {-1, -1};
    /// whether to use the serial mode
    bool is_serial = false;

  protected:

    /// map from global index to local index
    std::vector<int> global2local_row_;
    std::vector<int> global2local_col_;

    /// map from local index to global index
    std::vector<int> local2global_row_;
    std::vector<int> local2global_col_;

#ifdef __MPI
    void _init_proc_grid(const MPI_Comm comm, const bool mode);
    void _set_dist_info(const int mg, const int ng, const int nb);
#endif
};
#endif
