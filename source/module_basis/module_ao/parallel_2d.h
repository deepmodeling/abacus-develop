#include <fstream>
#ifndef _PARALLEL_2D_H_
#define _PARALLEL_2D_H_
#include "module_base/global_function.h"
#include "module_base/global_variable.h"

#ifdef __MPI
#include <mpi.h>
#endif

/// @brief  This class packs the basic information of
/// 2D-block-cyclic parallel distribution of an arbitrary matrix.
class Parallel_2D
{
public:
    Parallel_2D() = default;
    ~Parallel_2D() = default;

    /// local number of columns
    int get_col_size()const { return ncol; };

    /// local number of rows
    int get_row_size()const { return nrow; };

    /// local number of matrix elements
    int64_t get_local_size() const { return nloc; };

    /// get the local index of a global index (row)
    int global2local_row(const int igr) const { return global2local_row_[igr]; }

    /// get the local index of a global index (col)
    int global2local_col(const int igc) const { return global2local_col_[igc]; }

    /// get the global index of a local index (row)
    int local2global_row(const int ilr) const { return local2global_row_[ilr]; }

    /// get the global index of a local index (col)
    int local2global_col(const int ilc) const { return local2global_col_[ilc]; }

    /// check whether a global index is in this process
    bool in_this_processor(const int iw1_all, const int iw2_all) const;

    /// side length of 2d square block
    int get_block_size() const { return nb; };

    //void set_block_size(const int& nb_in) { this->nb = nb_in; };

    /// Set the 2D-structure of processors in each dimension.
    /// dim0 and dim1 will be set as close to sqrt(nproc) as possible.
    /// For example: nproc = 12,
    /// if mode==0, d dim0 = 3, dim1 = 4; else, dim0 = 3, dim1 = 3.
    //void set_proc_dim(const int& dsize, bool mode = 0);

#ifndef __MPI
    void set_serial(int mg, int ng);
#else
    int set(
        const int mg,
        const int ng,
        const int nb, // block is assumed to be square
        const MPI_Comm comm
    );

    /// BLACS context
    int blacs_ctxt = -1;

    /// ScaLAPACK descriptor
    int desc[9] = {};

    /// 2D Cartesian MPI communicator
    MPI_Comm comm_2D = MPI_COMM_NULL;
#endif

    ///// set the map from local index to global index,
    ///// and set local sizes (nrow, ncol, nloc) by the way
    //int set_local2global(const int& M_A/**< global row size*/,
    //    const int& N_A/**< global col size*/,
    //    std::ofstream& ofs_running,
    //    std::ofstream& ofs_warning);

    /////@brief set the desc[9] of the 2D-block-cyclic distribution
    //void set_desc(const int& gr/**< global row size*/,
    //    const int& gc/**< global col size*/,
    //    const int& lld/**< leading local dimension*/,
    //    bool first_time = true/**< true: call `Cblacs_get`; false: use `this->blacs_ctxt`*/);

    //void set_global2local(const int& M_A,
    //    const int& N_A,
    //    const bool& div_2d,
    //    std::ofstream& ofs_running);

    // FIXME the following variables should be private, but they are
    // currently widely used in the code. Public visibility is kept
    // for now, and might be changed in the future.

    /// local size (nloc = nrow * ncol)
    int nrow = 0;
    int ncol = 0;
    int64_t nloc = 0;

    /// block size
    int nb = 1;

    /// number of processes in each dimension of MPI_Cart structure
    int dim0 = 0;
    int dim1 = 0;

    /// process coordinate in the MPI_Cart structure
    int coord[2] = {-1, -1};

    /// test parameter
    int testpb = 0;


protected:

    /// map from global index to local index
    std::vector<int> global2local_row_;
    std::vector<int> global2local_col_;

    /// map from local index to global index
    std::vector<int> local2global_row_;
    std::vector<int> local2global_col_;
    // Peize Lin change int* to vector 2022.08.03

    /// set the map from local index to global index
    //void init_global2local(const int& M_A/**< global row size*/,
    //    const int& N_A/**< global col size*/,
    //    std::ofstream& ofs_running);

    /// factorizes n = p * q such that p, q are closest and p <= q
    static void _fact(int n, int& p, int& q);
};
#endif
