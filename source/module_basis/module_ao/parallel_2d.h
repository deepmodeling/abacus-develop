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
    Parallel_2D();
    ~Parallel_2D();

    /// local size (nloc = nrow * ncol)
    int nrow;
    int ncol;
    long nloc;

    /// block size,
    /// the default value of nb is 1,
    /// but can change to larger value from input.
    int nb = 1;

    /// the number of processors in each dimension of MPI_Cart structure
    int dim0;
    int dim1;
    /// the coordinate of current processor in each dimension of MPI_Cart structure
    int coord[2];

    /// test parameter
    int testpb;

    /// total number of columns of matrix in this processor
    int get_col_size()const { return this->ncol; };

    /// total number of rows of matrix in this processor
    int get_row_size()const { return this->nrow; };

    /// total number of elements of matrix in this processor
    int get_local_size()const { return this->nloc; };

    /// get the local index of a global index (row)
    int global2local_row(const int& igr) const
    {
        return this->global2local_row_[igr];
    }

    /// get the local index of a global index (col)
    int global2local_col(const int& igc) const
    {
        return this->global2local_col_[igc];
    }

    /// get the global index of a local index (row)
    int local2global_row(const int& ilr) const
    {
        return this->local2global_row_[ilr];
    }

    /// get the global index of a local index (col)
    int local2global_col(const int& ilc) const
    {
        return this->local2global_col_[ilc];
    }

    /// check whether an element is in this processor
    /// (check whether local-index > 0 )
    bool in_this_processor(const int& iw1_all, const int& iw2_all) const;

    void set_block_size(const int& nb_in) { this->nb = nb_in; };
    int get_block_size()const { return this->nb; };

    /// Set the 2D-structure of processors in each dimension.
    /// dim0 and dim1 will be set as close to sqrt(nproc) as possible.
    /// For example: nproc = 12,
    /// if mode==0, d dim0 = 3, dim1 = 4; else, dim0 = 3, dim1 = 3.
    void set_proc_dim(const int& dsize, bool mode = 0);

#ifdef __MPI
    int blacs_ctxt;    ///< blacs info
    int desc[9];    ///<for matrix, nlocal*nlocal    
    MPI_Comm comm_2D;   ///<communicator for 2D-block

    /// create the 'comm_2D' stratege.
    void mpi_create_cart(const MPI_Comm& diag_world);

    /// set the map from local index to global index,
    /// and set local sizes (nrow, ncol, nloc) by the way
    int set_local2global(const int& M_A/**< global row size*/,
        const int& N_A/**< global col size*/,
        std::ofstream& ofs_running,
        std::ofstream& ofs_warning);

    ///@brief set the desc[9] of the 2D-block-cyclic distribution
    void set_desc(const int& gr/**< global row size*/,
        const int& gc/**< global col size*/,
        const int& lld/**< leading local dimension*/);
#else
    void set_serial(const int& M_A/**< global row size*/,
        const int& N_A/**< global col size*/);
#endif

    void set_global2local(const int& M_A,
        const int& N_A,
        const bool& div_2d,
        std::ofstream& ofs_running);

protected:

    /// map from global index to local index
    int* global2local_row_ = nullptr;
    int* global2local_col_ = nullptr;

    /// map from local index to global index
    std::vector<int> local2global_row_;				// Peize Lin change int* to vector 2022.08.03
    std::vector<int> local2global_col_;

    /// set the map from local index to global index
    void init_global2local(const int& M_A/**< global row size*/,
        const int& N_A/**< global col size*/,
        std::ofstream& ofs_running);
};
#endif