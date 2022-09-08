#ifndef PARALLEL_ORBITALS_H
#define PARALLEL_ORBITALS_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "src_pdiag/pdiag_common.h"

/// These stucture packs the information of 2D-block-cyclic 
/// parallel distribution of basis, wavefunction and matrix.
struct Parallel_Orbitals
{

    Parallel_Orbitals();
    ~Parallel_Orbitals();
    
    /// map from global-index to local-index
    int* trace_loc_row;
    int* trace_loc_col;

    /// local size (nloc = nrow * ncol)
    int nrow;
	int ncol;
    long nloc;

    /// local size of bands, used for 2d wavefunction
    /// must divided on dim1 because of elpa interface
    int ncol_bands;
    
    /// ncol_bands*nrow
    long nloc_wfc;

    //ncol_bands*ncol_bands
    long nloc_Eij;

    /// block size
    int nb;

    /// the number of processors in each dimension of MPI_Cart structure
    int dim0;
    int dim1;
    
    int lastband_in_proc;
	int lastband_number; 

    ///---------------------------------------
    /// number of elements(basis-pairs) in this processon
    /// on all adjacent atoms-pairs(2D division)
    ///---------------------------------------
    int nnr;
	int *nlocdim;
	int *nlocstart;
    
#ifdef __MPI
    /// blacs info
    int blacs_ctxt;
    int desc[9];    //for matrix, nlocal*nlocal
    int desc_wfc[9]; //for wfc, nlocal*nbands
    int desc_Eij[9]; // for Eij in TDDFT, nbands*nbands
    int desc_wfc1[9]; // for wfc^T in TDDFT, nbands*nlocal
    /// communicator for 2D-block
    MPI_Comm comm_2D;
#endif

    /// only used in hpseps-diago
    int nspin = 1;
    int* loc_sizes;
    int loc_size;
    bool alloc_Z_LOC; //xiaohui add 2014-12-22
    double** Z_LOC; //xiaohui add 2014-06-19

    /// used in hpseps, dftu and exx
    LocalMatrix MatrixInfo;

    // test parameter
    int testpb;

    /// check whether a basis element is in this processor
    /// (check whether local-index > 0 )
    bool in_this_processor(const int& iw1_all, const int& iw2_all) const;

};


#endif
