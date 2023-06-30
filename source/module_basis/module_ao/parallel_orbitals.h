#ifndef PARALLEL_ORBITALS_H
#define PARALLEL_ORBITALS_H

#include "module_base/global_function.h"
#include "module_base/global_variable.h"

#ifdef __MPI
#include <mpi.h>
#endif

/// @brief  This structure packs the basic information of 2D-block-cyclic
/// parallel distribution of an arbitrary 2D Tensor.
struct Parallel_2D
{

    Parallel_2D();
    ~Parallel_2D();
    
    /// map from global-index to local-index
    int* trace_loc_row = nullptr;
    int* trace_loc_col = nullptr;

    /// map from local index to global index
    std::vector<int> row_set;				// Peize Lin change int* to vector 2022.08.03
    std::vector<int> col_set;

    /// local size (nloc = nrow * ncol)
    int nrow;
	int ncol;
    long nloc;
    
    /// block size
    /// default value of nb is 1,
    /// but can change to larger value from input.
    int nb = 1;

    /// the number of processors in each dimension of MPI_Cart structure
    int dim0;
    int dim1;
    /// the coordinate of current processor in each dimension of MPI_Cart structure
    int coord[2];

    /// test parameter
    int testpb;

#ifdef __MPI
    int blacs_ctxt;    ///< blacs info
    int desc[9];    ///<for matrix, nlocal*nlocal    
    MPI_Comm comm_2D;   ///<communicator for 2D-block
    /// create the 'comm_2D' stratege.
    void mpi_create_cart();

    /// set the map from local index to global index
    int set_local2global(const int& M_A, ///< global row size
        const int& N_A, ///< global col size
        std::ofstream& ofs_running,
        std::ofstream& ofs_warning);
#endif

    /// set the 2D-structure of processors in each dimension.
    /// dim0 and dim1 will be set as close to sqrt(nproc) as possible, 
    /// for example, if nproc = 12, dim0 = 3, dim1 = 4. 
    /// if mode==0, dim0 <= dim1; else, dim0 >= dim1.
    void set_proc_dim(const int& dsize, bool mode = 0);

    /// check whether a basis element is in this processor
    /// (check whether local-index > 0 )
    bool in_this_processor(const int& iw1_all, const int& iw2_all) const;

    void set_global2local(const int& M_A,
        const int& N_A,
        bool& div_2d,
        std::ofstream& ofs_running);

protected:
    /// set the map from local index to global index
    void init_global2local(const int& M_A, ///< global row size
        const int& N_A, ///< global col size
        std::ofstream& ofs_running);
};

/// These stucture packs the information of 2D-block-cyclic 
/// parallel distribution of basis, wavefunction and matrix.
struct Parallel_Orbitals : public Parallel_2D
{

    Parallel_Orbitals();
    ~Parallel_Orbitals();

    /// local size of bands, used for 2d wavefunction
    /// must divided on dim1 because of elpa interface
    int ncol_bands;
    int nrow_bands;
    
    /// ncol_bands*nrow
    long nloc_wfc;

    //ncol_bands*ncol_bands
    long nloc_Eij;

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
    int desc_wfc[9]; //for wfc, nlocal*nbands
    int desc_Eij[9]; // for Eij in TDDFT, nbands*nbands
    int desc_wfc1[9]; // for wfc^T in TDDFT, nbands*nlocal
    /// set the local size of wavefunction and Eij
    int set_nloc_wfc_Eij(const int& N_A,  ///< global col size
        std::ofstream& ofs_running,
        std::ofstream& ofs_warning);
#endif

    int nspin = 1;
    int* loc_sizes;
    int loc_size;

    // orbital index for each atom
    std::vector<int> atom_begin_row;
    std::vector<int> atom_begin_col;

    // set row and col begin index for each atom
    void set_atomic_trace(const int* iat2iwt, const int &nat, const int &nlocal);

    /**
     * @brief dimension getters for 2D-block-cyclic division of Hamiltonian matrix
     * get_col_size() : total number of columns of Hamiltonian matrix in this processor
     * get_row_size() : total number of rows of Hamiltonian matrix in this processor
     * get_col_size(iat) : number of columns of Hamiltonian matrix in atom iat
     * get_row_size(iat) : number of rows of Hamiltonian matrix in atom iat
    */
    int get_col_size()const;
    int get_row_size()const;
    int get_col_size(int iat) const;
    int get_row_size(int iat) const;

};


#endif
