#ifndef _PARALLEL_ORBITALS_H_
#define _PARALLEL_ORBITALS_H_
#include "parallel_2d.h"

/// This class packs the information of 2D-block-cyclic for LCAO code:
/// parallel distribution of basis, wavefunction and matrix.
class Parallel_Orbitals : public Parallel_2D
{
public:
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
    int set_nloc_wfc_Eij(const int& N_A/**< global row size*/,
        std::ofstream& ofs_running,
        std::ofstream& ofs_warning);

    ///@brief set the desc[9] of the 2D-block-cyclic distribution of wavefunction and Eij
    void set_desc_wfc_Eij(const int& nbasis,
        const int& nbands,
        const int& lld);
#endif

    int nspin = 1;
    int* loc_sizes;
    int loc_size;

    /**
     * @brief set row and col begin index for each atom
     * it should be called after:
     * 1. nrow and ncol are set;
     * 2. global2local_row_ and global2local_col_ are set;
     * @param iat2iwt : the map from atom index to global oribtal indexes
     * @param nat : number of atoms
     * @param nlocal : number of global orbitals
     */
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

    /**
     * @brief gather global indexes of orbitals in this processor
     * get_indexes_row() : global indexes (~NLOCAL) of rows of Hamiltonian matrix in this processor
     * get_indexes_col() : global indexes (~NLOCAL) of columns of Hamiltonian matrix in this processor
     * get_indexes_row(iat) : global indexes (~nw) of rows of Hamiltonian matrix in atom iat
     * get_indexes_col(iat) : global indexes (~nw) of columns of Hamiltonian matrix in atom iat
    */
    std::vector<int> get_indexes_row() const;
    std::vector<int> get_indexes_col() const;
    std::vector<int> get_indexes_row(int iat) const;
    std::vector<int> get_indexes_col(int iat) const;

    // private:
        // orbital index for each atom
    std::vector<int> atom_begin_row;
    std::vector<int> atom_begin_col;

    const int* iat2iwt_ = nullptr;

};
#endif