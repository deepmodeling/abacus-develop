#ifndef _PARALLEL_ORBITALS_H_
#define _PARALLEL_ORBITALS_H_
#include "source_base/block_cyclic_distribution.h"
#include <fstream>
#include <vector>

/// This class packs the information of 2D-block-cyclic for LCAO code:
/// parallel distribution of basis, wavefunction and matrix.
///
/// R4 (single storage): this class must never add a copy of any BlockCyclicDistribution
/// state field (nrow/ncol/nloc/nb/dim0/dim1/desc/blacs_ctxt).  The only
/// authoritative copy lives in the para_2d_ member below; re-deriving the
/// distribution from a duplicated field would silently desynchronize the
/// LCAO band/atom layers that are computed from it.
class Parallel_Orbitals
{
public:
    Parallel_Orbitals();
    ~Parallel_Orbitals();

    /// Access to the block-cyclic distribution information.  The non-const
    /// overload exists only for the two cases a const getter cannot serve:
    /// non-const mutation by setup_2d_division, and test-side descriptor writes.
    BlockCyclicDistribution& para_2d() { return para_2d_; }

    const BlockCyclicDistribution& para_2d() const { return para_2d_; }

    /// Block-cyclic distribution of the wavefunction (nbasis x nbands).  This is
    /// a distinct layout from the square Hamiltonian distribution, shared with
    /// the same process grid; it is a temporary adapter seam, not a duplicate of
    /// para_2d_.
    const BlockCyclicDistribution& wfc_layout() const { return wfc_layout_; }

    /// Block-cyclic layout of the MO representation: S and H represented in the MOs
    /// defined by the DFT reference determinant, i.e. Psi^dagger * S * Psi and
    /// Psi^dagger * H * Psi.  The square global extent of this layout follows from
    /// the number nmo of MOs.  Shares the process grid of para_2d_.
    const BlockCyclicDistribution& mo_repr_layout() const { return mo_repr_layout_; }

    /// local size of bands, used for 2d wavefunction
    /// must divided on dim1 because of elpa interface
    int ncol_bands;
    int nrow_bands;
    
    /// ncol_bands*nrow
    long nloc_wfc;

    int lastband_in_proc;
    int lastband_number; 

    ///-------------------------------------
    /// number of elements(basis-pairs) in this processon
    /// on all adjacent atoms-pairs(2D division)
    ///-------------------------------------
    int nnr=1;
    std::vector<int> nlocdim;
    std::vector<int> nlocstart;
    
#ifdef __MPI
    int desc_wfc[9]; //for wfc, nlocal*nbands

    /// set the local sizes of the wavefunction and the MO representation
    int set_wfc_and_mo_sizes(const int& nmo/**< number of MOs*/,
        std::ofstream& ofs_running,
        std::ofstream& ofs_warning);

    ///@brief set the 2D-block-cyclic layouts of the wavefunction and the MO representation
    void set_wfc_and_mo_layouts(const int& nbasis,
        const int& nmo);
#endif

    int get_wfc_global_nbands () const;
    int get_wfc_global_nbasis () const;

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
     * ncloc() : total number of columns of Hamiltonian matrix in this processor
     * nrloc() : total number of rows of Hamiltonian matrix in this processor
     * get_ncol_atom(iat) : number of columns of Hamiltonian matrix in atom iat
     * get_nrow_atom(iat) : number of rows of Hamiltonian matrix in atom iat
    */
    int get_ncol_atom(int iat) const;
    int get_nrow_atom(int iat) const;

    bool is_invalid_atom_pair(int iat1, int iat2) const;

    int get_nbands() const;

    int nbands;

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

private:
    /// Block-cyclic distribution of the Hamiltonian matrix.  Single source of
    /// truth for nrow/ncol/nloc/nb/dim0/dim1/desc/blacs_ctxt (R4).
    BlockCyclicDistribution para_2d_;

    /// Block-cyclic distribution of the wavefunction (nbasis x nbands), sharing
    /// the process grid of para_2d_.  Temporary adapter seam for runtime
    /// consumers; not a duplicate of the Hamiltonian distribution.
    BlockCyclicDistribution wfc_layout_;

    /// Block-cyclic layout of S and H represented in the MOs defined by the DFT
    /// reference determinant.  Its square global extent follows from the number nmo
    /// of MOs; it shares the process grid of para_2d_.
    BlockCyclicDistribution mo_repr_layout_;

};
#endif
