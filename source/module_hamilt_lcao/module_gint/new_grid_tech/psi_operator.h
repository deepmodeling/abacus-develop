#pragma once

#include <memory>
#include <vector>
#include <module_hamilt_lcao/module_hcontainer/hcontainer.h>
#include "biggrid.h"

namespace Gint
{

/**
 * @brief The class PsiOperator is used to perform operations on the wave function matrix psi, dpsi, etc.
 * 
 * In fact, the variables and functions of this class could be placed in the BigGrid class, but the lifecycle of the BigGrid class is relatively long. 
 * We do not want the BigGrid to contain too many member variables, as this could lead to excessive memory usage. 
 * Therefore, we separate this class out, so it can be destroyed after use.
 */
class PsiOperator
{
    public:
    // constructor
    PsiOperator(std::shared_ptr<BigGrid> biggrid);

    void psi_times_dm(
        const hamilt::HContainer<double>& DM, 
        const double* const* psi, 
        double** result,
        const bool is_symm) const;

    void psi_times_vldr3(
        const double* vldr3,
        const double* const* psi,
        double** result) const;
    
    void psi_times_psi_vldr3(
        const double* const* psi,
        const double* const* psi_vldr3,
        hamilt::HContainer<double>* hr) const;
    

    private:

    // get the index of the first meshgrid that both atom a and atom b affect
    int atom_pair_startidx_(int a, int b) const;

    // get the index of the last meshgrid that both atom a and atom b affect
    int atom_pair_endidx_(int a, int b) const;

    // the row number of the psi matrix
    // rows_ = biggrid_->get_meshgrid_num()
    int rows_;
    
    // the column number of the psi matrix
    // cols_ = biggrid_->get_mgrid_psi_len()
    int cols_;

    // the big grid that the psi matrix is associated with
    std::shared_ptr<const BigGrid> biggrid_;

    // record whether the atom affects the meshgrid
    // is_atom_on_mgrid_[i][j] = true if the ith atom affects the jth meshgrid, otherwise false
    std::vector<std::vector<bool>> is_atom_on_mgrid_;

    // the start index of the psi of each atom
    std::vector<int> atom_startidx_;

    // the length of psi of each atom
    // atom_psi_len_[i] = biggrid_->get_atoms()[i]->get_nw()
    std::vector<int> atom_psi_len_;
}

}