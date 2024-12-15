#pragma once

#include <memory>
#include <vector>
#include <module_hamilt_lcao/module_hcontainer/hcontainer.h>
#include "biggrid.h"

namespace Gint
{

/**
 * @brief The class phiOperator is used to perform operations on the wave function matrix phi, dphi, etc.
 * 
 * In fact, the variables and functions of this class could be placed in the BigGrid class, but the lifecycle of the BigGrid class is relatively long. 
 * We do not want the BigGrid to contain too many member variables, as this could lead to excessive memory usage. 
 * Therefore, we separate this class out, so it can be destroyed after use.
 */
class PhiOperator
{
    public:
    // constructor
    PhiOperator(std::shared_ptr<BigGrid> biggrid);

    // getter
    int get_rows() const {return rows_;};
    int get_cols() const {return cols_;};

    // get phi of the big grid
    // the dimension of phi is num_mgrids * (\sum_{i=0}^{atoms_->size()} atoms_[i]->nw)
    void set_phi(double* phi) const;

    // get phi and the gradient of phi of the big grid
    // the dimension of phi and dphi is num_mgrids * (\sum_{i=0}^{atoms_->size()} atoms_[i]->nw)
    void set_phi_dphi(double* phi, double* dphi_x, double* dphi_y, double* dphi_z) const;

    // get the hessian of the wave function values of the big grid
    // the dimension of ddphi is num_mgrids * (\sum_{i=0}^{atoms_->size()} atoms_[i]->nw)
    void set_ddphi(
        double* ddphi_xx, double* ddphi_xy, double* ddphi_xz,
        double* ddphi_yy, double* ddphi_yz, double* ddphi_zz) const;

    void phi_mul_dm(
        const hamilt::HContainer<double>& DM, 
        const double* const* phi, 
        double** result,
        const bool is_symm) const;

    void phi_mul_vldr3(
        const double* vldr3,
        const double* const* phi,
        double** result) const;
    
    void phi_mul_phi_vldr3(
        const double* const* phi,
        const double* const* phi_vldr3,
        hamilt::HContainer<double>* hr) const;
    

    private:

    // get the index of the first meshgrid that both atom a and atom b affect
    int atom_pair_startidx_(int a, int b) const;

    // get the index of the last meshgrid that both atom a and atom b affect
    int atom_pair_endidx_(int a, int b) const;

    // the row number of the phi matrix
    // rows_ = biggrid_->get_meshgrid_num()
    int rows_;
    
    // the column number of the phi matrix
    // cols_ = biggrid_->get_mgrid_phi_len()
    int cols_;

    // the coordinates of the meshgrids
    std::vector<Vec3d> meshgrid_coords_;

    // the big grid that the phi matrix is associated with
    std::shared_ptr<const BigGrid> biggrid_;

    // record whether the atom affects the meshgrid
    // is_atom_on_mgrids_[i][j] = true if the ith atom affects the jth meshgrid, otherwise false
    std::vector<std::vector<bool>> is_atom_on_mgrids_;

    // the start index of the phi of each atom
    std::vector<int> atom_startidx_;

    // the length of phi of each atom
    // atom_phi_len_[i] = biggrid_->get_atoms()[i]->get_nw()
    std::vector<int> atom_phi_len_;
};

}