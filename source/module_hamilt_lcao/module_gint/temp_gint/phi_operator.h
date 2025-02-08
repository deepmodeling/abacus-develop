#pragma once

#include <memory>
#include <vector>
#include <utility>
#include <module_hamilt_lcao/module_hcontainer/hcontainer.h>
#include "big_grid.h"

namespace ModuleGint
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
    PhiOperator()=default;

    // set the big grid that the phiOperator is associated with
    void set_bgrid(std::shared_ptr<const BigGrid> biggrid);

    // getter
    int get_rows() const {return rows_;};
    int get_cols() const {return cols_;};

    // get phi of the big grid
    // the dimension of phi is num_mgrids * (\sum_{i=0}^{atoms_->size()} atoms_[i]->nw)
    void set_phi(double* phi) const;

    // get phi and the gradient of phi of the big grid
    // the dimension of phi and dphi is num_mgrids * (\sum_{i=0}^{atoms_->size()} atoms_[i]->nw)
    // if you do not need phi, you can set phi to nullptr.
    void set_phi_dphi(double* phi, double* dphi_x, double* dphi_y, double* dphi_z) const;

    // get the hessian of the wave function values of the big grid
    // the dimension of ddphi is num_mgrids * (\sum_{i=0}^{atoms_->size()} atoms_[i]->nw)
    void set_ddphi(
        double* ddphi_xx, double* ddphi_xy, double* ddphi_xz,
        double* ddphi_yy, double* ddphi_yz, double* ddphi_zz) const;

    void phi_mul_dm(
        const double* phi, 
        const HContainer<double>& dm, 
        const bool is_symm, double* phi_dm) const;

    void phi_mul_vldr3(
        const double* vl,
        const double dr3,
        const double* phi,
        double* result) const;
    
    void phi_mul_phi_vldr3(
        const double* phi,
        const double* phi_vldr3,
        HContainer<double>* hr) const;
    
    void phi_dot_phi_dm(
        const double* phi,
        const double* phi_dm,
        double* rho) const;

    void phi_dot_dphi(
        const double* phi,
        const double* dphi_x,
        const double* dphi_y,
        const double* dphi_z,
        ModuleBase::matrix *fvl) const;
    
    void phi_dot_dphi_r(
        const double* phi,
        const double* dphi_x,
        const double* dphi_y,
        const double* dphi_z,
        ModuleBase::matrix *svl) const;

    private:
    void init_atom_pair_start_end_idx_();

    // get the index of the first and the last meshgrid that both atom a and atom b affect
    // Note that atom_pair_start_end_idx_ only stores the cases where a <= b, so this function is needed to retrieve the value
    const std::pair<int, int>& get_atom_pair_start_end_idx_(int a, int b) const
    {
        int x = std::min(a, b);
        int y = std::abs(a - b);
        return atom_pair_start_end_idx_[(2 * biggrid_->get_atoms_num() - x + 1) * x / 2 + y];
    };

    // the row number of the phi matrix
    // rows_ = biggrid_->get_mgrids_num()
    int rows_;
    
    // the column number of the phi matrix
    // cols_ = biggrid_->get_mgrid_phi_len()
    int cols_;
    
    // the local index of the meshgrids
    std::vector<int> meshgrids_local_idx_;

    // the big grid that the phi matrix is associated with
    std::shared_ptr<const BigGrid> biggrid_;

    // the relative coordinates of the atoms and the meshgrids
    // atoms_relative_coords_[i][j] is the relative coordinate of the jth meshgrid and the ith atom
    std::vector<std::vector<Vec3d>> atoms_relative_coords_;

    // record whether the atom affects the meshgrid
    // is_atom_on_mgrid_[i][j] = true if the ith atom affects the jth meshgrid, otherwise false
    // FIXME,std::vector<std::vector<bool>> is not a efficient data structure, we can use a 1D array to replace it. 
    std::vector<std::vector<bool>> is_atom_on_mgrid_;

    // the start index of the phi of each atom
    std::vector<int> atoms_startidx_;

    // the length of phi of each atom
    // atoms_phi_len_[i] = biggrid_->get_atom(i)->get_nw()
    // TODO: remove it
    std::vector<int> atoms_phi_len_;

    // This data structure is used to store the index of the first and last meshgrid affected by each atom pair
    std::vector<std::pair<int, int>> atom_pair_start_end_idx_; 
};

}