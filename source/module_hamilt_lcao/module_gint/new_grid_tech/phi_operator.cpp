#include "phi_operator.h"
#include "module_base/blas_connector.h"
#include "module_base/matrix.h"

namespace ModuleGint
{

void PhiOperator::set_bgrid(std::shared_ptr<const BigGrid> biggrid)
{
    biggrid_ = biggrid;
    rows_ = biggrid_->get_meshgrid_num();
    cols_ = biggrid_->get_mgrid_phi_len();

    biggrid_->set_atoms_startidx(atoms_startidx_);
    biggrid_->set_atoms_phi_len(atoms_phi_len_);
    biggrid_->set_mgrids_local_idx(meshgrids_local_idx_);

    int atoms_num = biggrid_->get_atom_num();
    atoms_relative_coords_.resize(atoms_num);
    is_atom_on_mgrid_.resize(atoms_num);
    for(int i = 0; i < atoms_num; ++i)
    {
        biggrid_->set_atom_relative_coords(biggrid_->get_atom(i), atoms_relative_coords_[i]);
        is_atom_on_mgrid_[i].resize(rows_);
        for(int j = 0; j < rows_; ++j)
        {
            is_atom_on_mgrid_[i][j] = atoms_relative_coords_[i][j].norm() <= biggrid_->get_atom(i)->get_rcut();
        }
    }
}

void PhiOperator::set_phi(double* phi) const
{
    for(int i = 0; i < biggrid_->get_atom_num(); ++i)
    {
        const auto atom = biggrid_->get_atom(i);
        atom->set_phi(atoms_relative_coords_[i], cols_, phi);
        phi += atom->get_nw();
    }
}

void PhiOperator::set_phi_dphi(double* phi, double* dphi_x, double* dphi_y, double* dphi_z) const
{
    for(int i = 0; i < biggrid_->get_atom_num(); ++i)
    {
        const auto atom = biggrid_->get_atom(i);
        atom->set_phi_dphi(atoms_relative_coords_[i], cols_, phi, dphi_x, dphi_y, dphi_z);
        phi += atom->get_nw();
        dphi_x += atom->get_nw();
        dphi_y += atom->get_nw();
        dphi_z += atom->get_nw();
    }
}

void PhiOperator::set_ddphi(
    double* ddphi_xx, double* ddphi_xy, double* ddphi_xz,
    double* ddphi_yy, double* ddphi_yz, double* ddphi_zz) const
{
    for(int i = 0; i < biggrid_->get_atom_num(); ++i)
    {
        const auto atom = biggrid_->get_atom(i);
        atom->set_ddphi(atoms_relative_coords_[i], cols_, ddphi_xx, ddphi_xy, ddphi_xz, ddphi_yy, ddphi_yz, ddphi_zz);
        ddphi_xx += atom->get_nw();
        ddphi_xy += atom->get_nw();
        ddphi_xz += atom->get_nw();
        ddphi_yy += atom->get_nw();
        ddphi_yz += atom->get_nw();
        ddphi_zz += atom->get_nw();
    }
}

void PhiOperator::phi_mul_dm(
    const double* const* phi, 
    const hamilt::HContainer<double>& DM, 
    const bool is_symm, double** phi_dm) const
{
    // parameters for lapack subroutines
    constexpr char side = 'L';
    constexpr char uplo = 'U';
    const char trans = 'N';
    const double alpha = 1.0;
    const double beta = 1.0;
    const double alpha1 = is_symm ? 2.0 : 1.0;

    for(int i = 0; i < biggrid_->get_atom_num(); ++i)
    {
        const auto atom_i = biggrid_->get_atom(i);
        const auto r_i = atom_i->get_R();

        if(is_symm)
        {
            const auto dm_mat = DM.find_matrix(atom_i->get_iat(), atom_i->get_iat(), 0, 0, 0);
            dsymm_(&side, &uplo, &atoms_phi_len_[i], &rows_, &alpha, dm_mat->get_pointer(), &atoms_phi_len_[i],
                &phi[0][atoms_startidx_[i]], &cols_, &beta, &phi_dm[0][atoms_startidx_[i]], &cols_);
        }

        const int start = is_symm ? i + 1 : 0;

        for(int j = start; j < biggrid_->get_atom_num(); ++j)
        {
            const auto atom_j = biggrid_->get_atom(j);
            const auto r_j = atom_j->get_R();
            // FIXME may be r = r_j - r_i
            const auto dm_mat = DM.find_matrix(atom_i->get_iat(), atom_j->get_iat(), r_i-r_j);

            // if dm_mat is nullptr, it means this atom pair does not affect any meshgrid in the unitcell
            if(dm_mat == nullptr)
            {
                continue;
            }

            const int start_idx = atom_pair_startidx_(i, j);
            const int end_idx = atom_pair_endidx_(i, j);
            const int len = end_idx - start_idx + 1;

            // if len<=0, it means this atom pair does not affect any meshgrid in this biggrid
            if(len <= 0)
            {
                continue;
            }

            dgemm_(&trans, &trans, &atoms_phi_len_[j], &len, &atoms_phi_len_[i], &alpha1, dm_mat->get_pointer(), &atoms_phi_len_[j],
                &phi[start_idx][atoms_startidx_[i]], &cols_, &beta, &phi_dm[start_idx][atoms_startidx_[j]], &cols_);
        }
    }
}

void PhiOperator::phi_mul_vldr3(const double* vl, const double dr3, const double* const* phi, double** result) const
{
    for(int i = 0; i < biggrid_->get_meshgrid_num(); i++)
    {
        double vldr3_mgrid = vl[meshgrids_local_idx_[i]] * dr3;
        for(int j = 0; j < cols_; j++)
        {
            result[i][j] = phi[i][j] * vldr3_mgrid;
        }
    }
}

void PhiOperator::phi_mul_phi_vldr3(
    const double* const* phi,
    const double* const* phi_vldr3,
    hamilt::HContainer<double>* hr) const
{
    const char transa='N', transb='T';
    const double alpha=1, beta=1;

    for(int i = 0; i < biggrid_->get_atom_num(); ++i)
    {
        const auto atom_i = biggrid_->get_atom(i);
        const auto& r_i = atom_i->get_R();
        const int iat_i = atom_i->get_iat();

        for(int j = 0; j < biggrid_->get_atom_num(); ++j)
        {
            const auto atom_j = biggrid_->get_atom(j);
            const auto& r_j = atom_j->get_R();
            const int iat_j = atom_j->get_iat();

            // only calculate the upper triangle matrix
            if(iat_i > iat_j)
            {
                continue;
            }

            // FIXME may be r = r_j - r_i
            const auto result = hr->find_matrix(iat_i, iat_j, r_i-r_j);

            if(result == nullptr)
            {
                continue;
            }

            const int start_idx = atom_pair_startidx_(i, j);
            const int end_idx = atom_pair_endidx_(i, j);
            const int len = end_idx - start_idx + 1;

            if(len <= 0)
            {
                continue;
            }

            dgemm_(&transa, &transb, &atoms_phi_len_[j], &atoms_phi_len_[i], &len, &alpha, &phi_vldr3[start_idx][atoms_startidx_[j]],
                &cols_,&phi[start_idx][atoms_startidx_[i]], &cols_, &beta, result->get_pointer(), &atoms_phi_len_[j]);
        }
    }
}

void PhiOperator::phi_dot_phi_dm(
    const double* const* phi,
    const double* const* phi_dm,
    double* rho) const
{
    const int inc = 1;
    for(int i = 0; i < biggrid_->get_meshgrid_num(); ++i)
    {
        rho[meshgrids_local_idx_[i]] += ddot_(&cols_, phi[i], &inc, phi_dm[i], &inc);
    }
}

void PhiOperator::phi_dot_dphi(
    const double* const* phi,
    const double* const* dphi_x,
    const double* const* dphi_y,
    const double* const* dphi_z,
    ModuleBase::matrix *fvl) const
{
    for(int i = 0; i < biggrid_->get_atom_num(); ++i)
    {
        const int start_idx = atoms_startidx_[i];
        const int phi_len = atoms_phi_len_[i];
        double rx = 0, ry = 0, rz = 0;
        for(int j = 0; j < biggrid_->get_meshgrid_num(); ++j)
        {
            for(int k = 0; k < phi_len; ++k)
            {
                const double phi_val = phi[j][start_idx + k];
                rx += phi_val * dphi_x[j][start_idx + k];
                ry += phi_val * dphi_y[j][start_idx + k];
                rz += phi_val * dphi_z[j][start_idx + k];
            }
        }
        fvl[0](i, 0) += rx * 2;
        fvl[0](i, 1) += ry * 2;
        fvl[0](i, 2) += rz * 2;
    }
}

void PhiOperator::phi_dot_dphi_r(
    const double* const *phi,
    const double* const *dphi_x,
    const double* const *dphi_y,
    const double* const *dphi_z,
    ModuleBase::matrix *svl) const
{
    double sxx = 0, sxy = 0, sxz = 0, syy = 0, syz = 0, szz = 0;
    for(int i = 0; i < biggrid_->get_meshgrid_num(); ++i)
    {
        for(int j = 0; j < biggrid_->get_atom_num(); ++j)
        {
            const int start_idx = atoms_startidx_[j];
            for(int k = 0; k < atoms_phi_len_[j]; ++k)
            {
                const int col_idx = start_idx + k;
                const double phi_val = phi[i][col_idx];
                const Vec3d& r3 = atoms_relative_coords_[j][i];
                sxx += phi_val * dphi_x[i][col_idx] * r3[0];
                sxy += phi_val * dphi_x[i][col_idx] * r3[1];
                sxz += phi_val * dphi_x[i][col_idx] * r3[2];
                syy += phi_val * dphi_y[i][col_idx] * r3[1];
                syz += phi_val * dphi_y[i][col_idx] * r3[2];
                szz += phi_val * dphi_z[i][col_idx] * r3[2];
            }
        }
    }
    svl[0](0, 0) += sxx * 2;
    svl[0](0, 1) += sxy * 2;
    svl[0](0, 2) += sxz * 2;
    svl[0](1, 1) += syy * 2;
    svl[0](1, 2) += syz * 2;
    svl[0](2, 2) += szz * 2;
}


//===============================
// private methods
//===============================
int PhiOperator::atom_pair_startidx_(int a, int b) const
{
    for(int i = 0; i < biggrid_->get_meshgrid_num(); ++i)
    {
        if(is_atom_on_mgrid_[a][i] && is_atom_on_mgrid_[b][i])
        {
            return i;
        }
    }
    return biggrid_->get_meshgrid_num();
}

int PhiOperator::atom_pair_endidx_(int a, int b) const
{
    for(int i = biggrid_->get_meshgrid_num() - 1; i >= 0; --i)
    {
        if(is_atom_on_mgrid_[a][i] && is_atom_on_mgrid_[b][i])
        {
            return i;
        }
    }
    return -1;
}

} // namespace ModuleGint