#include "phi_operator.h"
#include "module_base/blas_connector.h"

namespace ModuleGint
{

PhiOperator::PhiOperator(std::shared_ptr<BigGrid> biggrid)
: biggrid_(biggrid)
{
    rows_ = biggrid_->get_biggrid_info()->get_nmxyz();
    cols_ = biggrid_->get_mgrid_phi_len();
    is_atom_on_mgrids_ = biggrid_->get_is_atom_on_mgrids();
    meshgrids_local_idx_ = biggrid_->get_mgrids_local_idx();
    atom_startidx_ = biggrid_->get_atom_startidx();
    atom_phi_len_ = biggrid_->get_atom_phi_len();
}

void PhiOperator::set_phi(double* phi) const
{
    for(const auto& atom : biggrid_->get_atoms())
    {
        atom->set_phi(biggrid_->get_atom_relative_coords(*atom), cols_, phi);
        phi += atom->get_nw();
    }
}

void PhiOperator::set_phi_dphi(double* phi, double* dphi_x, double* dphi_y, double* dphi_z) const
{
    for(const auto& atom : biggrid_->get_atoms())
    {
        atom->set_phi_dphi(biggrid_->get_atom_relative_coords(*atom), cols_, phi, dphi_x, dphi_y, dphi_z);
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
    for(const auto& atom : biggrid_->get_atoms())
    {
        atom->set_ddphi(biggrid_->get_atom_relative_coords(*atom), cols_, ddphi_xx, ddphi_xy, ddphi_xz, ddphi_yy, ddphi_yz, ddphi_zz);
        ddphi_xx += atom->get_nw();
        ddphi_xy += atom->get_nw();
        ddphi_xz += atom->get_nw();
        ddphi_yy += atom->get_nw();
        ddphi_yz += atom->get_nw();
        ddphi_zz += atom->get_nw();
    }
}

void PhiOperator::phi_mul_dm(const hamilt::HContainer<double>& DM, const double* const* phi, double** result, const bool is_symm) const
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
        const auto atom_i = biggrid_->get_atoms()[i];
        const auto r_i = atom_i->get_r();

        if(is_symm)
        {
            const auto dm_mat = DM.find_matrix(atom_i->get_iat(), atom_i->get_iat(), 0, 0, 0);
            dsymm_(&side, &uplo, &atom_phi_len_[i], &rows_, &alpha, dm_mat->get_pointer(), &atom_phi_len_[i],
                &phi[0][atom_startidx_[i]], &cols_, &beta, &result[0][atom_startidx_[i]], &cols_);
        }

        const int start = is_symm ? i + 1 : 0;

        for(int j = start; j < biggrid_->get_atom_num(); ++j)
        {
            const auto atom_j = biggrid_->get_atoms()[j];
            const auto r_j = atom_j->get_r();
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

            dgemm_(&trans, &trans, &atom_phi_len_[j], &len, &atom_phi_len_[i], &alpha1, dm_mat->get_pointer(), &atom_phi_len_[j],
                &phi[start_idx][atom_startidx_[i]], &cols_, &beta, &result[start_idx][atom_startidx_[j]], &cols_);
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
        const auto atom_i = biggrid_->get_atoms()[i];
        const auto& r_i = atom_i->get_r();
        const int iat_i = atom_i->get_iat();

        for(int j = 0; j < biggrid_->get_atom_num(); ++j)
        {
            const auto atom_j = biggrid_->get_atoms()[j];
            const auto& r_j = atom_j->get_r();
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

            dgemm_(&transa, &transb, &atom_phi_len_[j], &atom_phi_len_[i], &len, &alpha, &phi_vldr3[start_idx][atom_startidx_[j]],
                &cols_,&phi[start_idx][atom_startidx_[i]], &cols_, &beta, result->get_pointer(), &atom_phi_len_[j]);
        }
    }
}

//===============================
// private methods
//===============================
int PhiOperator::atom_pair_startidx_(int a, int b) const
{
    for(int i = 0; i < biggrid_->get_meshgrid_num(); ++i)
    {
        if(is_atom_on_mgrids_[a][i] && is_atom_on_mgrids_[b][i])
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
        if(is_atom_on_mgrids_[a][i] && is_atom_on_mgrids_[b][i])
        {
            return i;
        }
    }
    return -1;
}

} // namespace ModuleGint