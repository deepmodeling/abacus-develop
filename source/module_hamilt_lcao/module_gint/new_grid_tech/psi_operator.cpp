#include "psi_operator.h"
#include "module_base/blas_connector.h"

namespace Gint
{

int PsiOperator::atom_pair_startidx_(int a, int b) const
{
    for(int i = 0; i < biggrid_->get_meshgrid_num(); ++i)
    {
        if(is_atom_on_mgrid_[a][i] && is_atom_on_mgrid_[b][i])
        {
            return i;
        }
    }
}

int PsiOperator::atom_pair_endidx_(int a, int b) const
{
    for(int i = biggrid_->get_meshgrid_num() - 1; i >= 0; --i)
    {
        if(is_atom_on_mgrid_[a][i] && is_atom_on_mgrid_[b][i])
        {
            return i;
        }
    }
}

PsiOperator::PsiOperator(std::shared_ptr<BigGrid> biggrid)
: biggrid_(biggrid)
{
    rows_ = biggrid_->get_biggrid_info()->get_nmxyz();
    cols_ = biggrid_->get_mgrid_psi_len();
    is_atom_on_mgrid_ = biggrid_->get_is_atom_on_mgrid();
    atom_startidx_ = biggrid_->get_atom_startidx();
    atom_psi_len_ = biggrid_->get_atom_psi_len();
}

void PsiOperator::psi_times_dm(const hamilt::HContainer<double>& DM, const double* const* psi, double** result, const bool is_symm) const
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
            dsymm_(&side, &uplo, &atom_psi_len_[i], &rows_, &alpha, dm_mat->get_pointer(), &atom_psi_len_[i],
                &psi[0][atom_startidx_[i]], &cols_, &beta, &result[0][atom_startidx_[i]], &cols_);
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

            // if len == 0, it means this atom pair does not affect any meshgrid in this biggrid
            if(len == 0)
            {
                continue;
            }

            dgemm_(&trans, &trans, &atom_psi_len_[j], &len, &atom_psi_len_[i], &alpha1, dm_mat->get_pointer(), &atom_psi_len_[j],
                &psi[start_idx][atom_startidx_[i]], &cols_, &beta, &result[start_idx][atom_startidx_[j]], &cols_);
        }
    }
}

void PsiOperator::psi_times_vldr3(const double* vldr3, const double* const* psi, double** result) const
{
    for(int i = 0; i < biggrid_->get_meshgrid_num(); i++)
    {
        for(int j = 0; j < biggrid_->get_atom_num(); j++)
        {
            for(int k = 0; k < atom_psi_len_[j]; k++)
            {
                const int idx = atom_startidx_[j];
                result[i][idx + k] = psi[i][idx + k] * vldr3[i];
            }
        }
    }
}

void PsiOperator::psi_times_psi_vldr3(
    const double* const* psi,
    const double* const* psi_vldr3,
    hamilt::HContainer<double>* hr) const
{
    const char transa='N', transb='T';
    const double alpha=1, beta=1;

    for(int i = 0; i < biggrid_->get_atom_num(); ++i)
    {
        const auto atom_i = biggrid_->get_atoms()[i];
        const auto r_i = atom_i->get_r();

        for(int j = i; j < biggrid_->get_atom_num(); ++j)
        {
            const auto atom_j = biggrid_->get_atoms()[j];
            const auto r_j = atom_j->get_r();
            // FIXME may be r = r_j - r_i
            const auto result = hr->find_matrix(atom_i->get_iat(), atom_j->get_iat(), r_i-r_j);

            if(result == nullptr)
            {
                continue;
            }

            const int start_idx = atom_pair_startidx_(i, j);
            const int end_idx = atom_pair_endidx_(i, j);
            const int len = end_idx - start_idx + 1;

            if(len == 0)
            {
                continue;
            }

            dgemm_(&transa, &transb, &atom_psi_len_[j], &atom_psi_len_[i], &len, &alpha, &psi_vldr3[start_idx][atom_startidx_[j]],
                &cols_,&psi[start_idx][atom_startidx_[i]], &cols_, &beta, result->get_pointer(), &atom_psi_len_[j]);
        }
    }
}


} // namespace Gint