#pragma once

#include "phi_operator.h"
#include "module_base/blas_connector.h"
#include "module_base/global_function.h"

namespace ModuleGint
{

template<typename T>
void PhiOperator::set_phi(T* phi) const
{
    for(int i = 0; i < biggrid_->get_atoms_num(); ++i)
    {
        const auto atom = biggrid_->get_atom(i);
        atom->set_phi(atoms_relative_coords_[i], cols_, phi);
        phi += atom->get_nw();
    }
}

template<typename T>
void PhiOperator::phi_mul_vldr3(
    const T*const vl,
    const T dr3,
    const T*const phi,
    T*const result) const
{
    int idx = 0;
    for(int i = 0; i < biggrid_->get_mgrids_num(); i++)
    {
        T vldr3_mgrid = vl[meshgrids_local_idx_[i]] * dr3;
        for(int j = 0; j < cols_; j++)
        {
            result[idx] = phi[idx] * vldr3_mgrid;
            idx++;
        }
    }
}

// this is a thread-safe function
template<typename T>
void PhiOperator::phi_mul_phi(
    const T*const phi_i,                // phi_i(igrid,iwt_i)
    const T*const phi_j,                // phi_j(igrid,iwt_j)
    HContainer<T>& hr,                  // hr(iwt_i,iwt_j)
    const Triangular_Matrix triangular_matrix) const
{
    const char transa='T', transb='N';
    const T alpha=1, beta=1;

    std::vector<T> tmp_hr;
    for(int i = 0; i < biggrid_->get_atoms_num(); ++i)
    {
        const auto atom_i = biggrid_->get_atom(i);
        const auto& r_i = atom_i->get_R();
        const int iat_i = atom_i->get_iat();
        const int n_i = atoms_phi_len_[i];

        for(int j = 0; j < biggrid_->get_atoms_num(); ++j)
        {
            const auto atom_j = biggrid_->get_atom(j);
            const auto& r_j = atom_j->get_R();
            const int iat_j = atom_j->get_iat();
            const int n_j = atoms_phi_len_[j];

            // only calculate the upper triangle matrix
            if(triangular_matrix==Triangular_Matrix::Upper && iat_i>iat_j)
            {
                continue;
            }
            // only calculate the upper triangle matrix
            else if(triangular_matrix==Triangular_Matrix::Lower && iat_i<iat_j)
            {
                continue;
            }

            // FIXME may be r = r_j - r_i
            const auto result = hr.find_matrix(iat_i, iat_j, r_i-r_j);

            if(result == nullptr)
            {
                continue;
            }

            const int start_idx = get_atom_pair_start_end_idx_(i, j).first;
            const int end_idx = get_atom_pair_start_end_idx_(i, j).second;
            const int len = end_idx - start_idx + 1;

            if(len <= 0)
            {
                continue;
            }

            tmp_hr.resize(n_i * n_j);
            ModuleBase::GlobalFunc::ZEROS(tmp_hr.data(), n_i*n_j);

            BlasConnector::gemm(
                transa, transb, n_i, n_j, len,
		        alpha, phi_i + start_idx * cols_ + atoms_startidx_[i], cols_,
                       phi_j + start_idx * cols_ + atoms_startidx_[j], cols_,
		        beta, tmp_hr.data(), n_j,
                base_device::AbacusDevice_t::CpuDevice);

            result->add_array_ts(tmp_hr.data());
        }
    }
}

} // namespace ModuleGint