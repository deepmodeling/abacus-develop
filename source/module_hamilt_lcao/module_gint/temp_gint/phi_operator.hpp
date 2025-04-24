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

// phi_dm(ir,iwt_2) = \sum_{iwt_1} phi(ir,iwt_1) * dm(iwt_1,iwt_2)
template<typename T>
void PhiOperator::phi_mul_dm(
    const T*const phi,                  // phi(ir,iwt)
    const HContainer<T>& dm,            // dm(iwt_1,iwt_2)
    const bool is_symm,
    T*const phi_dm) const               // phi_dm(ir,iwt)
{
    ModuleBase::GlobalFunc::ZEROS(phi_dm, rows_ * cols_);

    for(int i = 0; i < biggrid_->get_atoms_num(); ++i)
    {
        const auto atom_i = biggrid_->get_atom(i);
        const auto r_i = atom_i->get_R();

        if(is_symm)
        {
            const auto dm_mat = dm.find_matrix(atom_i->get_iat(), atom_i->get_iat(), 0, 0, 0);
            constexpr T alpha = 1.0;
            constexpr T beta = 1.0;
            BlasConnector::symm_cm(
                'L', 'U',
                atoms_phi_len_[i], rows_,
                alpha, dm_mat->get_pointer(), atoms_phi_len_[i],
                       &phi[0 * cols_ + atoms_startidx_[i]], cols_,
                beta, &phi_dm[0 * cols_ + atoms_startidx_[i]], cols_);
        }

        const int start = is_symm ? i + 1 : 0;

        for(int j = start; j < biggrid_->get_atoms_num(); ++j)
        {
            const auto atom_j = biggrid_->get_atom(j);
            const auto r_j = atom_j->get_R();
            // FIXME may be r = r_j - r_i
            const auto dm_mat = dm.find_matrix(atom_i->get_iat(), atom_j->get_iat(), r_i-r_j);

            // if dm_mat is nullptr, it means this atom pair does not affect any meshgrid in the unitcell
            if(dm_mat == nullptr)
            {
                continue;
            }

            const int start_idx = get_atom_pair_start_end_idx_(i, j).first;
            const int end_idx = get_atom_pair_start_end_idx_(i, j).second;
            const int len = end_idx - start_idx + 1;

            // if len<=0, it means this atom pair does not affect any meshgrid in this biggrid
            if(len <= 0)
            {
                continue;
            }

            const T alpha = is_symm ? 2.0 : 1.0;
            constexpr T beta = 1.0;
            BlasConnector::gemm(
                'N', 'N',
                len, atoms_phi_len_[j], atoms_phi_len_[i],
                alpha, &phi[start_idx * cols_ + atoms_startidx_[i]], cols_,
                       dm_mat->get_pointer(), atoms_phi_len_[j],
                beta, &phi_dm[start_idx * cols_ + atoms_startidx_[j]], cols_);
        }
    }
}

// result(ir) = phi(ir) * vl(ir)
template<typename T>
void PhiOperator::phi_mul_vldr3(
    const T*const vl,                   // vl(ir)
    const T dr3,
    const T*const phi,                  // phi(ir,iwt)
    T*const result) const               // result(ir,iwt)
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

// hr(iwt_i,iwt_j) = \sum_{ir} phi_i(ir,iwt_i) * phi_i(ir,iwt_j)
// this is a thread-safe function
template<typename T>
void PhiOperator::phi_mul_phi(
    const T*const phi_i,                // phi_i(ir,iwt)
    const T*const phi_j,                // phi_j(ir,iwt)
    HContainer<T>& hr,                  // hr(iwt_i,iwt_j)
    const Triangular_Matrix triangular_matrix) const
{
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

            constexpr T alpha=1, beta=1;
            BlasConnector::gemm(
                'T', 'N', n_i, n_j, len,
		        alpha, phi_i + start_idx * cols_ + atoms_startidx_[i], cols_,
                       phi_j + start_idx * cols_ + atoms_startidx_[j], cols_,
		        beta, tmp_hr.data(), n_j,
                base_device::AbacusDevice_t::CpuDevice);

            result->add_array_ts(tmp_hr.data());
        }
    }
}

// rho(ir) = \sum_{iwt} \phi_i(ir,iwt) * \phi_j^*(ir,iwt)
template<typename T>
void PhiOperator::phi_dot_phi(
    const T*const phi_i,           // phi_i(ir,iwt)
    const T*const phi_j,           // phi_j(ir,iwt)
    T*const rho) const             // rho(ir)
{
    constexpr int inc = 1;
    for(int i = 0; i < biggrid_->get_mgrids_num(); ++i)
    {
        rho[meshgrids_local_idx_[i]] += BlasConnector::dotc(cols_, phi_j+i*cols_, inc, phi_i+i*cols_, inc);
    }
}

} // namespace ModuleGint