#pragma once

#include "phi_operator.h"
#include "source_base/module_external/blas_connector.h"
#include "source_base/global_function.h"

namespace ModuleGint
{

template<typename T>
void PhiOperator::set_phi(T* phi) const
{
    for(int i = 0; i < biggrid_->get_atoms_num(); ++i)
    {
        const auto atom = biggrid_->get_atom(i);
        atom->set_phi(atom_rcoords_[i], cols_, phi);
        phi += atom->get_nw();
    }
}

// Helpers to dispatch a Tin-typed BLAS-GEMM target buffer:
//  - For Tin=double, write the GEMM result directly into the caller's
//    double* phi_dm (no scratch, no cast).
//  - For Tin=float, allocate a fp32 scratch buffer and cast the result into
//    phi_dm at the end. K of each individual GEMM is small (~atom nw) so
//    accumulating in fp32 inside the scratch is fine.
inline double* phi_mul_dm_scratch_(double* phi_dm, std::vector<double>& /*scratch*/, int /*size*/)
{
    return phi_dm;
}

template<typename Tin>
inline Tin* phi_mul_dm_scratch_(double* /*phi_dm*/, std::vector<Tin>& scratch, int size)
{
    scratch.assign(size, Tin(0));
    return scratch.data();
}

inline void phi_mul_dm_finalize_(double* /*phi_dm*/, const std::vector<double>& /*scratch*/, int /*size*/) {}

template<typename Tin>
inline void phi_mul_dm_finalize_(double* phi_dm, const std::vector<Tin>& scratch, int size)
{
    for (int k = 0; k < size; ++k)
    {
        phi_dm[k] = static_cast<double>(scratch[k]);
    }
}

// phi_dm(ir,iwt_2) = \sum_{iwt_1} phi(ir,iwt_1) * dm(iwt_1,iwt_2)
template<typename Tin>
void PhiOperator::phi_mul_dm(
    const Tin*const phi,                // phi(ir,iwt)
    const HContainer<Tin>& dm,          // dm(iwt_1,iwt_2)
    const bool is_symm,
    double*const phi_dm) const          // phi_dm(ir,iwt)
{
    std::vector<Tin> scratch;
    Tin* target = phi_mul_dm_scratch_(phi_dm, scratch, rows_ * cols_);
    ModuleBase::GlobalFunc::ZEROS(target, rows_ * cols_);

    for(int i = 0; i < biggrid_->get_atoms_num(); ++i)
    {
        const auto atom_i = biggrid_->get_atom(i);
        const auto r_i = atom_i->get_R();

        if(is_symm)
        {
            const auto dm_mat = dm.find_matrix(atom_i->get_iat(), atom_i->get_iat(), 0, 0, 0);
            constexpr Tin alpha = 1.0;
            constexpr Tin beta = 1.0;
            BlasConnector::symm_cm(
                'L', 'U',
                atoms_phi_len_[i], rows_,
                alpha, dm_mat->get_pointer(), atoms_phi_len_[i],
                       &phi[0 * cols_ + atoms_startidx_[i]], cols_,
                beta, &target[0 * cols_ + atoms_startidx_[i]], cols_);
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

            const Tin alpha = is_symm ? 2.0 : 1.0;
            constexpr Tin beta = 1.0;
            BlasConnector::gemm(
                'N', 'N',
                len, atoms_phi_len_[j], atoms_phi_len_[i],
                alpha, &phi[start_idx * cols_ + atoms_startidx_[i]], cols_,
                       dm_mat->get_pointer(), atoms_phi_len_[j],
                beta, &target[start_idx * cols_ + atoms_startidx_[j]], cols_);
        }
    }

    phi_mul_dm_finalize_(phi_dm, scratch, rows_ * cols_);
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
        T vldr3_mgrid = vl[mgrid_lidx_[i]] * dr3;
        for(int j = 0; j < cols_; j++)
        {
            result[idx] = phi[idx] * vldr3_mgrid;
            idx++;
        }
    }
}

// hr(iwt_i,iwt_j) += \sum_{ir} phi_i(ir,iwt_i) * phi_i(ir,iwt_j)
// this is a thread-safe function.
// The per-biggrid GEMM accumulator tmp_hr stays in Tin (small K, no significant
// precision loss); only the global add into HContainer<double> is widened.
template<typename Tin>
void PhiOperator::phi_mul_phi(
    const Tin*const phi_i,              // phi_i(ir,iwt)
    const Tin*const phi_j,              // phi_j(ir,iwt)
    HContainer<double>& hr,             // hr(iwt_i,iwt_j)
    const TriPart part) const
{
    std::vector<Tin> tmp_hr;
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
            if(part==TriPart::Upper && iat_i>iat_j)
            {
                continue;
            }
            // only calculate the upper triangle matrix
            else if(part==TriPart::Lower && iat_i<iat_j)
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

            constexpr Tin alpha=1, beta=1;
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

// Mixed-precision dotc wrapper. Accepts (double, double) or (double, float);
// when y is fp32 it is upcast into the caller-provided fp64 scratch buffer
// before dispatching to BlasConnector::dotc (which requires uniform types).
// The buffer is grown on demand and meant to be reused across many calls.
inline double dotc_mixed(int n, const double* x, const double* y,
                         std::vector<double>& /*buf*/)
{
    return BlasConnector::dotc(n, x, 1, y, 1);
}

inline double dotc_mixed(int n, const double* x, const float* y,
                         std::vector<double>& buf)
{
    if (static_cast<int>(buf.size()) < n) { buf.resize(n); }
    for (int k = 0; k < n; ++k) { buf[k] = static_cast<double>(y[k]); }
    return BlasConnector::dotc(n, x, 1, buf.data(), 1);
}

// rho(ir) = \sum_{iwt} \phi_i(ir,iwt) * \phi_j^*(ir,iwt)
// phi_j is always double (output of phi_mul_dm); phi_i may be fp32. dotc_mixed
// keeps the inner product in fp64 in either case.
template<typename Tin>
void PhiOperator::phi_dot_phi(
    const Tin*const phi_i,         // phi_i(ir,iwt)
    const double*const phi_j,      // phi_j(ir,iwt)
    double*const rho) const        // rho(ir)
{
    std::vector<double> buf;
    for(int i = 0; i < biggrid_->get_mgrids_num(); ++i)
    {
        rho[mgrid_lidx_[i]] += dotc_mixed(
            cols_, phi_j + i * cols_, phi_i + i * cols_, buf);
    }
}

} // namespace ModuleGint
