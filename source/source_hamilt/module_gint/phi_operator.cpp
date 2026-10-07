#include "phi_operator.h"
#include "source_base/global_function.h"
#include "source_base/matrix.h"
#include "source_base/module_external/blas_connector.h"

#include <cassert>

namespace ModuleGint
{

void PhiOperator::set_bgrid(std::shared_ptr<const BigGrid> biggrid)
{
    biggrid_ = biggrid;
    rows_ = biggrid_->get_mgrids_num();
    cols_ = biggrid_->get_phi_len();

    biggrid_->set_atoms_startidx(atoms_startidx_);
    biggrid_->set_atoms_phi_len(atoms_phi_len_);
    biggrid_->set_mgrids_local_idx(mgrid_lidx_);

    // init is_atom_on_mgrid_ and atom_rcoords_
    const int atoms_num = biggrid_->get_atoms_num();
    atom_rcoords_.resize(atoms_num);
    is_atom_on_mgrid_.resize(biggrid_->get_mgrids_num() * atoms_num);
    for(int i = 0; i < atoms_num; ++i)
    {
        biggrid_->set_atom_relative_coords(biggrid_->get_atom(i), atom_rcoords_[i]);
        for(int j = 0; j < rows_; ++j)
        {
            is_atom_on_mgrid_[i * rows_ + j] = atom_rcoords_[i][j].norm() <= biggrid_->get_atom(i)->get_rcut();
        }
    }

    // init atom_pair_range_
    init_atom_pair_idx_();
}

void PhiOperator::set_phi_dphi(double* phi, double* dphi_x, double* dphi_y, double* dphi_z) const
{
    for(int i = 0; i < biggrid_->get_atoms_num(); ++i)
    {
        const auto atom = biggrid_->get_atom(i);
        atom->set_phi_dphi(atom_rcoords_[i], cols_, phi, dphi_x, dphi_y, dphi_z);
        if(phi != nullptr)
        {
            phi += atom->get_nw();
        }
        dphi_x += atom->get_nw();
        dphi_y += atom->get_nw();
        dphi_z += atom->get_nw();
    }
}

void PhiOperator::set_ddphi(
    double* ddphi_xx, double* ddphi_xy, double* ddphi_xz,
    double* ddphi_yy, double* ddphi_yz, double* ddphi_zz) const
{
    for(int i = 0; i < biggrid_->get_atoms_num(); ++i)
    {
        const auto atom = biggrid_->get_atom(i);
        atom->set_ddphi(atom_rcoords_[i], cols_, ddphi_xx, ddphi_xy, ddphi_xz, ddphi_yy, ddphi_yz, ddphi_zz);
        ddphi_xx += atom->get_nw();
        ddphi_xy += atom->get_nw();
        ddphi_xz += atom->get_nw();
        ddphi_yy += atom->get_nw();
        ddphi_yz += atom->get_nw();
        ddphi_zz += atom->get_nw();
    }
}

void PhiOperator::phi_dot_dphi(
    const double* phi,
    const double* dphi_x,
    const double* dphi_y,
    const double* dphi_z,
    ModuleBase::matrix *fvl) const
{
    for(int i = 0; i < biggrid_->get_atoms_num(); ++i)
    {
        const int iat = biggrid_->get_atom(i)->get_iat();
        const int start_idx = atoms_startidx_[i];
        const int phi_len = atoms_phi_len_[i];
        double rx = 0, ry = 0, rz = 0;
        for(int j = 0; j < biggrid_->get_mgrids_num(); ++j)
        {
            for(int k = 0; k < phi_len; ++k)
            {
                int idx = j * cols_ + start_idx + k;
                const double phi_val = phi[idx];
                rx += phi_val * dphi_x[idx];
                ry += phi_val * dphi_y[idx];
                rz += phi_val * dphi_z[idx];
            }
        }
        fvl[0](iat, 0) += rx * 2;
        fvl[0](iat, 1) += ry * 2;
        fvl[0](iat, 2) += rz * 2;
    }
}

void PhiOperator::phi_dot_dphi_r(
    const double* phi,
    const double* dphi_x,
    const double* dphi_y,
    const double* dphi_z,
    ModuleBase::matrix *svl) const
{
    double sxx = 0, sxy = 0, sxz = 0, syy = 0, syz = 0, szz = 0;
    for(int i = 0; i < biggrid_->get_mgrids_num(); ++i)
    {
        for(int j = 0; j < biggrid_->get_atoms_num(); ++j)
        {
            const int start_idx = atoms_startidx_[j];
            const Vec3d& r3 = atom_rcoords_[j][i];
            for(int k = 0; k < atoms_phi_len_[j]; ++k)
            {
                const int idx = i * cols_ + start_idx + k;
                const double phi_val = phi[idx];
                sxx += phi_val * dphi_x[idx] * r3[0];
                sxy += phi_val * dphi_x[idx] * r3[1];
                sxz += phi_val * dphi_x[idx] * r3[2];
                syy += phi_val * dphi_y[idx] * r3[1];
                syz += phi_val * dphi_y[idx] * r3[2];
                szz += phi_val * dphi_z[idx] * r3[2];
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

void PhiOperator::cal_env_gamma(
    const double* phi,
    const double* wfc,
    const vector<int>& trace_lo,
    double* rho) const
{
    for(int i = 0; i < biggrid_->get_atoms_num(); ++i)
    {
        const auto atom = biggrid_->get_atom(i);
        const int iw_start = atom->get_start_iw();
        const int start_idx = atoms_startidx_[i];
        for(int j = 0; j < biggrid_->get_mgrids_num(); ++j)
        {
            if(is_atom_on_mgrid(i, j))
            {   
                double tmp = 0.0;
                int iw_lo = trace_lo[iw_start];
                for(int iw = 0; iw < atom->get_nw(); ++iw, ++iw_lo)
                {
                    tmp += phi[j * cols_ + start_idx + iw] * wfc[iw_lo];
                }
                rho[mgrid_lidx_[j]] += tmp;
            }
        }
    }
}

void PhiOperator::cal_env_k(
    const double* phi,
    const std::complex<double>* wfc,
    const vector<int>& trace_lo,
    const int ik,
    const int npol,
    const std::vector<Vec3d>& kvec_d,
    const int grid_size,
    std::complex<double>* wfc_r) const
{
    for(int i = 0; i < biggrid_->get_atoms_num(); ++i)
    {
        const auto atom = biggrid_->get_atom(i);
        const int iw_start = atom->get_start_iw();
        // GintAtom::get_R() stores R_Gint = I_bgrid - I_atom, which is the negative
        // of the actual AO-image translation R_AO. The Bloch sum requires
        // exp(+i 2pi k_d dot R_AO), so this is exp(-i 2pi k_d dot R_Gint).
        const Vec3d r_gint(atom->get_R());
        const double arg = -(kvec_d[ik] * r_gint) * ModuleBase::TWO_PI;
        const std::complex<double> kphase = std::complex<double>(cos(arg), sin(arg));
        const int start_idx = atoms_startidx_[i];
        for(int j = 0; j < biggrid_->get_mgrids_num(); ++j)
        {
            if(is_atom_on_mgrid(i, j))
            {
                const int phi_start_idx = j * cols_ + start_idx;
                for (int ipol = 0; ipol < npol; ++ipol)
                {
                    std::complex<double> tmp{0.0, 0.0};
                    for (int iw = 0; iw < atom->get_nw(); ++iw)
                    {
                        const int iw_lo = trace_lo[iw_start + iw * npol + ipol];
                        assert(iw_lo >= 0);
                        tmp += phi[phi_start_idx + iw] * wfc[iw_lo];
                    }
                    wfc_r[ipol * grid_size + mgrid_lidx_[j]] += tmp * kphase;
                }
            }
        }
    }
}


//===============================
// private methods
//===============================
void PhiOperator::init_atom_pair_idx_()
{
    int atoms_num = biggrid_->get_atoms_num();
    atom_pair_range_.resize(atoms_num * (atoms_num + 1) / 2);
    int mgrids_num = biggrid_->get_mgrids_num();
    int atom_pair_idx = 0;
    for(int i = 0; i < atoms_num; ++i)
    {
        // only calculate the upper triangle matrix
        for(int j = i; j < atoms_num; ++j)
        {
            int start_idx = mgrids_num;
            int end_idx = -1;
            for(int mgrid_idx = 0; mgrid_idx < mgrids_num; ++mgrid_idx)
            {
                if(is_atom_on_mgrid(i, mgrid_idx) && is_atom_on_mgrid(j, mgrid_idx))
                {
                    start_idx = mgrid_idx;
                    break;
                }
            }
            for(int mgrid_idx = mgrids_num - 1; mgrid_idx >= 0; --mgrid_idx)
            {
                if(is_atom_on_mgrid(i, mgrid_idx) && is_atom_on_mgrid(j, mgrid_idx))
                {
                    end_idx = mgrid_idx;
                    break;
                }
            }
            atom_pair_range_[atom_pair_idx].first = start_idx;
            atom_pair_range_[atom_pair_idx].second = end_idx;
            atom_pair_idx++;
        }
    }
}

} // namespace ModuleGint

//============================================================
// Template member function implementations (moved from .hpp)
//============================================================

namespace ModuleGint
{

namespace {

// Helper: dispatch a Tin-typed BLAS-GEMM target buffer.
// For Tin=double, write directly into phi_dm (no scratch).
// For Tin=float, allocate fp32 scratch and cast at the end.
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

} // namespace

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

template<typename Tin>
void PhiOperator::phi_mul_dm(
    const Tin*const phi,
    const HContainer<Tin>& dm,
    const bool is_symm,
    double*const phi_dm) const
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
            const auto dm_mat = dm.find_matrix(atom_i->get_iat(), atom_j->get_iat(), r_i-r_j);

            if(dm_mat == nullptr)
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
        T vldr3_mgrid = vl[mgrid_lidx_[i]] * dr3;
        for(int j = 0; j < cols_; j++)
        {
            result[idx] = phi[idx] * vldr3_mgrid;
            idx++;
        }
    }
}

template<typename Tin>
void PhiOperator::phi_mul_phi(
    const Tin*const phi_i,
    const Tin*const phi_j,
    HContainer<double>& hr,
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

            if(part==TriPart::Upper && iat_i>iat_j)
            {
                continue;
            }
            else if(part==TriPart::Lower && iat_i<iat_j)
            {
                continue;
            }

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
// when y is fp32 it is upcast into the caller-provided fp64 scratch buffer.
namespace {

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

} // namespace

template<typename Tin>
void PhiOperator::phi_dot_phi(
    const Tin*const phi_i,
    const double*const phi_j,
    double*const rho) const
{
    std::vector<double> buf;
    for(int i = 0; i < biggrid_->get_mgrids_num(); ++i)
    {
        rho[mgrid_lidx_[i]] += dotc_mixed(
            cols_, phi_j + i * cols_, phi_i + i * cols_, buf);
    }
}

//============================================================
// Explicit template instantiation
//============================================================

template void PhiOperator::set_phi<double>(double* phi) const;
template void PhiOperator::set_phi<float>(float* phi) const;

template void PhiOperator::phi_mul_dm<double>(
    const double*const phi, const HContainer<double>& dm,
    const bool is_symm, double*const phi_dm) const;
template void PhiOperator::phi_mul_dm<float>(
    const float*const phi, const HContainer<float>& dm,
    const bool is_symm, double*const phi_dm) const;

template void PhiOperator::phi_mul_vldr3<double>(
    const double*const vl, const double dr3,
    const double*const phi, double*const result) const;
template void PhiOperator::phi_mul_vldr3<float>(
    const float*const vl, const float dr3,
    const float*const phi, float*const result) const;

template void PhiOperator::phi_mul_phi<double>(
    const double*const phi_i, const double*const phi_j,
    HContainer<double>& hr, const TriPart part) const;
template void PhiOperator::phi_mul_phi<float>(
    const float*const phi_i, const float*const phi_j,
    HContainer<double>& hr, const TriPart part) const;

template void PhiOperator::phi_dot_phi<double>(
    const double*const phi_i, const double*const phi_j,
    double*const rho) const;
template void PhiOperator::phi_dot_phi<float>(
    const float*const phi_i, const double*const phi_j,
    double*const rho) const;

} // namespace ModuleGint
