#include "phi_operator.h"
#include "source_base/global_function.h"
#include "source_base/matrix.h"

namespace ModuleGint
{

void PhiOperator::set_bgrid(std::shared_ptr<const BigGrid> biggrid)
{
    biggrid_ = biggrid;
    rows_ = biggrid_->get_mgrids_num();
    cols_ = biggrid_->get_phi_len();

    biggrid_->set_atoms_startidx(atoms_startidx_);
    biggrid_->set_atoms_phi_len(atoms_phi_len_);
    biggrid_->set_mgrids_local_idx(meshgrids_local_idx_);

    // init is_atom_on_mgrid_ and atoms_relative_coords_
    const int atoms_num = biggrid_->get_atoms_num();
    atoms_relative_coords_.resize(atoms_num);
    is_atom_on_mgrid_.resize(biggrid_->get_mgrids_num() * atoms_num);
    for(int i = 0; i < atoms_num; ++i)
    {
        biggrid_->set_atom_relative_coords(biggrid_->get_atom(i), atoms_relative_coords_[i]);
        for(int j = 0; j < rows_; ++j)
        {
            is_atom_on_mgrid_[i * rows_ + j] = atoms_relative_coords_[i][j].norm() <= biggrid_->get_atom(i)->get_rcut();
        }
    }

    // init atom_pair_start_end_idx_
    init_atom_pair_start_end_idx_();
}

void PhiOperator::set_phi_dphi(double* phi, double* dphi_x, double* dphi_y, double* dphi_z) const
{
    for(int i = 0; i < biggrid_->get_atoms_num(); ++i)
    {
        const auto atom = biggrid_->get_atom(i);
        atom->set_phi_dphi(atoms_relative_coords_[i], cols_, phi, dphi_x, dphi_y, dphi_z);
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
        atom->set_ddphi(atoms_relative_coords_[i], cols_, ddphi_xx, ddphi_xy, ddphi_xz, ddphi_yy, ddphi_yz, ddphi_zz);
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
            const Vec3d& r3 = atoms_relative_coords_[j][i];
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
                rho[meshgrids_local_idx_[j]] += tmp;
            }
        }
    }
}

void PhiOperator::cal_env_k(
    const double* phi,
    const std::complex<double>* wfc,
    const vector<int>& trace_lo,
    const int ik,
    const int nspin,
    const int npol,
    const int lgd,
    const std::vector<Vec3d>& kvec_c,
    const std::vector<Vec3d>& kvec_d,
    double* rho) const
{
    for(int i = 0; i < biggrid_->get_atoms_num(); ++i)
    {
        const auto atom = biggrid_->get_atom(i);
        const int iw_start = atom->get_start_iw();
        const Vec3d R(atom->get_unitcell_idx());
        const double arg = (kvec_d[ik] * R) * ModuleBase::TWO_PI;
        const std::complex<double> kphase = std::complex<double>(cos(arg), sin(arg));
        const int start_idx = atoms_startidx_[i];
        for(int j = 0; j < biggrid_->get_mgrids_num(); ++j)
        {
            if(is_atom_on_mgrid(i, j))
            {   
                std::complex<double> tmp{0.0, 0.0};
                int phi_start_idx = j * cols_ + start_idx;

                int iw_lo = 0;
                if (nspin == 4) // is it a simple add of 2 spins?
                {
                    for (int is = 0; is < 2; ++is)
                    {
                        iw_lo = trace_lo[iw_start] / npol + lgd / npol * is;
                        for (int iw = 0; iw < atom->get_nw(); ++iw, ++iw_lo)
                        {
                            tmp += std::complex<double>(phi[phi_start_idx + iw], 0.0) * wfc[iw_lo] * kphase;
                        }
                    }
                }
                else
                {
                    iw_lo = trace_lo[iw_start];
                    for (int iw = 0; iw < atom->get_nw(); ++iw, ++iw_lo)
                    {
                        tmp += std::complex<double>(phi[phi_start_idx + iw], 0.0) * wfc[iw_lo] * kphase;
                    }
                }
                rho[meshgrids_local_idx_[j]] += tmp.real();
            }
        }
    }
}


//===============================
// private methods
//===============================
void PhiOperator::init_atom_pair_start_end_idx_()
{
    int atoms_num = biggrid_->get_atoms_num();
    atom_pair_start_end_idx_.resize(atoms_num * (atoms_num + 1) / 2);
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
            atom_pair_start_end_idx_[atom_pair_idx].first = start_idx;
            atom_pair_start_end_idx_[atom_pair_idx].second = end_idx;
            atom_pair_idx++;
        }
    }
}

} // namespace ModuleGint