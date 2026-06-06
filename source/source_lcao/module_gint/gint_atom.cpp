#include "source_base/ylm.h"
#include "gint_atom.h"
#include "source_cell/unitcell.h"
#include "gint_helper.h"

namespace ModuleGint
{
GintAtom::GintAtom(
    const Atom* atom,
    int it, int ia, int iat,
    Vec3i biggrid_idx,
    Vec3i unitcell_idx,
    Vec3d tau_in_biggrid,
    const Numerical_Orbital* orb,
    const UnitCell* ucell)
: atom_(atom), it_(it), ia_(ia), iat_(iat), biggrid_idx_(biggrid_idx),
  unitcell_idx_(unitcell_idx), tau_in_biggrid_(tau_in_biggrid),
  orb_(orb), ucell_(ucell)
{
    p_psi_uniform_.resize(atom_->nw);
    p_dpsi_uniform_.resize(atom_->nw);
    p_ddpsi_uniform_.resize(atom_->nw);
    radial_blocks_.reserve(atom_->nw);
    for (int iw=0; iw < atom_->nw; ++iw)
    {
        if ( atom_->iw2_new[iw] )
        {
            int l = atom_->iw2l[iw];
            int n = atom_->iw2n[iw];
            const auto& phi_ln = orb_->PhiLN(l, n);
            p_psi_uniform_[iw] = phi_ln.psi_uniform.data();
            p_dpsi_uniform_[iw] = phi_ln.dpsi_uniform.data();
            p_ddpsi_uniform_[iw] = phi_ln.ddpsi_uniform.data();

            RadialBlock block;
            block.begin_iw = iw;
            block.size = 2 * l + 1;
            block.ylm_begin = atom_->iw2_ylm[iw];
            block.psi_uniform = p_psi_uniform_[iw];
            block.dpsi_uniform = p_dpsi_uniform_[iw];
            radial_blocks_.push_back(block);
        }
    }
}

template <typename T>
void GintAtom::set_phi(const std::vector<Vec3d>& coords, const int stride, T* phi) const
{
    const int num_mgrids = coords.size();
    const double dr_uniform = orb_->PhiLN(0, 0).dr_uniform;

    std::vector<double> ylma;
    const auto* blocks = radial_blocks_.data();
    const int num_blocks = radial_blocks_.size();

    for(int im = 0; im < num_mgrids; im++)
    {
        const Vec3d& coord = coords[im];
        const double dist = coord.norm() < 1e-9 ? 1e-9 : coord.norm();
        if(dist > orb_->getRcut())
        {   
            ModuleBase::GlobalFunc::ZEROS(phi + im * stride, atom_->nw);
        }
        else
        {
            ModuleBase::Ylm::sph_harm(atom_->nwl, coord.x/dist, coord.y/dist, coord.z/dist, ylma);

            const double position = dist / dr_uniform;
            const int ip = static_cast<int>(position);
            const double dx = position - ip;
            const double dx2 = dx * dx;
            const double dx3 = dx2 * dx;

            const double c3 = 3.0 * dx2 - 2.0 * dx3;
            const double c1 = 1.0 - c3;
            const double c2 = (dx - 2.0 * dx2 + dx3) * dr_uniform;
            const double c4 = (dx3 - dx2) * dr_uniform;

            T* phi_row = phi + im * stride;
            for (int ib = 0; ib < num_blocks; ++ib)
            {
                const auto& block = blocks[ib];
                const double* psi_uniform = block.psi_uniform;
                const double* dpsi_uniform = block.dpsi_uniform;
                const double psi = c1 * psi_uniform[ip] + c2 * dpsi_uniform[ip]
                    + c3 * psi_uniform[ip + 1] + c4 * dpsi_uniform[ip + 1];

                const int begin_iw = block.begin_iw;
                const int end_iw = begin_iw + block.size;
                int idx_lm = block.ylm_begin;
                for (int iw = begin_iw; iw < end_iw; ++iw, ++idx_lm)
                {
                    phi_row[iw] = psi * ylma[idx_lm];
                }
            }
        }
    }
}

template <typename T>
void GintAtom::set_phi_dphi(
    const std::vector<Vec3d>& coords, const int stride,
    T* phi, T* dphi_x, T* dphi_y, T* dphi_z) const
{
    if (phi != nullptr) {
        phi = (T*)__builtin_assume_aligned(phi, 64);
    }
    if (dphi_x != nullptr) {
        dphi_x = (T*)__builtin_assume_aligned(dphi_x, 64);
        dphi_y = (T*)__builtin_assume_aligned(dphi_y, 64);
        dphi_z = (T*)__builtin_assume_aligned(dphi_z, 64);
    }
    const int num_mgrids = coords.size();
    const double dr_uniform = orb_->PhiLN(0, 0).dr_uniform;
    
    const int nylm = std::pow(atom_->nwl + 1, 2);
    std::vector<double> rly(nylm);
    
    // 展平为一维连续内存
    std::vector<double> grly_data(nylm * 3);
    
    // 构造代理二维指针数组以适配底层接口
    std::vector<double*> grly_ptrs(nylm);
    for(int i = 0; i < nylm; ++i) {
        grly_ptrs[i] = &grly_data[i * 3];
    }
    
    for(int im = 0; im < num_mgrids; im++)
    {
        const Vec3d& coord = coords[im];
        const double dist = coord.norm() < 1e-9 ? 1e-9 : coord.norm();

        if(dist > orb_->getRcut())
        {
            if(phi != nullptr)
            {
                ModuleBase::GlobalFunc::ZEROS(phi + im * stride, atom_->nw);
            }
            ModuleBase::GlobalFunc::ZEROS(dphi_x + im * stride, atom_->nw);
            ModuleBase::GlobalFunc::ZEROS(dphi_y + im * stride, atom_->nw);
            ModuleBase::GlobalFunc::ZEROS(dphi_z + im * stride, atom_->nw);
        }
        else
        {
            // 使用代理指针数组传入，底层函数会将结果写入 grly_data 中
            ModuleBase::Ylm::grad_rl_sph_harm(atom_->nwl, coord.x, coord.y, coord.z, rly.data(), grly_ptrs.data());

            const double position = dist / dr_uniform;
            const int ip = static_cast<int>(position);
            const double x0 = position - ip;
            const double x1 = 1.0 - x0;
            const double x2 = 2.0 - x0;
            const double x3 = 3.0 - x0;
            const double x12 = x1 * x2 / 6;
            const double x03 = x0 * x3 / 2;

            double tmp, dtmp;
            
            // 对每个轨道 iw 进行计算
            for(int iw = 0; iw < atom_->nw; ++iw)
            {
                if(atom_->iw2_new[iw])
                {
                    auto psi_uniform = p_psi_uniform_[iw];
                    auto dpsi_uniform = p_dpsi_uniform_[iw];

                    tmp = x12 * (psi_uniform[ip] * x3 + psi_uniform[ip + 3] * x0)
                        + x03 * (psi_uniform[ip + 1] * x2 - psi_uniform[ip + 2] * x1);

                    dtmp = x12 * (dpsi_uniform[ip] * x3 + dpsi_uniform[ip + 3] * x0)
                        + x03 * (dpsi_uniform[ip + 1] * x2 - dpsi_uniform[ip + 2] * x1);
                } 

                const int ll = atom_->iw2l[iw];
                const int idx_lm = atom_->iw2_ylm[iw];

                double rl = 1.0;
                switch (ll) {
                    case 4: rl = dist * dist * dist * dist; break;
                    case 3: rl = dist * dist * dist; break;
                    case 2: rl = dist * dist; break;
                    case 1: rl = dist; break;
                    case 0: rl = 1.0; break;
                    default: rl = pow_int(dist, ll); 
                }
                
                const double tmprl = tmp / rl;
                const double tmpdphi_rly = (dtmp - tmp * ll / dist) / rl / dist; 

                // 移除错误的内部 im 循环，直接对当前网格点(im)和当前轨道(iw)赋值
                if(phi != nullptr)
                {
                    phi[im * stride + iw] = tmprl * rly[idx_lm];
                }
                
                if(dphi_x != nullptr)
                {
                    double tmpdphi_rly_val = tmpdphi_rly * rly[idx_lm];
                    
                    // 使用一维数组偏移寻址
                    dphi_x[im * stride + iw] = tmpdphi_rly_val * coord.x + tmprl * grly_data[idx_lm * 3 + 0];
                    dphi_y[im * stride + iw] = tmpdphi_rly_val * coord.y + tmprl * grly_data[idx_lm * 3 + 1];
                    dphi_z[im * stride + iw] = tmpdphi_rly_val * coord.z + tmprl * grly_data[idx_lm * 3 + 2];
                }
            }
        }
    }
}

// explicit instantiation
template void GintAtom::set_phi(const std::vector<Vec3d>& coords, const int stride, float* phi) const;
template void GintAtom::set_phi(const std::vector<Vec3d>& coords, const int stride, double* phi) const;
template void GintAtom::set_phi(const std::vector<Vec3d>& coords, const int stride, std::complex<double>* phi) const;
template void GintAtom::set_phi_dphi(const std::vector<Vec3d>& coords, const int stride, float* phi, float* dphi_x, float* dphi_y, float* dphi_z) const;
template void GintAtom::set_phi_dphi(const std::vector<Vec3d>& coords, const int stride, double* phi, double* dphi_x, double* dphi_y, double* dphi_z) const;

} // namespace ModuleGint