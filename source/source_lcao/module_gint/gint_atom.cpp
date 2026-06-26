#include "source_base/ylm.h"
#include "gint_atom.h"
#include "source_cell/unitcell.h"
#include "gint_helper.h"

#ifdef GINT_SET_PHI_OPENMP
#ifdef _OPENMP
#include <omp.h>
#endif
#endif

namespace ModuleGint
{

namespace
{

template <typename T>
void set_phi_at_meshgrid(
    const Atom* atom,
    const int nw,
    const int nwl,
    const double dr_uniform,
    const double rcut,
    const GintAtom::RadialBlock* blocks,
    const int num_blocks,
    const Vec3d& coord,
    T* phi_row,
    std::vector<double>& ylma)
{
    double dist = coord.norm();
    if (dist < 1e-9)
    {
        dist = 1e-9;
    }

    if (dist > rcut)
    {
        ModuleBase::GlobalFunc::ZEROS(phi_row, nw);
        return;
    }

    const double inv_dist = 1.0 / dist;
    ModuleBase::Ylm::sph_harm(
        nwl,
        coord.x * inv_dist,
        coord.y * inv_dist,
        coord.z * inv_dist,
        ylma);

    const double position = dist / dr_uniform;
    const int ip = static_cast<int>(position);
    const double dx = position - ip;
    const double dx2 = dx * dx;
    const double dx3 = dx2 * dx;

    const double c3 = 3.0 * dx2 - 2.0 * dx3;
    const double c1 = 1.0 - c3;
    const double c2 = (dx - 2.0 * dx2 + dx3) * dr_uniform;
    const double c4 = (dx3 - dx2) * dr_uniform;

    for (int ib = 0; ib < num_blocks; ++ib)
    {
        const auto& block = blocks[ib];
        const double* psi_uniform = block.psi_uniform;
        const double* dpsi_uniform = block.dpsi_uniform;
        const double psi = c1 * psi_uniform[ip] + c2 * dpsi_uniform[ip]
            + c3 * psi_uniform[ip + 1] + c4 * dpsi_uniform[ip + 1];

        const int begin_iw = block.begin_iw;
        const int end_iw = begin_iw + block.size;
        const int ylm_begin = block.ylm_begin;
#if defined(GINT_SET_PHI_SIMD)
#pragma omp simd
#endif
        for (int iw = begin_iw; iw < end_iw; ++iw)
        {
            phi_row[iw] = psi * ylma[iw - begin_iw + ylm_begin];
        }
    }
}

} // namespace
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
            // The first orbital in each radial block always starts from m = 0.
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
    if (num_mgrids == 0)
    {
        return;
    }

    const double dr_uniform = orb_->PhiLN(0, 0).dr_uniform;
    const double rcut = orb_->getRcut();
    const int nw = atom_->nw;
    const int nwl = atom_->nwl;
    const auto* blocks = radial_blocks_.data();
    const int num_blocks = radial_blocks_.size();

#if defined(GINT_SET_PHI_OPENMP) && defined(_OPENMP)
    // cal_gint_rho/vl 已在 BigGrid 层做了 OpenMP 并行；默认关闭嵌套并行，
    // 仅在非并行上下文中才启动格点级并行（例如单元测试或独立 profiling）。
    if (omp_in_parallel())
    {
        std::vector<double> ylma;
        for (int im = 0; im < num_mgrids; ++im)
        {
            set_phi_at_meshgrid(
                atom_,
                nw,
                nwl,
                dr_uniform,
                rcut,
                blocks,
                num_blocks,
                coords[im],
                phi + im * stride,
                ylma);
        }
    }
    else
    {
#pragma omp parallel if(num_mgrids > 1)
        {
            std::vector<double> ylma;
#pragma omp for schedule(static)
            for (int im = 0; im < num_mgrids; ++im)
            {
                set_phi_at_meshgrid(
                    atom_,
                    nw,
                    nwl,
                    dr_uniform,
                    rcut,
                    blocks,
                    num_blocks,
                    coords[im],
                    phi + im * stride,
                    ylma);
            }
        }
    }
#else
    std::vector<double> ylma;
    for (int im = 0; im < num_mgrids; ++im)
    {
        set_phi_at_meshgrid(
            atom_,
            nw,
            nwl,
            dr_uniform,
            rcut,
            blocks,
            num_blocks,
            coords[im],
            phi + im * stride,
            ylma);
    }
#endif
}

template <typename T>
void GintAtom::set_phi_dphi(
    const std::vector<Vec3d>& coords, const int stride,
    T* phi, T* dphi_x, T* dphi_y, T* dphi_z) const
{
    const int num_mgrids = coords.size();
    
    // orb_ does not have the member variable dr_uniform
    const double dr_uniform = orb_->PhiLN(0, 0).dr_uniform;
    
    const int nylm = std::pow(atom_->nwl + 1, 2);
    std::vector<double> rly(nylm);
    std::vector<double> grly(nylm * 3);
    
    for(int im = 0; im < num_mgrids; im++)
    {
        const Vec3d& coord = coords[im];
        // 1e-9 is to avoid division by zero
        const double dist = coord.norm() < 1e-9 ? 1e-9 : coord.norm();

        if(dist > orb_->getRcut())
        {
            // if the distance is larger than the cutoff radius,
            // the wave function values are all zeros
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
            // spherical harmonics
            // TODO: vectorize the sph_harm function, 
            // the vectorized function can be called once for all meshgrids in a biggrid
            ModuleBase::Ylm::grad_rl_sph_harm(atom_->nwl, coord.x, coord.y, coord.z, rly.data(), grly.data());

            // interpolation
            const double position = dist / dr_uniform;
            const int ip = static_cast<int>(position);
            const double x0 = position - ip;
            const double x1 = 1.0 - x0;
            const double x2 = 2.0 - x0;
            const double x3 = 3.0 - x0;
            const double x12 = x1 * x2 / 6;
            const double x03 = x0 * x3 / 2;

            double tmp, dtmp;
            for(int iw = 0; iw < atom_->nw; ++iw)
            {
                // this is a new 'l', we need 1D orbital wave
                // function from interpolation method.
                if(atom_->iw2_new[iw])
                {
                    auto psi_uniform = p_psi_uniform_[iw];
                    auto dpsi_uniform = p_dpsi_uniform_[iw];
                    // use Polynomia Interpolation method to get the
                    // wave functions

                    tmp = x12 * (psi_uniform[ip] * x3 + psi_uniform[ip + 3] * x0)
                        + x03 * (psi_uniform[ip + 1] * x2 - psi_uniform[ip + 2] * x1);

                    dtmp = x12 * (dpsi_uniform[ip] * x3 + dpsi_uniform[ip + 3] * x0)
                        + x03 * (dpsi_uniform[ip + 1] * x2 - dpsi_uniform[ip + 2] * x1);
                } // new l is used.

                // get the 'l' of this localized wave function
                const int ll = atom_->iw2l[iw];
                const int idx_lm = atom_->iw2_ylm[iw];

                const double rl = pow_int(dist, ll);
                const double tmprl = tmp / rl;

                // 3D wave functions
                if(phi != nullptr)
                {
                    phi[im * stride + iw] = tmprl * rly[idx_lm];
                }
                
                // derivative of wave functions with respect to atom positions.
                const double tmpdphi_rly = (dtmp - tmp * ll / dist) / rl * rly[idx_lm] / dist;

                dphi_x[im * stride + iw] =  tmpdphi_rly * coord.x + tmprl * grly[idx_lm*3];
                dphi_y[im * stride + iw] =  tmpdphi_rly * coord.y + tmprl * grly[idx_lm*3 + 1];
                dphi_z[im * stride + iw] =  tmpdphi_rly * coord.z + tmprl * grly[idx_lm*3 + 2];
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
}
