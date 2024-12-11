#include "gint_atom.h"
#include "module_base/timer.h"
#include "module_base/ylm.h"
#include "module_base/array_pool.h"

namespace Gint
{

template <typename T>
void GintAtom::get_psi(const std::vector<Vec3d> coords, const int stride, T* psi)
{
    ModuleBase::timer::tick("GintAtom", "get_psi");

    const int num_mgrids = coords.size();

    // orb_ does not have the member variable dr_uniform
    const double dr_uniform = orb_->PhiLN(0, 0).dr_uniform;

    // store the pointer to reduce repeated address fetching
    std::vector<const double*> p_psi_uniform(atom_->nw);
    std::vector<const double*> p_dpsi_uniform(atom_->nw);
    for(int iw = 0; iw < atom_->nw; iw++)
    {
        if(atom_->iw2_new[iw])
        {
            int l = atom_->iw2l[iw];
            int n = atom_->iw2n[iw];
            p_psi_uniform[iw] = orb_->PhiLN(l, n).psi_uniform.data();
            p_dpsi_uniform[iw] = orb_->PhiLN(l, n).dpsi_uniform.data();
        }
    }

    // store the spherical harmonics
    // it's outside the loop to reduce the vector allocation overhead
    std::vector<double> ylma;

    for(int im = 0; im < num_mgrids; im++)
    {
        const Vec3d& coord = coords[im];
        const double dist = coord.norm();

        // avoid division by zero
        if(dist < 1e-9)
        {
            dist = 1e-9;
        }

        if(dist > orb_->getRcut())
        {   
            // if the distance is larger than the cutoff radius,
            // the wave function values are all zeros
            ModuleBase::GlobalFunc::ZEROS(psi + im * stride, atom_->nw);
        }
        else
        {
            // spherical harmonics
            // TODO: vectorize the sph_harm function, 
            // the vectorized function can be called once for all meshgrids in a biggrid
            ModuleBase::Ylm::sph_harm(atom_->nwl, coord.x/dist, coord.y/dist, coord.z/dist, ylma);

            // interpolation

            // these parameters are related to interpolation
            // because once the distance from atom to grid point is known,
            // we can obtain the parameters for interpolation and
            // store them first! these operations can save lots of efforts.
            const double position = dist / dr_uniform;
            const int ip = static_cast<int>(position);
            const double dx = position - ip;
            const double dx2 = dx * dx;
            const double dx3 = dx2 * dx;

            const double c3 = 3.0 * dx2 - 2.0 * dx3;
            const double c1 = 1.0 - c3;
            const double c2 = (dx - 2.0 * dx2 + dx3) * dr_uniform;
            const double c4 = (dx3 - dx2) * dr_uniform;

            double phi = 0;
            
            for(int iw = 0; iw < atom_->nw; iw++)
            {
                if(atom_->iw2_new[iw])
                {
                    auto psi_uniform = p_psi_uniform[iw];
                    auto dpsi_uniform = p_dpsi_uniform[iw];
                    phi = c1 * psi_uniform[ip] + c2 * dpsi_uniform[ip]
                        + c3 * psi_uniform[ip + 1] + c4 * dpsi_uniform[ip + 1];
                }
                psi[im * stride + iw] = phi * ylma[atom_->iw2_ylm[iw]];
            }
        }
    }

    ModuleBase::timer::tick("GintAtom", "get_psi");
}

template <typename T>
void GintAtom::get_psi_dpsir(
    const std::vector<Vec3d> coords, const int stride,
    T* psi, T* dpsi_x, T* dpsi_y, T* dpsi_z)
{
    ModuleBase::timer::tick("GintAtom", "get_psi_dpsir");

    const int num_mgrids = coords.size();
    
    // orb_ does not have the member variable dr_uniform
    const double dr_uniform = orb_->PhiLN(0, 0).dr_uniform;

    // store the pointer to reduce repeated address fetching
    std::vector<const double*> p_psi_uniform(atom_->nw);
    std::vector<const double*> p_dpsi_uniform(atom_->nw);
    std::vector<int> psi_nr_uniform(atom_->nw);
    for (int iw=0; iw< atom->nw; ++iw)
    {
        if ( atom->iw2_new[iw] )
        {
            int l = atom_->iw2l[iw];
            int n = atom_->iw2n[iw];
            p_psi_uniform[iw] = orb_->PhiLN(l, n).psi_uniform.data();
            p_dpsi_uniform[iw] = orb_->PhiLN(l, n).dpsi_uniform.data();
            psi_nr_uniform[iw] = orb_->PhiLN(l, n).nr_uniform;
        }
    }
    
    std::vector<double> rly(std::pow(atom_->nwl + 1, 2));
    ModuleBase::Array_Pool<double> grly(std::pow(atom_->nwl + 1, 2), 3);
    
    for(int im = 0; im < num_mgrids; im++)
    {
        const Vec3d& coord = coords[im];
        const double dist = coord.norm();

        // avoid division by zero
        if(dist < 1e-9)
        {
            dist = 1e-9;
        }

        if(dist > orb_->getRcut())
        {
            // if the distance is larger than the cutoff radius,
            // the wave function values are all zeros
            ModuleBase::GlobalFunc::ZEROS(dpsi_x + im * stride, atom_->nw);
            ModuleBase::GlobalFunc::ZEROS(dpsi_y + im * stride, atom_->nw);
            ModuleBase::GlobalFunc::ZEROS(dpsi_z + im * stride, atom_->nw);
        }
        else
        {
            // spherical harmonics
            // TODO: vectorize the sph_harm function, 
            // the vectorized function can be called once for all meshgrids in a biggrid
            ModuleBase::Ylm::grad_rl_sph_harm(atom_->nwl, coord.x, coord.y, coord.z, rly, grly.get_ptr_2D());

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
                    auto psi_uniform = p_psi_uniform[iw];
                    auto dpsi_uniform = p_dpsi_uniform[iw];

                    if(ip >= psi_nr_uniform[iw] - 4)
                    {
                        tmp = dtmp = 0.0;
                    }
                    else
                    {
                        // use Polynomia Interpolation method to get the
                        // wave functions

                        tmp = x12 * (psi_uniform[ip] * x3 + psi_uniform[ip + 3] * x0)
                            + x03 * (psi_uniform[ip + 1] * x2 - psi_uniform[ip + 2] * x1);

                        dtmp = x12 * (dpsi_uniform[ip] * x3 + dpsi_uniform[ip + 3] * x0)
                            + x03 * (dpsi_uniform[ip + 1] * x2 - dpsi_uniform[ip + 2] * x1);
                    }
                } // new l is used.

                // get the 'l' of this localized wave function
                const int ll = atom_->iw2l[iw];
                const int idx_lm = atom_->iw2_ylm[iw];

                const double rl = pow_int(dist, ll);
                const double tmprl = tmp / rl;

                // 3D wave functions
                psi[im * stride + iw] = tmprl * rly[idx_lm];

                // derivative of wave functions with respect to atom positions.
                const double tmpdphi_rly = (dtmp - tmp * ll / dist) / rl * rly[idx_lm] / dist;

                dpsi_x[im * stride + iw] =  tmpdphi_rly * coord.x + tmprl * grly[idx_lm][0];
                dpsi_y[im * stride + iw] =  tmpdphi_rly * coord.y + tmprl * grly[idx_lm][1];
                dpsi_z[im * stride + iw] =  tmpdphi_rly * coord.z + tmprl * grly[idx_lm][2];
            }
        }
    }

    ModuleBase::timer::tick("GintAtom", "get_psi_dpsir");
}

// explicit instantiation
template void GintAtom::get_psi(const std::vector<Vec3d> coords, const int stride, double* psi);
template void GintAtom::get_psi_dpsir(const std::vector<Vec3d> coords, const int stride, double* psi, double* dpsi_x, double* dpsi_y, double* dpsi_z);
}