#include "module_base/array_pool.h"
#include "module_base/timer.h"
#include "module_base/ylm.h"
#include "gint_atom.h"
#include "gint_helper.h"

namespace ModuleGint
{

template <typename T>
void GintAtom::set_ddphi(
    const std::vector<Vec3d>& coords, const int stride,
    T* ddphi_xx, T* ddphi_xy, T* ddphi_xz,
    T* ddphi_yy, T* ddphi_yz, T* ddphi_zz) const
{
    ModuleBase::timer::tick("GintAtom", "set_ddphi");
    
    const int num_mgrids = coords.size();

    // orb_ does not have the member variable dr_uniform
    const double dr_uniform = orb_->PhiLN(0, 0).dr_uniform;

    // store the pointer to reduce repeated address fetching
    std::vector<const double*> p_psi_uniform(atom_->nw);
    std::vector<const double*> p_dpsi_uniform(atom_->nw);
    std::vector<const double*> p_ddpsi_uniform(atom_->nw);
    std::vector<int> phi_nr_uniform(atom_->nw);
    for (int iw=0; iw< atom_->nw; ++iw)
    {
        if ( atom_->iw2_new[iw] )
        {
            int l = atom_->iw2l[iw];
            int n = atom_->iw2n[iw];
            p_psi_uniform[iw] = orb_->PhiLN(l, n).psi_uniform.data();
            p_dpsi_uniform[iw] = orb_->PhiLN(l, n).dpsi_uniform.data();
            p_ddpsi_uniform[iw] = orb_->PhiLN(l, n).ddpsi_uniform.data();
            phi_nr_uniform[iw] = orb_->PhiLN(l, n).nr_uniform;
        }
    }

    std::vector<double> rly(std::pow(atom_->nwl + 1, 2));
    ModuleBase::Array_Pool<double> grly(std::pow(atom_->nwl + 1, 2), 3);
    // TODO: A better data structure such as a 3D tensor can be used to store dphi
    std::vector<std::vector<std::vector<double>>> dphi(atom_->nw, std::vector<std::vector<double>>(6, std::vector<double>(3)));
    Vec3d coord1;
    ModuleBase::Array_Pool<double> displ(6, 3);
    displ[0][0] = 0.0001; // in x direction
    displ[1][0] = -0.0001;
    displ[2][1] = 0.0001; // in y direction
    displ[3][1] = -0.0001;
    displ[4][2] = 0.0001; // in z direction
    displ[5][2] = -0.0001;

    for(int im = 0; im < num_mgrids; im++)
    {
        const Vec3d& coord = coords[im];
        // 1e-9 is to avoid division by zero
        const double dist = coord.norm() < 1e-9 ? 1e-9 : coord.norm();

        if(dist > orb_->getRcut())
        {
            // if the distance is larger than the cutoff radius,
            // the wave function values are all zeros
            ModuleBase::GlobalFunc::ZEROS(ddphi_xx + im * stride, atom_->nw);
            ModuleBase::GlobalFunc::ZEROS(ddphi_xy + im * stride, atom_->nw);
            ModuleBase::GlobalFunc::ZEROS(ddphi_xz + im * stride, atom_->nw);
            ModuleBase::GlobalFunc::ZEROS(ddphi_yy + im * stride, atom_->nw);
            ModuleBase::GlobalFunc::ZEROS(ddphi_yz + im * stride, atom_->nw);
            ModuleBase::GlobalFunc::ZEROS(ddphi_zz + im * stride, atom_->nw);
            continue;
        }

        for(int i = 0; i < 6; i++)
        {
            coord1[0] = coord[0] + displ[i][0];
            coord1[1] = coord[1] + displ[i][1];
            coord1[2] = coord[2] + displ[i][2];

            // sphereical harmonics
            ModuleBase::Ylm::grad_rl_sph_harm(atom_->nwl, coord1[0], coord1[1], coord1[2], rly.data(), grly.get_ptr_2D());

            const double dist1 = coord1.norm() < 1e-9 ? 1e-9 : coord1.norm();

            const double position = dist1 / dr_uniform;
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
                if(atom_->iw2_new[iw])
                {
                    auto psi_uniform = p_psi_uniform[iw];
                    auto dpsi_uniform = p_dpsi_uniform[iw];

                    if(ip >= phi_nr_uniform[iw] - 4)
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
                }

                // get the 'l' of this localized wave function
                const int ll = atom_->iw2l[iw];
                const int idx_lm = atom_->iw2_ylm[iw];

                const double rl = pow_int(dist1, ll);

                // derivative of wave functions with respect to atom positions.
                const double tmpdphi_rly = (dtmp - tmp * ll / dist1) / rl * rly[idx_lm] / dist1;
                const double tmprl = tmp / rl;

                dphi[iw][i][0] =  tmpdphi_rly * coord1[0] + tmprl * grly[idx_lm][0];
                dphi[iw][i][1] =  tmpdphi_rly * coord1[1] + tmprl * grly[idx_lm][1];
                dphi[iw][i][2] =  tmpdphi_rly * coord1[2] + tmprl * grly[idx_lm][2];
            } // end iw
        }  // end i

        for(int iw = 0; iw < atom_->nw; iw++)
        {
            int idx = im * stride + iw;
            ddphi_xx[idx] = (dphi[iw][0][0] - dphi[iw][1][0]) / 0.0002;
            ddphi_xy[idx]
                = ((dphi[iw][2][0] - dphi[iw][3][0]) + (dphi[iw][0][1] - dphi[iw][1][1])) / 0.0004;
            ddphi_xz[idx]
                = ((dphi[iw][4][0] - dphi[iw][5][0]) + (dphi[iw][0][2] - dphi[iw][1][2])) / 0.0004;
            ddphi_yy[idx] = (dphi[iw][2][1] - dphi[iw][3][1]) / 0.0002;
            ddphi_yz[idx]
                = ((dphi[iw][4][1] - dphi[iw][5][1]) + (dphi[iw][2][2] - dphi[iw][3][2])) / 0.0004;
            ddphi_zz[idx] = (dphi[iw][4][2] - dphi[iw][5][2]) / 0.0002;
        }

        // else
        //     // the analytical method for evaluating 2nd derivatives
        //     // it is not used currently
        //     {
        //         // Add it here, but do not run it. If there is a need to run this code 
        //         // in the future, include it in the previous initialization process.
        //         for (int iw=0; iw< atom->nw; ++iw)
        //         {
        //             if ( atom->iw2_new[iw] )
        //             {
        //                 it_ddpsi_uniform[iw] = gt.d2phi_u[it*gt.nwmax + iw].data();
        //             }
        //         }
        //         // End of code addition section.

        //         std::vector<std::vector<double>> hrly;
        //         ModuleBase::Ylm::grad_rl_sph_harm(ucell.atoms[it].nwl, dr[0], dr[1], dr[2], rly, grly.data());
        //         ModuleBase::Ylm::hes_rl_sph_harm(ucell.atoms[it].nwl, dr[0], dr[1], dr[2], hrly);
        //         const double position = distance / delta_r;

        //         const double iq = static_cast<int>(position);
        //         const int ip = static_cast<int>(position);
        //         const double x0 = position - iq;
        //         const double x1 = 1.0 - x0;
        //         const double x2 = 2.0 - x0;
        //         const double x3 = 3.0 - x0;
        //         const double x12 = x1 * x2 / 6;
        //         const double x03 = x0 * x3 / 2;

        //         double tmp, dtmp, ddtmp;

        //         for (int iw = 0; iw < atom->nw; ++iw)
        //         {
        //             // this is a new 'l', we need 1D orbital wave
        //             // function from interpolation method.
        //             if (atom->iw2_new[iw])
        //             {
        //                 auto psi_uniform = it_psi_uniform[iw];
        //                 auto dpsi_uniform = it_dpsi_uniform[iw];
        //                 auto ddpsi_uniform = it_ddpsi_uniform[iw];

        //                 // if ( iq[id] >= philn.nr_uniform-4)
        //                 if (iq >= it_phi_nr_uniform[iw]-4)
        //                 {
        //                     tmp = dtmp = ddtmp = 0.0;
        //                 }
        //                 else
        //                 {
        //                     // use Polynomia Interpolation method to get the
        //                     // wave functions

        //                     tmp = x12 * (psi_uniform[ip] * x3 + psi_uniform[ip + 3] * x0)
        //                             + x03 * (psi_uniform[ip + 1] * x2 - psi_uniform[ip + 2] * x1);

        //                     dtmp = x12 * (dpsi_uniform[ip] * x3 + dpsi_uniform[ip + 3] * x0)
        //                             + x03 * (dpsi_uniform[ip + 1] * x2 - dpsi_uniform[ip + 2] * x1);

        //                     ddtmp = x12 * (ddpsi_uniform[ip] * x3 + ddpsi_uniform[ip + 3] * x0)
        //                             + x03 * (ddpsi_uniform[ip + 1] * x2 - ddpsi_uniform[ip + 2] * x1);
        //                 }
        //             } // new l is used.

        //             // get the 'l' of this localized wave function
        //             const int ll = atom->iw2l[iw];
        //             const int idx_lm = atom->iw2_ylm[iw];

        //             const double rl = pow_int(distance, ll);
        //             const double r_lp2 =rl * distance * distance;

        //             // d/dr (R_l / r^l)
        //             const double tmpdphi = (dtmp - tmp * ll / distance) / rl;
        //             const double term1 = ddtmp / r_lp2;
        //             const double term2 = (2 * ll + 1) * dtmp / r_lp2 / distance;
        //             const double term3 = ll * (ll + 2) * tmp / r_lp2 / distance / distance;
        //             const double term4 = tmpdphi / distance;
        //             const double term5 = term1 - term2 + term3;

        //             // hessian of (R_l / r^l)
        //             const double term_xx = term4 + dr[0] * dr[0] * term5;
        //             const double term_xy = dr[0] * dr[1] * term5;
        //             const double term_xz = dr[0] * dr[2] * term5;
        //             const double term_yy = term4 + dr[1] * dr[1] * term5;
        //             const double term_yz = dr[1] * dr[2] * term5;
        //             const double term_zz = term4 + dr[2] * dr[2] * term5;

        //             // d/dr (R_l / r^l) * alpha / r
        //             const double term_1x = dr[0] * term4;
        //             const double term_1y = dr[1] * term4;
        //             const double term_1z = dr[2] * term4;

        //             p_ddphi_xx[iw]
        //                 = term_xx * rly[idx_lm] + 2.0 * term_1x * grly[idx_lm][0] + tmp / rl * hrly[idx_lm][0];
        //             p_ddphi_xy[iw] = term_xy * rly[idx_lm] + term_1x * grly[idx_lm][1] + term_1y * grly[idx_lm][0]
        //                                 + tmp / rl * hrly[idx_lm][1];
        //             p_ddphi_xz[iw] = term_xz * rly[idx_lm] + term_1x * grly[idx_lm][2] + term_1z * grly[idx_lm][0]
        //                                 + tmp / rl * hrly[idx_lm][2];
        //             p_ddphi_yy[iw]
        //                 = term_yy * rly[idx_lm] + 2.0 * term_1y * grly[idx_lm][1] + tmp / rl * hrly[idx_lm][3];
        //             p_ddphi_yz[iw] = term_yz * rly[idx_lm] + term_1y * grly[idx_lm][2] + term_1z * grly[idx_lm][1]
        //                                 + tmp / rl * hrly[idx_lm][4];
        //             p_ddphi_zz[iw]
        //                 = term_zz * rly[idx_lm] + 2.0 * term_1z * grly[idx_lm][2] + tmp / rl * hrly[idx_lm][5];

        //         } // iw
        //     }     // end if
    }

    ModuleBase::timer::tick("GintAtom", "set_ddphi");
}

// explicit instantiation
template void GintAtom::set_ddphi(const std::vector<Vec3d>& coords, const int stride,
                                    double* ddphi_xx, double* ddphi_xy, double* ddphi_xz,
                                    double* ddphi_yy, double* ddphi_yz, double* ddphi_zz) const;

}