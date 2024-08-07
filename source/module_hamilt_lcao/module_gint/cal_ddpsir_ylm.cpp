#include "gint_tools.h"
#include "module_base/global_function.h"
#include "module_base/timer.h"
#include "module_base/ylm.h"
namespace Gint_Tools{
void cal_ddpsir_ylm(
    const Grid_Technique& gt, 
    const int bxyz,
    const int na_grid,                 // number of atoms on this grid
    const int grid_index,              // 1d index of FFT index (i,j,k)
    const double delta_r,              // delta_r of the uniform FFT grid
    const int* const block_index,      // block_index[na_grid+1], count total number of atomis orbitals
    const int* const block_size,       // block_size[na_grid],	number of columns of a band
    const bool* const* const cal_flag, // cal_flag[bxyz][na_grid],	whether the atom-grid distance is larger than cutoff
    double* const* const ddpsir_ylm_xx, 
    double* const* const ddpsir_ylm_xy, 
    double* const* const ddpsir_ylm_xz,
    double* const* const ddpsir_ylm_yy, 
    double* const* const ddpsir_ylm_yz, 
    double* const* const ddpsir_ylm_zz)
{
    ModuleBase::timer::tick("Gint_Tools", "cal_ddpsir_ylm");

    int it=0;
    double distance1=0.0;
    double distance=0.0;
    std::array<double,3> dr={0.0,0.0,0.0};
    std::array<double,3> mt={0.0,0.0,0.0};
    std::array<double,3> dr1={0.0,0.0,0.0};
                    
    Atom* atom;
    const UnitCell& ucell = *gt.ucell;
    std::vector<const double*> it_psi_uniform(gt.nwmax);
    std::vector<const double*> it_dpsi_uniform(gt.nwmax);
    std::vector<const double*> it_d2psi_uniform(gt.nwmax);

    // array to store spherical harmonics and its derivatives
    // the first dimension equals 36 because the maximum nwl is 5.
    std::array<double, 36> rly;
    ModuleBase::Array_Pool<double> grly(36, 3);
    ModuleBase::Array_Pool<std::array<double,3>> dpsi(gt.nwmax,6);
    ModuleBase::Array_Pool<double> displ(6,3);
    displ[0][0] = 0.0001; 
    displ[1][0] = -0.0001; // in x direction
    displ[2][1] = 0.0001; 
    displ[3][1] = -0.0001; // in y direction
    displ[4][2] = 0.0001; 
    displ[5][2] = -0.0001; // in z direction

    const int bcell_start = gt.bcell_start[grid_index];
    for (int id = 0; id < na_grid; id++)
    {
        const int mcell_index = bcell_start + id;
        get_grid_bigcell_distance(gt, mcell_index ,it, mt);
        Atom* atom = &ucell.atoms[it];
        get_psi_dpsi(gt, it, atom, it_psi_uniform, it_dpsi_uniform);
                
        for (int ib = 0; ib < bxyz; ib++)
        {
            double* const p_ddpsi_xx = &ddpsir_ylm_xx[ib][block_index[id]];
            double* const p_ddpsi_xy = &ddpsir_ylm_xy[ib][block_index[id]];
            double* const p_ddpsi_xz = &ddpsir_ylm_xz[ib][block_index[id]];
            double* const p_ddpsi_yy = &ddpsir_ylm_yy[ib][block_index[id]];
            double* const p_ddpsi_yz = &ddpsir_ylm_yz[ib][block_index[id]];
            double* const p_ddpsi_zz = &ddpsir_ylm_zz[ib][block_index[id]];
            if (!cal_flag[ib][id])
            {
                ModuleBase::GlobalFunc::ZEROS(p_ddpsi_xx, block_size[id]);
                ModuleBase::GlobalFunc::ZEROS(p_ddpsi_xy, block_size[id]);
                ModuleBase::GlobalFunc::ZEROS(p_ddpsi_xz, block_size[id]);
                ModuleBase::GlobalFunc::ZEROS(p_ddpsi_yy, block_size[id]);
                ModuleBase::GlobalFunc::ZEROS(p_ddpsi_yz, block_size[id]);
                ModuleBase::GlobalFunc::ZEROS(p_ddpsi_zz, block_size[id]);
            }
            else
            {

                cal_grid_atom_distance(distance,dr,mt,gt.meshcell_pos[ib].data());
                // for some unknown reason, the finite difference between dpsi and ddpsi
                // using analytical expression is always wrong; as a result,
                // I switch to explicit finite difference method for evaluating
                // the second derivatives of the orbitals
                if (true)
                {
                    for (int i = 0; i < 6; i++)
                    {
                        cal_grid_atom_distance(distance1,
                                                dr1,
                                                dr,
                                                displ.get_ptr_2D()[i]);

                        ModuleBase::Ylm::grad_rl_sph_harm(atom->nwl, 
                                                          dr1[0],
                                                          dr1[1], 
                                                          dr1[2],
                                                          rly.data(), 
                                                          grly.get_ptr_2D());

                        dpsi_spl_intrp(distance1,
                                        dr1,
                                        delta_r,
                                        i,
                                        atom,
                                        rly.data(),
                                        grly.get_ptr_2D(),
                                        it_psi_uniform,
                                        it_dpsi_uniform,
                                        dpsi);
                    }

                    for (int iw = 0; iw < atom->nw; iw++)
                    {
                        p_ddpsi_xx[iw] = (dpsi[iw][0][0] - dpsi[iw][1][0]) / 0.0002;
                        p_ddpsi_xy[iw] = ((dpsi[iw][2][0] - dpsi[iw][3][0]) + 
                                          (dpsi[iw][0][1] - dpsi[iw][1][1])) / 0.0004;
                        p_ddpsi_xz[iw] = ((dpsi[iw][4][0] - dpsi[iw][5][0]) + 
                                          (dpsi[iw][0][2] - dpsi[iw][1][2])) / 0.0004;
                        p_ddpsi_yy[iw] = (dpsi[iw][2][1] - dpsi[iw][3][1]) / 0.0002;
                        p_ddpsi_yz[iw] = ((dpsi[iw][4][1] - dpsi[iw][5][1]) + 
                                          (dpsi[iw][2][2] - dpsi[iw][3][2])) / 0.0004;
                        p_ddpsi_zz[iw] = (dpsi[iw][4][2] - dpsi[iw][5][2]) / 0.0002;
                    }
                }
                else
                // the analytical method for evaluating 2nd derivatives
                // it is not used currently
                {
                    // Add it here, but do not run it. If there is a need to run this code 
                    // in the future, include it in the previous initialization process.
                    for (int iw=0; iw< atom->nw; ++iw)
                    {
                        if ( atom->iw2_new[iw] )
                        {
                            it_d2psi_uniform[iw] = gt.d2psi_u[it*gt.nwmax + iw].data();
                        }
                    }
                    // End of code addition section.

                    std::vector<std::vector<double>> hrly;
                    ModuleBase::Ylm::grad_rl_sph_harm(ucell.atoms[it].nwl, dr[0], dr[1], dr[2], rly.data(), grly.get_ptr_2D());
                    ModuleBase::Ylm::hes_rl_sph_harm(ucell.atoms[it].nwl, dr[0], dr[1], dr[2], hrly);
                    const double position = distance / delta_r;

                    const double iq = static_cast<int>(position);
                    const int ip = static_cast<int>(position);
                    const double x0 = position - iq;
                    const double x1 = 1.0 - x0;
                    const double x2 = 2.0 - x0;
                    const double x3 = 3.0 - x0;
                    const double x12 = x1 * x2 / 6;
                    const double x03 = x0 * x3 / 2;

                    double tmp, dtmp, ddtmp;

                    for (int iw = 0; iw < atom->nw; ++iw)
                    {
                        // this is a new 'l', we need 1D orbital wave
                        // function from interpolation method.
                        if (atom->iw2_new[iw])
                        {
                            auto psi_uniform = it_psi_uniform[iw];
                            auto dpsi_uniform = it_dpsi_uniform[iw];
                            auto ddpsi_uniform = it_d2psi_uniform[iw];

                            // if ( iq[id] >= philn.nr_uniform-4)
                            // if (iq >= it_psi_nr_uniform[iw]-4)
                            // {
                            //     tmp = dtmp = ddtmp = 0.0;
                            // }
                            // else
                            // {
                                // use Polynomia Interpolation method to get the
                                // wave functions

                                tmp = x12 * (psi_uniform[ip] * x3 + psi_uniform[ip + 3] * x0)
                                      + x03 * (psi_uniform[ip + 1] * x2 - psi_uniform[ip + 2] * x1);

                                dtmp = x12 * (dpsi_uniform[ip] * x3 + dpsi_uniform[ip + 3] * x0)
                                       + x03 * (dpsi_uniform[ip + 1] * x2 - dpsi_uniform[ip + 2] * x1);

                                ddtmp = x12 * (ddpsi_uniform[ip] * x3 + ddpsi_uniform[ip + 3] * x0)
                                        + x03 * (ddpsi_uniform[ip + 1] * x2 - ddpsi_uniform[ip + 2] * x1);
                            // }
                        } // new l is used.

                        // get the 'l' of this localized wave function
                        const int ll = atom->iw2l[iw];
                        const int idx_lm = atom->iw2_ylm[iw];

                        const double rl = pow_int(distance, ll);
                        const double r_lp2 =rl * distance * distance;

                        // d/dr (R_l / r^l)
                        const double tmpdphi = (dtmp - tmp * ll / distance) / rl;
                        const double term1 = ddtmp / r_lp2;
                        const double term2 = (2 * ll + 1) * dtmp / r_lp2 / distance;
                        const double term3 = ll * (ll + 2) * tmp / r_lp2 / distance / distance;
                        const double term4 = tmpdphi / distance;
                        const double term5 = term1 - term2 + term3;

                        // hessian of (R_l / r^l)
                        const double term_xx = term4 + dr[0] * dr[0] * term5;
                        const double term_xy = dr[0] * dr[1] * term5;
                        const double term_xz = dr[0] * dr[2] * term5;
                        const double term_yy = term4 + dr[1] * dr[1] * term5;
                        const double term_yz = dr[1] * dr[2] * term5;
                        const double term_zz = term4 + dr[2] * dr[2] * term5;

                        // d/dr (R_l / r^l) * alpha / r
                        const double term_1x = dr[0] * term4;
                        const double term_1y = dr[1] * term4;
                        const double term_1z = dr[2] * term4;

                        p_ddpsi_xx[iw] = term_xx * rly[idx_lm] + 2.0 * term_1x * grly[idx_lm][0] + 
                                        tmp / rl * hrly[idx_lm][0];
                        p_ddpsi_xy[iw] = term_xy * rly[idx_lm] + term_1x * grly[idx_lm][1] + 
                                        term_1y * grly[idx_lm][0] + tmp / rl * hrly[idx_lm][1];
                        p_ddpsi_xz[iw] = term_xz * rly[idx_lm] + term_1x * grly[idx_lm][2] + 
                                        term_1z * grly[idx_lm][0] + tmp / rl * hrly[idx_lm][2];
                        p_ddpsi_yy[iw] = term_yy * rly[idx_lm] + 2.0 * term_1y * grly[idx_lm][1] + 
                                        tmp / rl * hrly[idx_lm][3];
                        p_ddpsi_yz[iw] = term_yz * rly[idx_lm] + term_1y * grly[idx_lm][2] + 
                                        term_1z * grly[idx_lm][1] + tmp / rl * hrly[idx_lm][4];
                        p_ddpsi_zz[iw] = term_zz * rly[idx_lm] + 2.0 * term_1z * grly[idx_lm][2] +
                                        tmp / rl * hrly[idx_lm][5];

                    } // iw
                }     // end if
            }         // else
        }             // end ib

    }                 // end id(atom)
    ModuleBase::timer::tick("Gint_Tools", "cal_ddpsir_ylm");
    return;
}
}