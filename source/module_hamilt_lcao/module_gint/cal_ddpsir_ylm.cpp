#include "gint_tools.h"

namespace Gint_Tools{
void cal_ddpsir_ylm(
    const Grid_Technique& gt, const int bxyz,
    const int na_grid,                 // number of atoms on this grid
    const int grid_index,              // 1d index of FFT index (i,j,k)
    const double delta_r,              // delta_r of the uniform FFT grid
    const int* const block_index,      // block_index[na_grid+1], count total number of atomis orbitals
    const int* const block_size,       // block_size[na_grid],	number of columns of a band
    const bool* const* const cal_flag, // cal_flag[bxyz][na_grid],	whether the atom-grid distance is larger than cutoff
    double* const* const ddpsir_ylm_xx, double* const* const ddpsir_ylm_xy, double* const* const ddpsir_ylm_xz,
    double* const* const ddpsir_ylm_yy, double* const* const ddpsir_ylm_yz, double* const* const ddpsir_ylm_zz)
{
    ModuleBase::timer::tick("Gint_Tools", "cal_ddpsir_ylm");
    const UnitCell& ucell = *gt.ucell;
    std::vector<const double*> it_psi_uniform(gt.nwmax);
    std::vector<const double*> it_dpsi_uniform(gt.nwmax);
    std::vector<const double*> it_d2psi_uniform(gt.nwmax);
    std::vector<int> it_psi_nr_uniform(gt.nwmax);
    // array to store spherical harmonics and its derivatives
    // the first dimension equals 36 because the maximum nwl is 5.
    double rly[36];
    ModuleBase::Array_Pool<double> grly(36, 3);

    for (int id = 0; id < na_grid; id++)
    {
        const int mcell_index = gt.bcell_start[grid_index] + id;
        const int imcell = gt.which_bigcell[mcell_index];
        int iat = gt.which_atom[mcell_index];
        const int it = ucell.iat2it[iat];
        const int ia = ucell.iat2ia[iat];
        Atom* atom = &ucell.atoms[it];

        const double mt[3] = {gt.meshball_positions[imcell][0] - gt.tau_in_bigcell[iat][0],
                              gt.meshball_positions[imcell][1] - gt.tau_in_bigcell[iat][1],
                              gt.meshball_positions[imcell][2] - gt.tau_in_bigcell[iat][2]};

        for (int iw=0; iw< atom->nw; ++iw)
        {
            if ( atom->iw2_new[iw] )
            {
                it_psi_uniform[iw]= gt.psi_u[it*gt.nwmax + iw].data();
                it_dpsi_uniform[iw] = gt.dpsi_u[it*gt.nwmax + iw].data();
                it_psi_nr_uniform[iw]= gt.psi_u[it*gt.nwmax + iw].size();
            }
        }

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
               Gint_Func::ZEROS(p_ddpsi_xx, block_size[id]);
                Gint_Func::ZEROS(p_ddpsi_xy, block_size[id]);
                Gint_Func::ZEROS(p_ddpsi_xz, block_size[id]);
                Gint_Func::ZEROS(p_ddpsi_yy, block_size[id]);
                Gint_Func::ZEROS(p_ddpsi_yz, block_size[id]);
                Gint_Func::ZEROS(p_ddpsi_zz, block_size[id]);
            }
            else
            {
                const double dr[3]
                    = {// vectors between atom and grid
                       gt.meshcell_pos[ib][0] + mt[0], gt.meshcell_pos[ib][1] + mt[1], gt.meshcell_pos[ib][2] + mt[2]};
                double distance = std::sqrt(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]);

                // for some unknown reason, the finite difference between dpsi and ddpsi
                // using analytical expression is always wrong; as a result,
                // I switch to explicit finite difference method for evaluating
                // the second derivatives of the orbitals
                if (/*distance < 1e-9*/ true)
                {
                    double*** dpsi = new double**[atom->nw];
                    for (int i = 0; i < atom->nw; i++)
                    {
                        dpsi[i] = new double*[6];
                        for (int j = 0; j < 6; j++)
                        {
                            dpsi[i][j] = new double[3];
                            Gint_Func::ZEROS(dpsi[i][j], 3);
                        }
                    }

                    double* dr1 = new double[3];

                    double** displ = new double*[6];
                    for (int i = 0; i < 6; i++)
                    {
                        displ[i] = new double[3];
                        Gint_Func::ZEROS(displ[i], 3);
                    }
                    displ[0][0] = 0.0001; // in x direction
                    displ[1][0] = -0.0001;
                    displ[2][1] = 0.0001; // in y direction
                    displ[3][1] = -0.0001;
                    displ[4][2] = 0.0001; // in z direction
                    displ[5][2] = -0.0001;

                    for (int i = 0; i < 6; i++)
                    {
                        dr1[0] = dr[0] + displ[i][0];
                        dr1[1] = dr[1] + displ[i][1];
                        dr1[2] = dr[2] + displ[i][2];

                        ModuleBase::Ylm::grad_rl_sph_harm(ucell.atoms[it].nwl, dr1[0], dr1[1], dr1[2], rly, grly.get_ptr_2D());

                        double distance1 = std::sqrt(dr1[0] * dr1[0] + dr1[1] * dr1[1] + dr1[2] * dr1[2]);
                        if (distance1 < 1e-9) {
                            distance1 = 1e-9;
}

                        const double position = distance1 / delta_r;

                        const int ip = static_cast<int>(position);
                        const double iq = static_cast<int>(position);
                        const double x0 = position - iq;
                        const double x1 = 1.0 - x0;
                        const double x2 = 2.0 - x0;
                        const double x3 = 3.0 - x0;
                        const double x12 = x1 * x2 / 6;
                        const double x03 = x0 * x3 / 2;

                        double tmp, dtmp;

                        for (int iw = 0; iw < atom->nw; ++iw)
                        {
                            // this is a new 'l', we need 1D orbital wave
                            // function from interpolation method.
                            if (atom->iw2_new[iw])
                            {
                                auto psi_uniform = it_psi_uniform[iw];
                                auto dpsi_uniform = it_dpsi_uniform[iw];

                                // if ( iq[id] >= philn.nr_uniform-4)
                                if (iq >= it_psi_nr_uniform[iw]-4)
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
                            const int ll = atom->iw2l[iw];
                            const int idx_lm = atom->iw2_ylm[iw];

                            const double rl = pow(distance1, ll);

                            // derivative of wave functions with respect to atom positions.
                            const double tmpdphi_rly = (dtmp - tmp * ll / distance1) / rl * rly[idx_lm] / distance1;
                            const double tmprl = tmp / rl;

                            dpsi[iw][i][0] = tmpdphi_rly * dr1[0] + tmprl * grly[idx_lm][0];
                            dpsi[iw][i][1] = tmpdphi_rly * dr1[1] + tmprl * grly[idx_lm][1];
                            dpsi[iw][i][2] = tmpdphi_rly * dr1[2] + tmprl * grly[idx_lm][2];
                        } // end iw
                    }     // end i = 0-6

                    for (int iw = 0; iw < atom->nw; iw++)
                    {
                        p_ddpsi_xx[iw] = (dpsi[iw][0][0] - dpsi[iw][1][0]) / 0.0002;
                        p_ddpsi_xy[iw]
                            = ((dpsi[iw][2][0] - dpsi[iw][3][0]) + (dpsi[iw][0][1] - dpsi[iw][1][1])) / 0.0004;
                        p_ddpsi_xz[iw]
                            = ((dpsi[iw][4][0] - dpsi[iw][5][0]) + (dpsi[iw][0][2] - dpsi[iw][1][2])) / 0.0004;
                        p_ddpsi_yy[iw] = (dpsi[iw][2][1] - dpsi[iw][3][1]) / 0.0002;
                        p_ddpsi_yz[iw]
                            = ((dpsi[iw][4][1] - dpsi[iw][5][1]) + (dpsi[iw][2][2] - dpsi[iw][3][2])) / 0.0004;
                        p_ddpsi_zz[iw] = (dpsi[iw][4][2] - dpsi[iw][5][2]) / 0.0002;
                    }

                    for (int i = 0; i < atom->nw; i++)
                    {
                        for (int j = 0; j < 6; j++)
                        {
                            delete[] dpsi[i][j];
                        }
                        delete[] dpsi[i];
                    }
                    delete[] dpsi;

                    delete[] dr1;
                    for (int i = 0; i < 6; i++)
                    {
                        delete[] displ[i];
                    }
                    delete[] displ;
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
                    ModuleBase::Ylm::grad_rl_sph_harm(ucell.atoms[it].nwl, dr[0], dr[1], dr[2], rly, grly.get_ptr_2D());
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
                            if (iq >= it_psi_nr_uniform[iw]-4)
                            {
                                tmp = dtmp = ddtmp = 0.0;
                            }
                            else
                            {
                                // use Polynomia Interpolation method to get the
                                // wave functions

                                tmp = x12 * (psi_uniform[ip] * x3 + psi_uniform[ip + 3] * x0)
                                      + x03 * (psi_uniform[ip + 1] * x2 - psi_uniform[ip + 2] * x1);

                                dtmp = x12 * (dpsi_uniform[ip] * x3 + dpsi_uniform[ip + 3] * x0)
                                       + x03 * (dpsi_uniform[ip + 1] * x2 - dpsi_uniform[ip + 2] * x1);

                                ddtmp = x12 * (ddpsi_uniform[ip] * x3 + ddpsi_uniform[ip + 3] * x0)
                                        + x03 * (ddpsi_uniform[ip + 1] * x2 - ddpsi_uniform[ip + 2] * x1);
                            }
                        } // new l is used.

                        // get the 'l' of this localized wave function
                        const int ll = atom->iw2l[iw];
                        const int idx_lm = atom->iw2_ylm[iw];

                        const double rl = pow(distance, ll);
                        const double r_lp2 = pow(distance, ll + 2);

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

                        p_ddpsi_xx[iw]
                            = term_xx * rly[idx_lm] + 2.0 * term_1x * grly[idx_lm][0] + tmp / rl * hrly[idx_lm][0];
                        p_ddpsi_xy[iw] = term_xy * rly[idx_lm] + term_1x * grly[idx_lm][1] + term_1y * grly[idx_lm][0]
                                         + tmp / rl * hrly[idx_lm][1];
                        p_ddpsi_xz[iw] = term_xz * rly[idx_lm] + term_1x * grly[idx_lm][2] + term_1z * grly[idx_lm][0]
                                         + tmp / rl * hrly[idx_lm][2];
                        p_ddpsi_yy[iw]
                            = term_yy * rly[idx_lm] + 2.0 * term_1y * grly[idx_lm][1] + tmp / rl * hrly[idx_lm][3];
                        p_ddpsi_yz[iw] = term_yz * rly[idx_lm] + term_1y * grly[idx_lm][2] + term_1z * grly[idx_lm][1]
                                         + tmp / rl * hrly[idx_lm][4];
                        p_ddpsi_zz[iw]
                            = term_zz * rly[idx_lm] + 2.0 * term_1z * grly[idx_lm][2] + tmp / rl * hrly[idx_lm][5];

                    } // iw
                }     // end if
            }         // else
        }             // end ib
    }                 // end id(atom)
    ModuleBase::timer::tick("Gint_Tools", "cal_ddpsir_ylm");
    return;
}
}