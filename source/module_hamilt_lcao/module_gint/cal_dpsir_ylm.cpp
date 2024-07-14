#include "gint_tools.h"
namespace Gint_Tools{
void cal_dpsir_ylm(
    const Grid_Technique& gt, const int bxyz,
    const int na_grid,                 // number of atoms on this grid
    const int grid_index,              // 1d index of FFT index (i,j,k)
    const double delta_r,              // delta_r of the uniform FFT grid
    const int* const block_index,      // block_index[na_grid+1], count total number of atomis orbitals
    const int* const block_size,       // block_size[na_grid],	number of columns of a band
    const bool* const* const cal_flag, // cal_flag[bxyz][na_grid],	whether the atom-grid distance is larger than cutoff
    double* const* const psir_ylm, double* const* const dpsir_ylm_x, double* const* const dpsir_ylm_y,
    double* const* const dpsir_ylm_z)
{
    ModuleBase::timer::tick("Gint_Tools", "cal_dpsir_ylm");
    const UnitCell& ucell = *gt.ucell;
    std::vector<const double*> it_psi_uniform(gt.nwmax);
    std::vector<const double*> it_dpsi_uniform(gt.nwmax);
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
        // preprocess index
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
            double* const p_psi = &psir_ylm[ib][block_index[id]];
            double* const p_dpsi_x = &dpsir_ylm_x[ib][block_index[id]];
            double* const p_dpsi_y = &dpsir_ylm_y[ib][block_index[id]];
            double* const p_dpsi_z = &dpsir_ylm_z[ib][block_index[id]];
            if (!cal_flag[ib][id])
            {
                Gint_Func::ZEROS(p_psi, block_size[id]);
                Gint_Func::ZEROS(p_dpsi_x, block_size[id]);
                Gint_Func::ZEROS(p_dpsi_y, block_size[id]);
                Gint_Func::ZEROS(p_dpsi_z, block_size[id]);
            }
            else
            {
                const double dr[3]
                    = {// vectors between atom and grid
                       gt.meshcell_pos[ib][0] + mt[0], gt.meshcell_pos[ib][1] + mt[1], gt.meshcell_pos[ib][2] + mt[2]};
                double distance = std::sqrt(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]);

                ModuleBase::Ylm::grad_rl_sph_harm(ucell.atoms[it].nwl, dr[0], dr[1], dr[2], rly, grly.get_ptr_2D());
                if (distance < 1e-9)
                    distance = 1e-9;

                const double position = distance / delta_r;

                const double iq = static_cast<int>(position);
                const int ip = static_cast<int>(position);
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
                        if (iq >= it_psi_nr_uniform[iw] - 4)
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

                    const double rl = pow(distance, ll);

                    // 3D wave functions
                    p_psi[iw] = tmp * rly[idx_lm] / rl;

                    // derivative of wave functions with respect to atom positions.
                    const double tmpdphi_rly = (dtmp - tmp * ll / distance) / rl * rly[idx_lm] / distance;
                    const double tmprl = tmp / rl;

                    p_dpsi_x[iw] = tmpdphi_rly * dr[0] + tmprl * grly[idx_lm][0];
                    p_dpsi_y[iw] = tmpdphi_rly * dr[1] + tmprl * grly[idx_lm][1];
                    p_dpsi_z[iw] = tmpdphi_rly * dr[2] + tmprl * grly[idx_lm][2];
                } // iw
            }     // else
        }
    }
    ModuleBase::timer::tick("Gint_Tools", "cal_dpsir_ylm");
    return;
}
}