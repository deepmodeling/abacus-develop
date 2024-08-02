#include "gint_tools.h"
#include "module_base/timer.h"
#include "module_base/ylm.h"
#include "module_base/array_pool.h"
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
    const int bcell_start = gt.bcell_start[grid_index];
    Atom* atom;
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
        const int mcell_index = bcell_start + id;
        const int imcell = gt.which_bigcell[mcell_index];
        int iat = gt.which_atom[mcell_index];
        const int it = ucell.iat2it[iat];
        const int ia = ucell.iat2ia[iat];
        
        const double mt[3] = {gt.meshball_positions[imcell][0] - gt.tau_in_bigcell[iat][0],
                              gt.meshball_positions[imcell][1] - gt.tau_in_bigcell[iat][1],
                              gt.meshball_positions[imcell][2] - gt.tau_in_bigcell[iat][2]};

        Atom* atom = &ucell.atoms[it];
        get_psi_dpsi(gt,atom->nw, it,
                atom->iw2_new,it_psi_uniform, it_dpsi_uniform);
        double distance;
        double dr[3];
        for (int ib = 0; ib < bxyz; ib++)
        {
            double* const p_psi = &psir_ylm[ib][block_index[id]];
            double* const p_dpsi_x = &dpsir_ylm_x[ib][block_index[id]];
            double* const p_dpsi_y = &dpsir_ylm_y[ib][block_index[id]];
            double* const p_dpsi_z = &dpsir_ylm_z[ib][block_index[id]];
            if (!cal_flag[ib][id])
            {
                ModuleBase::GlobalFunc::ZEROS(p_psi, block_size[id]);
                ModuleBase::GlobalFunc::ZEROS(p_dpsi_x, block_size[id]);
                ModuleBase::GlobalFunc::ZEROS(p_dpsi_y, block_size[id]);
                ModuleBase::GlobalFunc::ZEROS(p_dpsi_z, block_size[id]);
            }
            else
            {
                cal_grid_atom_distance(distance,dr,mt,gt.meshcell_pos[ib].data());

                ModuleBase::Ylm::grad_rl_sph_harm(ucell.atoms[it].nwl,
                                                  dr[0], dr[1], dr[2], 
                                                  rly, grly.get_ptr_2D());

                dpsi_spline_interpolation(distance,
                                          dr,
                                          delta_r,
                                          atom->nw,
                                          atom->iw2_new,
                                          atom->iw2l,
                                          atom->iw2_ylm,
                                          rly,
                                          grly.get_ptr_2D(),
                                          it_psi_uniform,
                                          it_dpsi_uniform,
                                          p_psi,
                                          p_dpsi_x,
                                          p_dpsi_y,
                                          p_dpsi_z);
            }     // else
        }
    }
    ModuleBase::timer::tick("Gint_Tools", "cal_dpsir_ylm");
    return;
}
}