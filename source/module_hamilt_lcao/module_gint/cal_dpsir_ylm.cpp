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
    int it;
    double distance;
    double dr[3];
    double mt[3];

    Atom* atom;
    const UnitCell& ucell = *gt.ucell;
    std::vector<const double*> it_psi_uniform(gt.nwmax);
    std::vector<const double*> it_dpsi_uniform(gt.nwmax);;
    // array to store spherical harmonics and its derivatives
    // the first dimension equals 36 because the maximum nwl is 5.
    double rly[36];
    ModuleBase::Array_Pool<double> grly(36, 3);
    const int bcell_start = gt.bcell_start[grid_index];
    for (int id = 0; id < na_grid; id++)
    {
        const int mcell_index = bcell_start + id;
        get_grid_bigcell_distance(gt, mcell_index ,it, mt);
        Atom* atom = &ucell.atoms[it];
        get_psi_dpsi(gt, it, atom, it_psi_uniform, it_dpsi_uniform);

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

                ModuleBase::Ylm::grad_rl_sph_harm(atom->nwl,
                                                  dr[0], 
                                                  dr[1], 
                                                  dr[2], 
                                                  rly, 
                                                  grly.get_ptr_2D());

                dpsi_spline_interpolation(distance,
                                          dr,
                                          delta_r,
                                          atom,
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