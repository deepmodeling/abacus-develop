#include "gint_tools.h"
#include "module_base/global_function.h"
#include "module_base/timer.h"
#include "module_base/ylm.h"
namespace Gint_Tools{
void cal_psir_ylm(
    const Grid_Technique& gt, const int bxyz,
    const int na_grid,            // number of atoms on this grid
    const int grid_index,         // 1d index of FFT index (i,j,k)
    const double delta_r,         // delta_r of the uniform FFT grid
    const int* const block_index, // block_index[na_grid+1], count total number of atomis orbitals
    const int* const block_size,  // block_size[na_grid],	number of columns of a band
    const bool* const* const cal_flag,
    double* const* const psir_ylm) // cal_flag[bxyz][na_grid],	whether the atom-grid distance is larger than cutoff
{
    ModuleBase::timer::tick("Gint_Tools", "cal_psir_ylm");

    int it=0;
    double distance=0.0;
    std::array<double, 3> dr{0.0, 0.0, 0.0};
    std::array<double, 3> mt{0.0, 0.0, 0.0};

    Atom* atom;
    std::vector<double> ylma;
    std::vector<const double*> it_psi_uniform(gt.nwmax);
    std::vector<const double*> it_dpsi_uniform(gt.nwmax);
    const UnitCell& ucell = *gt.ucell;
    const int bcell_start = gt.bcell_start[grid_index];
    for (int id = 0; id < na_grid; id++)
    {
        // there are two parameters we want to know here:
        // in which bigcell of the meshball the atom is in?
        // what's the cartesian coordinate of the bigcell?
        // meshball_positions should be the bigcell position in meshball
        // to the center of meshball.
        // calculated in cartesian coordinates
        // the std::vector from the grid which is now being operated to the atom position.
        // in meshball language, is the std::vector from imcell to the center cel, plus
        // tau_in_bigcell.
        const int mcell_index = bcell_start + id;

        get_grid_bigcell_distance(gt,
                                 mcell_index ,
                                 it, 
                                 mt);

        atom = &ucell.atoms[it];
        get_psi_dpsi(gt, 
                     it, 
                     atom, 
                     it_psi_uniform, 
                     it_dpsi_uniform);

        // loop over the grids in the big cell
        for (int ib = 0; ib < bxyz; ib++)
        {
            double* p = &psir_ylm[ib][block_index[id]];
            if (!cal_flag[ib][id])
            {
                ModuleBase::GlobalFunc::ZEROS(p, block_size[id]);
            }
            else
            {
                cal_grid_atom_distance(distance,dr,mt,gt.meshcell_pos[ib].data());
                //------------------------------------------------------
                // spherical harmonic functions Ylm
                //------------------------------------------------------
                ModuleBase::Ylm::sph_harm(atom->nwl, 
                                            dr[0] / distance,
                                            dr[1] / distance, 
                                            dr[2] / distance,
                                            ylma);

                // these parameters are related to interpolation
                // because once the distance from atom to grid point is known,
                // we can obtain the parameters for interpolation and
                // store them first! these operations can save lots of efforts.

                spl_intrp(distance,
                            delta_r, 
                            atom, 
                            ylma,
                            it_psi_uniform, 
                            it_dpsi_uniform, 
                            p);

            }     
        }         // end ib
    }             // end id
    ModuleBase::timer::tick("Gint_Tools", "cal_psir_ylm");
    return;
}
}