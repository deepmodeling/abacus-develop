#include "gint_tools.h"
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

    int it;
    int ip;
    double distance;
    double dr[3];
    double mt[3];

    Atom* atom;
    std::vector<double> coeffs(4,0.0);
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
	const int iat = gt.which_atom[mcell_index]; 
	const int it = gt.ucell->iat2it[iat];  
	const int imcell = gt.which_bigcell[mcell_index];
	const double mt[3] = {gt.meshball_positions[imcell][0] - gt.tau_in_bigcell[iat][0],
							gt.meshball_positions[imcell][1] - gt.tau_in_bigcell[iat][1],
							gt.meshball_positions[imcell][2] - gt.tau_in_bigcell[iat][2]};
        // get_grid_bigcell_distance(gt, bcell_start, id ,it, mt);

        atom = &ucell.atoms[it];
        get_psi_dpsi(gt, atom->nw, it, atom->iw2_new, it_psi_uniform, it_dpsi_uniform);

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

                interp_coeff(distance, delta_r, ip, coeffs.data());
                // these parameters are related to interpolation
                // because once the distance from atom to grid point is known,
                // we can obtain the parameters for interpolation and
                // store them first! these operations can save lots of efforts.
                double phi = 0;
                for (int iw = 0; iw < atom->nw; ++iw)
                {
                    if (atom->iw2_new[iw])
                    {
                        auto psi_uniform = it_psi_uniform[iw];
                        auto dpsi_uniform = it_dpsi_uniform[iw];
                        phi = coeffs[0] * psi_uniform[ip] + coeffs[1] * dpsi_uniform[ip] // radial wave functions
                                + coeffs[2] * psi_uniform[ip + 1] + coeffs[3] * dpsi_uniform[ip + 1];
                    }
                    p[iw] = phi * ylma[atom->iw2_ylm[iw]];
                } // end iw

            }     
        }         // end ib
    }             // end id
    ModuleBase::timer::tick("Gint_Tools", "cal_psir_ylm");
    return;
}
}