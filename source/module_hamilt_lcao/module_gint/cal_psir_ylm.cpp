#include "gint_tools.h"
#include "module_base/timer.h"
#include "module_base/ylm.h"
namespace Gint_Tools{
void cal_psir_ylm(
    const Grid_Technique& gt,
    const int bxyz,
    const int na_grid,            // number of atoms on this grid
    const int grid_index,         // 1d index of FFT index (i,j,k)
    const double delta_r,         // delta_r of the uniform FFT grid
    const int* const block_index, // block_index[na_grid+1], count total number of atomis orbitals
    const int* const block_size,  // block_size[na_grid],	number of columns of a band
    const bool* const* const cal_flag,
    double* const* const psir_ylm) // cal_flag[bxyz][na_grid],	whether the atom-grid distance is larger than cutoff
{
    ModuleBase::timer::tick("Gint_Tools", "cal_psir_ylm");
    std::vector<double> ylma;
    const UnitCell& ucell = *gt.ucell;
    std::vector<const double*> it_psi_uniform(gt.nwmax);
    std::vector<const double*> it_dpsi_uniform(gt.nwmax);

    for (int id = 0; id < na_grid; id++)
    {
        // there are two parameters we want to know here:
        // in which bigcell of the meshball the atom is in?
        // what's the cartesian coordinate of the bigcell?
        const int mcell_index = gt.bcell_start[grid_index] + id;

        const int iat = gt.which_atom[mcell_index]; // index of atom
        const int it = ucell.iat2it[iat];           // index of atom type
        const Atom* const atom = &ucell.atoms[it];
        std::vector<const double*> it_psi_uniform(atom->nw);
        std::vector<const double*> it_dpsi_uniform(atom->nw);
        // preprocess index
        for (int iw = 0; iw < atom->nw; ++iw)
        {
            if (atom->iw2_new[iw])
            {
                it_psi_uniform[iw]= gt.psi_u[it*gt.nwmax + iw].data();
                it_dpsi_uniform[iw] = gt.dpsi_u[it*gt.nwmax + iw].data();
            }
        }

        // meshball_positions should be the bigcell position in meshball
        // to the center of meshball.
        // calculated in cartesian coordinates
        // the std::vector from the grid which is now being operated to the atom position.
        // in meshball language, is the std::vector from imcell to the center cel, plus
        // tau_in_bigcell.
        const int imcell = gt.which_bigcell[mcell_index];
        const double mt[3] = {gt.meshball_positions[imcell][0] - gt.tau_in_bigcell[iat][0],
                              gt.meshball_positions[imcell][1] - gt.tau_in_bigcell[iat][1],
                              gt.meshball_positions[imcell][2] - gt.tau_in_bigcell[iat][2]};

        // number of grids in each big cell (bxyz)
        for (int ib = 0; ib < bxyz; ib++)
        {
            double* p = &psir_ylm[ib][block_index[id]];
            if (!cal_flag[ib][id])
            {
                ModuleBase::GlobalFunc::ZEROS(p, block_size[id]);
            }
            else
            {
                // meshcell_pos: z is the fastest
                const double dr[3]
                    = {gt.meshcell_pos[ib][0] + mt[0], gt.meshcell_pos[ib][1] + mt[1], gt.meshcell_pos[ib][2] + mt[2]};
                double distance
                    = std::sqrt(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]); // distance between atom and grid
                // if(distance[id] > gt.orbital_rmax) continue;
                if (distance < 1.0E-9)
                    distance += 1.0E-9;

                //------------------------------------------------------
                // spherical harmonic functions Ylm
                //------------------------------------------------------
                //	Ylm::get_ylm_real(this->nnn[it], this->dr[id], ylma);
                ModuleBase::Ylm::sph_harm(ucell.atoms[it].nwl, dr[0] / distance, dr[1] / distance, dr[2] / distance,
                                          ylma);
                // these parameters are related to interpolation
                // because once the distance from atom to grid point is known,
                // we can obtain the parameters for interpolation and
                // store them first! these operations can save lots of efforts.
                const double position = distance / delta_r;
                const int ip = static_cast<int>(position);
                const double dx = position - ip;
                const double dx2 = dx * dx;
                const double dx3 = dx2 * dx;

                const double c3 = 3.0 * dx2 - 2.0 * dx3;
                const double c1 = 1.0 - c3;
                const double c2 = (dx - 2.0 * dx2 + dx3) * delta_r;
                const double c4 = (dx3 - dx2) * delta_r;

                double phi = 0;
                for (int iw = 0; iw < atom->nw; ++iw)
                {
                    if (atom->iw2_new[iw])
                    {
                        auto psi_uniform = it_psi_uniform[iw];
                        auto dpsi_uniform = it_dpsi_uniform[iw];
                        phi = c1 * psi_uniform[ip] + c2 * dpsi_uniform[ip] // radial wave functions
                              + c3 * psi_uniform[ip + 1] + c4 * dpsi_uniform[ip + 1];
                    }
                    p[iw] = phi * ylma[atom->iw2_ylm[iw]];
                } // end iw
            }     // end distance<=(rcuts[it]-1.0e-15)
        }         // end ib
    }             // end id
    ModuleBase::timer::tick("Gint_Tools", "cal_psir_ylm");
    return;
}
}