#include "gint_tools.h"
#include "module_base/global_function.h"
#include "module_base/timer.h"
#include "module_base/ylm.h"
namespace Gint_Tools{
void cal_dpsirr_ylm(
    const Grid_Technique& gt, 
	const int bxyz,
    const int na_grid,                 // number of atoms on this grid
    const int grid_index,              // 1d index of FFT index (i,j,k)
    const int* const block_index,      // block_index[na_grid+1], count total number of atomis orbitals
    const int* const block_size,       // block_size[na_grid],	number of columns of a band
    const bool* const* const cal_flag, // cal_flag[bxyz][na_grid],	whether the atom-grid distance is larger than cutoff
    double* const* const dpsir_ylm_x, 
	double* const* const dpsir_ylm_y,
	double* const* const dpsir_ylm_z,
    double* const* const dpsirr_ylm)
{
    ModuleBase::timer::tick("Gint_Tools", "cal_dpsirr_ylm");
	int it=0;
    double distance=0.0;
    std::array<double, 3> dr{0.0, 0.0, 0.0};
    std::array<double, 3> mt{0.0, 0.0, 0.0};

    Atom* atom;
    const UnitCell& ucell = *gt.ucell;
	const int bcell_start = gt.bcell_start[grid_index];
    for (int id = 0; id < na_grid; id++)
    {
        const int mcell_index = bcell_start + id;
		get_grid_bigcell_distance(gt, mcell_index ,it, mt);
        atom = &ucell.atoms[it];

		for(int ib=0; ib<bxyz; ib++)
		{
			double*const p_dpsi_x=&dpsir_ylm_x[ib][block_index[id]];
			double*const p_dpsi_y=&dpsir_ylm_y[ib][block_index[id]];
			double*const p_dpsi_z=&dpsir_ylm_z[ib][block_index[id]];
			double*const p_dpsirr=&dpsirr_ylm[ib][block_index[id] * 6];
			if(!cal_flag[ib][id])
			{
				ModuleBase::GlobalFunc::ZEROS(p_dpsirr, block_size[id] * 6);
			}
			else
			{
				cal_grid_atom_distance(distance,dr,mt,gt.meshcell_pos[ib].data());

				for (int iw=0; iw< atom->nw; ++iw)
				{
					p_dpsirr[iw * 6] = p_dpsi_x[iw]*dr[0];
					p_dpsirr[iw * 6 + 1] = p_dpsi_x[iw]*dr[1];
					p_dpsirr[iw * 6 + 2] = p_dpsi_x[iw]*dr[2];
					p_dpsirr[iw * 6 + 3] = p_dpsi_y[iw]*dr[1];
					p_dpsirr[iw * 6 + 4] = p_dpsi_y[iw]*dr[2];
					p_dpsirr[iw * 6 + 5] = p_dpsi_z[iw]*dr[2];
				}//iw
			}//else
		}
		}
		ModuleBase::timer::tick("Gint_Tools", "cal_dpsirr_ylm");
		return;
	}
}//namespace Gint_Tools