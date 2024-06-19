//=========================================================
//REFACTOR : Peize Lin, 2021.06.28
//=========================================================
#include "gint_tools.h"

#include <cmath>

#include "module_base/timer.h"
#include "module_base/ylm.h"
#include "module_basis/module_ao/ORB_read.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

namespace Gint_Tools{
int* get_vindex(const int bxyz, const int bx, const int by, const int bz, const int nplane, const int start_ind,
                const int ncyz)
{
    int* vindex = new int[bxyz];
    int bindex = 0;

		for(int ii=0; ii<bx; ii++)
		{
			const int ipart = ii*ncyz;
			for(int jj=0; jj<by; jj++)
			{
				const int jpart = jj*nplane + ipart;
				for(int kk=0; kk<bz; kk++)
				{
					vindex[bindex] = start_ind + kk + jpart;
					++bindex;
				}
			}
		}
		return vindex;
	}

	// here vindex refers to local potentials
    int* get_vindex(
        const int bxyz,
        const int bx,
        const int by,
        const int bz,
        const int nplane,
        const int ncyz,
		const int ibx,
		const int jby,
		const int kbz)
	{
		int *vindex = new int[bxyz];
		int bindex=0;
		// z is the fastest,
		// ipart can be obtained by using a previously stored array
		for(int ii=0; ii<bx; ii++)
		{
			const int ipart = (ibx+ii)*ncyz;
			for(int jj=0; jj<by; jj++)
			{
				// jpart can be obtained by using a previously stored array
				const int jpart = (jby+jj)*nplane + ipart;
				for(int kk=0; kk<bz; kk++)
				{
					vindex[bindex] = kbz+kk + jpart;
					++bindex;
				}
			}
		}
		return vindex;
	}

	// extract the local potentials.
	double* get_vldr3(
        const double* const vlocal,		// vlocal[ir]
        const int bxyz,
        const int bx,
        const int by,
        const int bz,
        const int nplane,
        const int ncyz,
		const int ibx,
		const int jby,
		const int kbz,
		const double dv)
	{
		// set the index for obtaining local potentials
		int* vindex = Gint_Tools::get_vindex(bxyz, bx, by, bz, nplane, ncyz, ibx, jby, kbz);
		double *vldr3 = new double[bxyz];
		for(int ib=0; ib<bxyz; ib++)
		{
			vldr3[ib]=vlocal[vindex[ib]] * dv;
		}
		delete[] vindex;
		return vldr3;
	}

	double* get_vldr3(
        const double* const vlocal,		// vlocal[ir]
        const int bxyz,
        const int bx,
        const int by,
        const int bz,
        const int nplane,
        const int start_ind,
		const int ncyz,
		const double dv)
	{
		// set the index for obtaining local potentials
		int* vindex = Gint_Tools::get_vindex(bxyz, bx, by, bz, nplane, start_ind, ncyz);
		double *vldr3 = new double[bxyz];
		for(int ib=0; ib<bxyz; ib++)
		{
			vldr3[ib]=vlocal[vindex[ib]] * dv;
		}
		delete[] vindex;
		return vldr3;
	}

void get_block_info(const Grid_Technique& gt, const int bxyz, const int na_grid, const int grid_index, int*& block_iw,
                    int*& block_index, int*& block_size, bool**& cal_flag)
{
    block_iw = new int[na_grid];
    block_index = new int[na_grid + 1];
    block_size = new int[na_grid];
    cal_flag = new bool*[bxyz];
    for (int ib = 0; ib < bxyz; ib++)
    {
        cal_flag[ib] = new bool[na_grid];
    }
    const UnitCell& ucell = *gt.ucell;
    block_index[0] = 0;
    for (int id = 0; id < na_grid; id++)
    {
        const int mcell_index = gt.bcell_start[grid_index] + id;
        const int iat = gt.which_atom[mcell_index];    // index of atom
        const int it = ucell.iat2it[iat];              // index of atom type
        const int ia = ucell.iat2ia[iat];              // index of atoms within each type
        const int start = ucell.itiaiw2iwt(it, ia, 0); // the index of the first wave function for atom (it,ia)
        block_iw[id] = gt.trace_lo[start];
        block_index[id + 1] = block_index[id] + ucell.atoms[it].nw;
        block_size[id] = ucell.atoms[it].nw;

			const int imcell=gt.which_bigcell[mcell_index];
			const double mt[3] = {
				gt.meshball_positions[imcell][0] - gt.tau_in_bigcell[iat][0],
				gt.meshball_positions[imcell][1] - gt.tau_in_bigcell[iat][1],
				gt.meshball_positions[imcell][2] - gt.tau_in_bigcell[iat][2]};

			for(int ib=0; ib<bxyz; ib++)
			{
				// meshcell_pos: z is the fastest
				const double dr[3] = {
					gt.meshcell_pos[ib][0] + mt[0],
					gt.meshcell_pos[ib][1] + mt[1],
					gt.meshcell_pos[ib][2] + mt[2]};
				const double distance = std::sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);	// distance between atom and grid

            if (distance > gt.rcuts[it] - 1.0e-10)
                cal_flag[ib][id] = false;
            else
                cal_flag[ib][id] = true;
        } // end ib
    }
}

void cal_dpsirr_ylm(
    const Grid_Technique& gt, const int bxyz,
    const int na_grid,                 // number of atoms on this grid
    const int grid_index,              // 1d index of FFT index (i,j,k)
    const int* const block_index,      // block_index[na_grid+1], count total number of atomis orbitals
    const int* const block_size,       // block_size[na_grid],	number of columns of a band
    const bool* const* const cal_flag, // cal_flag[bxyz][na_grid],	whether the atom-grid distance is larger than cutoff
    double* const* const dpsir_ylm_x, double* const* const dpsir_ylm_y, double* const* const dpsir_ylm_z,
    double* const* const dpsir_ylm_xx, double* const* const dpsir_ylm_xy, double* const* const dpsir_ylm_xz,
    double* const* const dpsir_ylm_yy, double* const* const dpsir_ylm_yz, double* const* const dpsir_ylm_zz)
{
    ModuleBase::timer::tick("Gint_Tools", "cal_dpsirr_ylm");
    const UnitCell& ucell = *gt.ucell;
    for (int id = 0; id < na_grid; id++)
    {
        const int mcell_index = gt.bcell_start[grid_index] + id;
        const int imcell = gt.which_bigcell[mcell_index];
        int iat = gt.which_atom[mcell_index];
        const int it = ucell.iat2it[iat];
        Atom* atom = &ucell.atoms[it];

			const double mt[3]={
				gt.meshball_positions[imcell][0] - gt.tau_in_bigcell[iat][0],
				gt.meshball_positions[imcell][1] - gt.tau_in_bigcell[iat][1],
				gt.meshball_positions[imcell][2] - gt.tau_in_bigcell[iat][2]};

			for(int ib=0; ib<bxyz; ib++)
			{
				double*const p_dpsi_x=&dpsir_ylm_x[ib][block_index[id]];
				double*const p_dpsi_y=&dpsir_ylm_y[ib][block_index[id]];
				double*const p_dpsi_z=&dpsir_ylm_z[ib][block_index[id]];
				double*const p_dpsi_xx=&dpsir_ylm_xx[ib][block_index[id]];
				double*const p_dpsi_xy=&dpsir_ylm_xy[ib][block_index[id]];
				double*const p_dpsi_xz=&dpsir_ylm_xz[ib][block_index[id]];
				double*const p_dpsi_yy=&dpsir_ylm_yy[ib][block_index[id]];
				double*const p_dpsi_yz=&dpsir_ylm_yz[ib][block_index[id]];
				double*const p_dpsi_zz=&dpsir_ylm_zz[ib][block_index[id]];
				if(!cal_flag[ib][id])
				{
					ModuleBase::GlobalFunc::ZEROS(p_dpsi_xx, block_size[id]);
					ModuleBase::GlobalFunc::ZEROS(p_dpsi_xy, block_size[id]);
					ModuleBase::GlobalFunc::ZEROS(p_dpsi_xz, block_size[id]);
					ModuleBase::GlobalFunc::ZEROS(p_dpsi_yy, block_size[id]);
					ModuleBase::GlobalFunc::ZEROS(p_dpsi_yz, block_size[id]);
					ModuleBase::GlobalFunc::ZEROS(p_dpsi_zz, block_size[id]);
				}
				else
				{
					const double dr[3]={						// vectors between atom and grid
						gt.meshcell_pos[ib][0] + mt[0],
						gt.meshcell_pos[ib][1] + mt[1],
						gt.meshcell_pos[ib][2] + mt[2]};

					for (int iw=0; iw< atom->nw; ++iw)
					{

						p_dpsi_xx[iw] = p_dpsi_x[iw]*dr[0];
						p_dpsi_xy[iw] = p_dpsi_x[iw]*dr[1];
						p_dpsi_xz[iw] = p_dpsi_x[iw]*dr[2];
						p_dpsi_yy[iw] = p_dpsi_y[iw]*dr[1];
						p_dpsi_yz[iw] = p_dpsi_y[iw]*dr[2];
						p_dpsi_zz[iw] = p_dpsi_z[iw]*dr[2];

					}//iw
				}//else
			}
		}
		ModuleBase::timer::tick("Gint_Tools", "cal_dpsirr_ylm");
		return;
	}

	// atomic basis sets
	// psir_vlbr3[bxyz][LD_pool]
    Gint_Tools::Array_Pool<double> get_psir_vlbr3(
        const int bxyz,
        const int na_grid,  					    // how many atoms on this (i,j,k) grid
		const int LD_pool,
		const int*const block_index,		    	// block_index[na_grid+1], count total number of atomis orbitals
		const bool*const*const cal_flag,	    	// cal_flag[bxyz][na_grid],	whether the atom-grid distance is larger than cutoff
		const double*const vldr3,			    	// vldr3[bxyz]
		const double*const*const psir_ylm)		    // psir_ylm[bxyz][LD_pool]
	{
		Gint_Tools::Array_Pool<double> psir_vlbr3(bxyz, LD_pool);
		for(int ib=0; ib<bxyz; ++ib)
		{
			for(int ia=0; ia<na_grid; ++ia)
			{
				if(cal_flag[ib][ia])
				{
					for(int i=block_index[ia]; i<block_index[ia+1]; ++i)
					{
						psir_vlbr3.ptr_2D[ib][i]=psir_ylm[ib][i]*vldr3[ib];
					}
				}
				else
				{
					for(int i=block_index[ia]; i<block_index[ia+1]; ++i)
					{
						psir_vlbr3.ptr_2D[ib][i]=0;
					}
				}

			}
		}
		return psir_vlbr3;
	}


// calculating (psi_DMR)_mu = sum_nu DMR_mu,nu psi_nu
// note : there is a difference between rho and force
// in calculating rho, due to symmetry, the summation over mu,nu
// can be done as sum_mu,mu + 2 sum_mu<nu, saving some time
// but for force, we cannot exchange the index

} // namespace Gint_Tools