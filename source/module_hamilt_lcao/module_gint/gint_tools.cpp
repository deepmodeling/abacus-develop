//=========================================================
//REFACTOR : Peize Lin, 2021.06.28
//=========================================================
#include "gint_tools.h"

#include <cmath>

#include "module_base/timer.h"
#include "module_base/ylm.h"
#include "module_base/array_pool.h"
#include "module_basis/module_ao/ORB_read.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

namespace Gint_Tools{
void get_vindex(const int bxyz, const int bx, const int by, const int bz, 
				const int nplane, const int start_ind,
				const int ncyz,int* vindex)
{
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
	}

	// here vindex refers to local potentials

	// extract the local potentials.
	void get_gint_vldr3(
		double* vldr3,
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
		std::vector<int> vindex(bxyz,0);
		Gint_Tools::get_vindex(bxyz, bx, by, bz, nplane, start_ind, ncyz,vindex.data());
		for(int ib=0; ib<bxyz; ib++)
		{
			vldr3[ib]=vlocal[vindex[ib]] * dv;
		}
	}

	void get_block_info(const Grid_Technique& gt, const int bxyz, const int na_grid, const int grid_index, int* block_iw,
						int* block_index, int* block_size, bool** cal_flag)
	{
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

			if (distance > gt.rcuts[it] - 1.0e-10) {
				cal_flag[ib][id] = false;
			} else {
				cal_flag[ib][id] = true;
			}
			} // end ib
		}
	}

void get_grid_bigcell_distance(const Grid_Technique& gt,
								const int mcell_index,
								int& it,
								std::array<double, 3>& mt)
{
	const int iat = gt.which_atom[mcell_index]; 
	const int imcell = gt.which_bigcell[mcell_index];
	it = gt.ucell->iat2it[iat];  
	mt[0] = gt.meshball_positions[imcell][0] - gt.tau_in_bigcell[iat][0];
	mt[1] =	gt.meshball_positions[imcell][1] - gt.tau_in_bigcell[iat][1];
	mt[2] =	gt.meshball_positions[imcell][2] - gt.tau_in_bigcell[iat][2];
}

void cal_grid_atom_distance(double &distance,
                            std::array<double, 3>& dr,
                            std::array<double, 3>& mt,
                            const double* meshcell_pos)
{
	dr[0] = meshcell_pos[0] + mt[0];
	dr[1] = meshcell_pos[1] + mt[1];
	dr[2] = meshcell_pos[2] + mt[2];
	distance = std::sqrt(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]); 
	if (distance < 1.0E-9) { distance += 1.0E-9;
}
}


void spl_intrp(const double distance,
							const double delta_r,
							Atom*& atom,
							std::vector<double>& ylma,
							std::vector<const double*>& it_psi_uniform,
							std::vector<const double*>& it_dpsi_uniform,
							double *p)
{
	double coeffs[4] = {0.0};
	const double position = distance / delta_r;
	int ip = static_cast<int>(position);
	const double dx = position - ip;
	const double dx2 = dx * dx;
	const double dx3 = dx2 * dx;
	coeffs[2] = 3.0 * dx2 - 2.0 * dx3;
	coeffs[0] = 1.0 - coeffs[2];
	coeffs[1] = (dx - 2.0 * dx2 + dx3) * delta_r;
	coeffs[3] = (dx3 - dx2) * delta_r;

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

void dpsi_spl_intrp(const double distance,
					std::array<double,3>& dr,
					const double delta_r,
					Atom*& atom,
					double* rly,
					double** grly,
					std::vector<const double*>& it_psi_uniform,
					std::vector<const double*>& it_dpsi_uniform,
					double *p_psi,
					double *p_dpsi_x,
					double *p_dpsi_y,
					double *p_dpsi_z)
{
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
			// use Polynomia Interpolation method to get the
			// wave functions

			tmp = x12 * (psi_uniform[ip] * x3 + psi_uniform[ip + 3] * x0)
					+ x03 * (psi_uniform[ip + 1] * x2 - psi_uniform[ip + 2] * x1);

			dtmp = x12 * (dpsi_uniform[ip] * x3 + dpsi_uniform[ip + 3] * x0)
					+ x03 * (dpsi_uniform[ip + 1] * x2 - dpsi_uniform[ip + 2] * x1);
                        // }
		} // new l is used.

		// get the 'l' of this localized wave function
			const int ll = atom->iw2l[iw];
			const int idx_lm = atom->iw2_ylm[iw];

			const double rl = pow_int(distance, ll);
			const double tmprl = tmp / rl;
                    
			// 3D wave functions
			p_psi[iw] = tmprl * rly[idx_lm];

		// derivative of wave functions with respect to atom positions.
			const double tmpdphi_rly = (dtmp - tmp * ll / distance) / rl * rly[idx_lm] / distance;

			p_dpsi_x[iw] = tmpdphi_rly * dr[0] + tmprl * grly[idx_lm][0];
			p_dpsi_y[iw] = tmpdphi_rly * dr[1] + tmprl * grly[idx_lm][1];
			p_dpsi_z[iw] = tmpdphi_rly * dr[2] + tmprl * grly[idx_lm][2];
		} // iw
}

void dpsi_spl_intrp(const double distance1,
					std::array<double,3>& dr1,
					const double delta_r,
					const int i,
					Atom*& atom,
					double* rly,
					double** grly,
					std::vector<const double*>& it_psi_uniform,
					std::vector<const double*>& it_dpsi_uniform,
					ModuleBase::Array_Pool<std::array<double,3>>& dpsi)
{
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
					// use Polynomia Interpolation method to get the
					// wave functions

				tmp = x12 * (psi_uniform[ip] * x3 + psi_uniform[ip + 3] * x0)
						+ x03 * (psi_uniform[ip + 1] * x2 - psi_uniform[ip + 2] * x1);

				dtmp = x12 * (dpsi_uniform[ip] * x3 + dpsi_uniform[ip + 3] * x0)
						+ x03 * (dpsi_uniform[ip + 1] * x2 - dpsi_uniform[ip + 2] * x1);
			} // new l is used.

			// get the 'l' of this localized wave function
			const int ll = atom->iw2l[iw];
			const int idx_lm = atom->iw2_ylm[iw];

			const double rl = pow_int(distance1, ll);

			// derivative of wave functions with respect to atom positions.
			const double tmpdphi_rly = (dtmp - tmp * ll / distance1) / rl * rly[idx_lm] / distance1;
			const double tmprl = tmp / rl;

			dpsi[iw][i][0] = tmpdphi_rly * dr1[0] + tmprl * grly[idx_lm][0];
			dpsi[iw][i][1] = tmpdphi_rly * dr1[1] + tmprl * grly[idx_lm][1];
			dpsi[iw][i][2] = tmpdphi_rly * dr1[2] + tmprl * grly[idx_lm][2];
		} // end iw
}     // end i = 0-6



	// atomic basis sets
	// psir_vlbr3[bxyz][LD_pool]
    ModuleBase::Array_Pool<double> get_psir_vlbr3(
        const int bxyz,
        const int na_grid,  					    // how many atoms on this (i,j,k) grid
		const int LD_pool,
		const int*const block_index,		    	// block_index[na_grid+1], count total number of atomis orbitals
		const bool*const*const cal_flag,	    	// cal_flag[bxyz][na_grid],	whether the atom-grid distance is larger than cutoff
		const double*const vldr3,			    	// vldr3[bxyz]
		const double*const*const psir_ylm)		    // psir_ylm[bxyz][LD_pool]
	{
		ModuleBase::Array_Pool<double> psir_vlbr3(bxyz, LD_pool);
		for(int ib=0; ib<bxyz; ++ib)
		{
			for(int ia=0; ia<na_grid; ++ia)
			{
				if(cal_flag[ib][ia])
				{
					for(int i=block_index[ia]; i<block_index[ia+1]; ++i)
					{
						psir_vlbr3[ib][i]=psir_ylm[ib][i]*vldr3[ib];
					}
				}
				else
				{
					for(int i=block_index[ia]; i<block_index[ia+1]; ++i)
					{
						psir_vlbr3[ib][i]=0;
					}
				}

			}
		}
		return psir_vlbr3;
	}

} // namespace Gint_Tools
