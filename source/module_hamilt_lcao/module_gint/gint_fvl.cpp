#include "gint_k.h"
#include "module_base/timer.h"
#include "module_base/ylm.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_base/array_pool.h"

// This function utilizes the cache more effectively than calling the ddot function, thus performing faster.
void Gint::cal_meshball_force(
    const int grid_index,
    const int na_grid,  					    // how many atoms on this (i,j,k) grid
	const int*const block_size, 			    // block_size[na_grid],	number of columns of a band
	const int*const block_index,		    	// block_index[na_grid+1], count total number of atomis orbitals
	const double*const*const psir_vlbr3_DMR,	    // psir_vlbr3[this->bxyz][LD_pool]
    const double*const*const dpsir_x,	    // psir_vlbr3[this->bxyz][LD_pool]
    const double*const*const dpsir_y,	    // psir_vlbr3[this->bxyz][LD_pool]
    const double*const*const dpsir_z,	    // psir_vlbr3[this->bxyz][LD_pool]
    ModuleBase::matrix *force)
{
    for(int ia1=0;ia1<na_grid;ia1++)
    {
        const int mcell_index=this->gridt->bcell_start[grid_index] + ia1;
        const int iat=this->gridt->which_atom[mcell_index]; // index of atom
		double rx = 0;
		double ry = 0;
		double rz = 0;
        for(int ib=0; ib<this->bxyz; ib++)
        {
            for(int iw=0; iw<block_size[ia1]; iw++)
			{
				double psir_vlbr3 = psir_vlbr3_DMR[ib][block_index[ia1]+iw];
				rx += psir_vlbr3 * dpsir_x[ib][block_index[ia1]+iw];
				ry += psir_vlbr3 * dpsir_y[ib][block_index[ia1]+iw];
				rz += psir_vlbr3 * dpsir_z[ib][block_index[ia1]+iw];
			}
        }
        force[0](iat,0) += rx * 2.0;
        force[0](iat,1) += ry * 2.0;
		force[0](iat,2) += rz * 2.0;  
    }
	return;
}

// This function utilizes the cache more effectively than calling the ddot function, thus performing faster.
void Gint::cal_meshball_stress(
    const int na_grid,  					    // how many atoms on this (i,j,k) grid
	const int*const block_index,		    	// block_index[na_grid+1], count total number of atomis orbitals
	const double*const psir_vlbr3_DMR,
    const double*const dpsirr,
    ModuleBase::matrix *stress)
{
	double rxx = 0;
	double rxy = 0;
	double rxz = 0;
	double ryy = 0;
	double ryz = 0;
	double rzz = 0;
	const int size = block_index[na_grid] * this->bxyz;

    for(int i=0; i<size; ++i)
    {
		double psir_vlbr3 = psir_vlbr3_DMR[i];
		rxx += psir_vlbr3 * dpsirr[i * 6];
		rxy += psir_vlbr3 * dpsirr[i * 6 + 1];
		rxz += psir_vlbr3 * dpsirr[i * 6 + 2];
		ryy += psir_vlbr3 * dpsirr[i * 6 + 3];
		ryz += psir_vlbr3 * dpsirr[i * 6 + 4];
		rzz += psir_vlbr3 * dpsirr[i * 6 + 5];
    }
	stress[0](0,0) += rxx*2;
    stress[0](0,1) += rxy*2;
	stress[0](0,2) += rxz*2;
    stress[0](1,1) += ryy*2;
    stress[0](1,2) += ryz*2;
    stress[0](2,2) += rzz*2;
    return;
}
