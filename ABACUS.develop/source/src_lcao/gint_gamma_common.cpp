#include "gint_gamma.h"
#include "src_pw/global.h"
#include "src_global/ylm.h"

// here vindex refers to local potentials
int* Gint_Gamma::get_vindex(
	const int ncyz,
	const int ibx,
	const int jby,
	const int kbz)
{
	int *vindex = (int*)malloc(pw.bxyz*sizeof(int));
    int bindex=0;
    // z is the fastest,
	// ipart can be obtained by using a previously stored array
    for(int ii=0; ii<pw.bx; ii++)
    {
        const int ipart = (ibx+ii)*ncyz;
        for(int jj=0; jj<pw.by; jj++)
        {
			// jpart can be obtained by using a previously stored array
            const int jpart = (jby+jj)*pw.nczp + ipart;
            for(int kk=0; kk<pw.bz; kk++)
            {
                vindex[bindex] = kbz+kk + jpart;
				++bindex;
            }
        }
    }
	return vindex;
}

// extract the local potentials.
double* Gint_Gamma::get_vldr3(
	const int ncyz,
	const int ibx,
	const int jby,
	const int kbz) const
{
	// set the index for obtaining local potentials
	int* vindex = Gint_Gamma::get_vindex(ncyz, ibx, jby, kbz);	
	double *vldr3 = (double*)malloc(pw.bxyz*sizeof(double));					
	for(int ib=0; ib<pw.bxyz; ib++)
	{
		vldr3[ib]=this->vlocal[vindex[ib]] * this->vfactor;
	}
	free(vindex);	vindex=nullptr;
	return vldr3;
}

// index of wave functions for each block
int* Gint_Gamma::get_block_iw(
	const int na_grid,  		// how many atoms on this (i,j,k) grid
	const int grid_index,		// 1d index of FFT index (i,j,k))
	const int max_size)
{
	int *block_iw = (int*)malloc(max_size*sizeof(int));
    for (int id=0; id<na_grid; id++)
    {
        const int mcell_index=GridT.bcell_start[grid_index] + id;
        const int iat=GridT.which_atom[mcell_index]; // index of atom
        const int it=ucell.iat2it[ iat ]; // index of atom type
        const int ia=ucell.iat2ia[ iat ]; // index of atoms within each type
        const int start=ucell.itiaiw2iwt(it, ia, 0); // the index of the first wave function for atom (it,ia)
        block_iw[id]=GridT.trace_lo[start];
	}
	return block_iw;
}

int* Gint_Gamma::get_block_index(
	const int na_grid,  		// how many atoms on this (i,j,k) grid
	const int grid_index)		// 1d index of FFT index (i,j,k)
{
	int* block_index = (int*)malloc((na_grid+1)*sizeof(int));
    block_index[0] = 0;
    for (int id=0; id<na_grid; id++)
	{
        const int mcell_index = GridT.bcell_start[grid_index] + id;
        const int iat = GridT.which_atom[mcell_index]; // index of atom
        const int it = ucell.iat2it[iat]; // index of atom type
        block_index[id+1] = block_index[id]+ucell.atoms[it].nw;
	}
	return block_index;
}

// band size: number of columns of a band
int* Gint_Gamma::get_block_size(
	const int na_grid,			// how many atoms on this (i,j,k) grid
	const int grid_index)		// 1d index of FFT index (i,j,k)
{
	int* block_size = (int*)malloc(na_grid*sizeof(int));
    for (int id=0; id<na_grid; id++)
	{
        const int mcell_index=GridT.bcell_start[grid_index] + id;
        const int iat=GridT.which_atom[mcell_index]; // index of atom
        const int it=ucell.iat2it[ iat ]; // index of atom type
        block_size[id]=ucell.atoms[it].nw;	
	}
	return block_size;
}

// whether the atom-grid distance is larger than cutoff
bool** Gint_Gamma::get_cal_flag(
	const int na_grid, 		// number of atoms on this grid 
	const int grid_index)
{
	bool** cal_flag = (bool**)malloc(pw.bxyz*sizeof(bool*));
	for(int ib=0; ib<pw.bxyz; ++ib)
		cal_flag[ib] = (bool*)malloc(na_grid*sizeof(bool));

	for (int id=0; id<na_grid; id++)
	{
		// there are two parameters we want to know here:
		// in which bigcell of the meshball the atom in?
		// what's the cartesian coordinate of the bigcell?
		const int mcell_index=GridT.bcell_start[grid_index] + id;
		const int iat=GridT.which_atom[mcell_index];		
		const int it=ucell.iat2it[iat];

		// meshball_positions should be the bigcell position in meshball
		// to the center of meshball.
		// calculated in cartesian coordinates
		// the vector from the grid which is now being operated to the atom position.
		// in meshball language, is the vector from imcell to the center cel, plus
		// tau_in_bigcell.
		const int imcell=GridT.which_bigcell[mcell_index];
		const double mt[3] = {
			GridT.meshball_positions[imcell][0] - GridT.tau_in_bigcell[iat][0],
			GridT.meshball_positions[imcell][1] - GridT.tau_in_bigcell[iat][1],
			GridT.meshball_positions[imcell][2] - GridT.tau_in_bigcell[iat][2]};

		for(int ib=0; ib<pw.bxyz; ib++)
		{
			// meshcell_pos: z is the fastest
			const double dr[3] = {
				GridT.meshcell_pos[ib][0] + mt[0],
				GridT.meshcell_pos[ib][1] + mt[1],
				GridT.meshcell_pos[ib][2] + mt[2]};
			const double distance = std::sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);	// distance between atom and grid

			if(distance > (ORB.Phi[it].getRcut()-1.0e-15))
				cal_flag[ib][id]=false;
			else
				cal_flag[ib][id]=true;
		}// end ib
	}// end id
	return cal_flag;
}


Gint_Gamma::Array_Pool<double> Gint_Gamma::cal_psir_ylm(
	const int na_grid, // number of atoms on this grid 
	const int LD_pool,
	const int grid_index, // 1d index of FFT index (i,j,k) 
	const double delta_r, // delta_r of the uniform FFT grid
	const int*const block_index,  // count total number of atomis orbitals
	const int*const block_size, 
	const bool*const*const cal_flag) // whether the atom-grid distance is larger than cutoff
{
	Array_Pool<double> psir_ylm(pw.bxyz, LD_pool);
	for (int id=0; id<na_grid; id++)
	{
		// there are two parameters we want to know here:
		// in which bigcell of the meshball the atom is in?
		// what's the cartesian coordinate of the bigcell?
		const int mcell_index=GridT.bcell_start[grid_index] + id;

		const int iat=GridT.which_atom[mcell_index]; // index of atom
		const int it=ucell.iat2it[iat]; // index of atom type
		const Atom*const atom=&ucell.atoms[it];

		// meshball_positions should be the bigcell position in meshball
		// to the center of meshball.
		// calculated in cartesian coordinates
		// the vector from the grid which is now being operated to the atom position.
		// in meshball language, is the vector from imcell to the center cel, plus
		// tau_in_bigcell.
		const int imcell=GridT.which_bigcell[mcell_index];
		const double mt[3] = {
			GridT.meshball_positions[imcell][0] - GridT.tau_in_bigcell[iat][0],
			GridT.meshball_positions[imcell][1] - GridT.tau_in_bigcell[iat][1],
			GridT.meshball_positions[imcell][2] - GridT.tau_in_bigcell[iat][2]};

		// number of grids in each big cell (bxyz)
		for(int ib=0; ib<pw.bxyz; ib++)
		{
			double *p=&psir_ylm.ptr_2D[ib][block_index[id]];
			if(!cal_flag[ib][id]) 
			{
				ZEROS(p, block_size[id]);
			}
			else
			{
				// meshcell_pos: z is the fastest
				const double dr[3] = {
					GridT.meshcell_pos[ib][0] + mt[0],
					GridT.meshcell_pos[ib][1] + mt[1],
					GridT.meshcell_pos[ib][2] + mt[2]};	
				double distance = std::sqrt( dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2] );	// distance between atom and grid
				//if(distance[id] > GridT.orbital_rmax) continue;
				if (distance < 1.0E-9) distance += 1.0E-9;
				
				//------------------------------------------------------
				// spherical harmonic functions Ylm
				//------------------------------------------------------
				std::vector<double> ylma;
				//	Ylm::get_ylm_real(this->nnn[it], this->dr[id], ylma);
				Ylm::sph_harm ( ucell.atoms[it].nwl,
						dr[0] / distance,
						dr[1] / distance,
						dr[2] / distance,
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

				const double c3 = 3.0*dx2-2.0*dx3;
				const double c1 = 1.0-c3;
				const double c2 = (dx-2.0*dx2+dx3)*delta_r;
				const double c4 = (dx3-dx2)*delta_r;

				double phi=0;
				for (int iw=0; iw< atom->nw; ++iw, ++p)
				{
					if ( atom->iw2_new[iw] )
					{
						const Numerical_Orbital_Lm &philn = ORB.Phi[it].PhiLN(
								atom->iw2l[iw],
								atom->iw2n[iw]);
						phi = c1*philn.psi_uniform[ip] + c2*philn.dpsi_uniform[ip]			 // radial wave functions
							+ c3*philn.psi_uniform[ip+1] + c4*philn.dpsi_uniform[ip+1];
					}
					*p=phi * ylma[atom->iw2_ylm[iw]];
				} // end iw
			}// end distance<=(ORB.Phi[it].getRcut()-1.0e-15)
		}// end ib
	}// end id
	return psir_ylm;
}
