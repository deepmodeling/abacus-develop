#include "gint_gamma.h"
#include "src_pw/global.h"

// here vindex refers to local potentials
int* get_vindex(
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
	int* vindex = get_vindex(ncyz, ibx, jby, kbz);	
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

int* Gint_Gamma::get_colidx(
	const int na_grid,  		// how many atoms on this (i,j,k) grid
	const int grid_index)		// 1d index of FFT index (i,j,k)
{
	int* colidx = (int*)malloc((na_grid+1)*sizeof(int));
    colidx[0] = 0;
    for (int id=0; id<na_grid; id++)
	{
        const int mcell_index = GridT.bcell_start[grid_index] + id;
        const int iat = GridT.which_atom[mcell_index]; // index of atom
        const int it = ucell.iat2it[iat]; // index of atom type
        colidx[id+1] = colidx[id]+ucell.atoms[it].nw;
	}
	return colidx;
}

// band size: number of columns of a band
int* Gint_Gamma::get_bsize(
	const int na_grid,			// how many atoms on this (i,j,k) grid
	const int grid_index)		// 1d index of FFT index (i,j,k)
{
	int* bsize = (int*)malloc(na_grid*sizeof(int));
    for (int id=0; id<na_grid; id++)
	{
        const int mcell_index=GridT.bcell_start[grid_index] + id;
        const int iat=GridT.which_atom[mcell_index]; // index of atom
        const int it=ucell.iat2it[ iat ]; // index of atom type
        bsize[id]=ucell.atoms[it].nw;	
	}
	return bsize;
}
