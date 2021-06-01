#include "gint_gamma.h"
#include "gint_tools.h"
#include "src_pw/global.h"

// extract the local potentials.
double* Gint_Gamma::get_vldr3(
	const int ncyz,
	const int ibx,
	const int jby,
	const int kbz) const
{
	// set the index for obtaining local potentials
	int* vindex = Gint_Tools::get_vindex(ncyz, ibx, jby, kbz);	
	double *vldr3 = (double*)malloc(pw.bxyz*sizeof(double));					
	for(int ib=0; ib<pw.bxyz; ib++)
	{
		vldr3[ib]=this->vlocal[vindex[ib]] * this->vfactor;
	}
	free(vindex);	vindex=nullptr;
	return vldr3;
}
