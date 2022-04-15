#ifndef PZHETRD_H
#define PZHETRD_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "pdiag_common.h"
#include "pzhetrd.h"

void pzhetrd(MPI_Comm comm2D, LocalMatrix loc_A,int N,int NB,
	    std::complex<double> *a,double *diag,double *off_diag,std::complex<double>  *norm,char uplo);

#endif
