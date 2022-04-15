#ifndef PDSYTRD_H
#define PDSYTRD_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "pdiag_common.h"

void pdsytrd(MPI_Comm comm2D, LocalMatrix loc_A,int N,int NB,
	    double *a, double *diag, double *off_diag, double *norm,char uplo);

#endif
