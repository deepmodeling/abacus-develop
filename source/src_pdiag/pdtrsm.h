#ifndef PDTRSM_H
#define PDTRSM_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "pdiag_common.h"

void pdtrsm(char isuplo,int b_n,MPI_Comm comm_2D,int NB,int N_A,
	    double *A,double *B,LocalMatrix loc_A,LocalMatrix loc_B);

#endif
