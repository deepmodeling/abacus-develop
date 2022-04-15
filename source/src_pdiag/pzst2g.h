#ifndef PZST2G_H
#define PZST2G_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "pdiag_common.h"

void pzst2g(MPI_Comm comm_2D,int NB,int N_A,std::complex<double> *A,std::complex<double> *B,
	     LocalMatrix loc_A,int loc_size,char isuplo);

#endif
