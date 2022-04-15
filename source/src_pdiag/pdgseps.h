#ifndef PDGSEPS_H
#define PDGSEPS_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "pdiag_common.h"

void pdgseps(
		MPI_Comm comm_2D,
		int n,
		int nb,
		double *A, 
		double *B,
		double *Z,
		double *eigen,
		LocalMatrix loc_A,
		char uplo, 
		int &loc_size,
		int &loc_pos);

#endif
