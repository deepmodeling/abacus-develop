#ifndef PDSTEBZ_H
#define PDSTEBZ_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "pdiag_common.h"

void pdstebz(MPI_Comm comm_2D, double *D,double *E,double *eigen,int N);

#endif
