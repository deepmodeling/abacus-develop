//==========================================================
// AUTHOR : Fang Wei, Mohan Chen
// DATE : 2008
// LAST UPDATE : 2009-3-23 mohan add GATHER_MINIMUM_DOUBLE
//==========================================================
#ifndef PARALLEL_GLOBAL_H
#define PARALLEL_GLOBAL_H
#include "../src_spillage/common.h"

#ifdef __MPI
#include <mpi.h>
extern MPI_Datatype mpicomplex;
extern MPI_Op myOp;
extern MPI_Comm POOL_WORLD;
#endif

//void myProd(complex<double> *in,complex<double> *inout,int *len,MPI_Datatype *dptr);

namespace Parallel_Global
{
	void read_mpi_parameters(int argc, char **argv);

#ifdef __MPI
	void myProd(complex<double> *in,complex<double> *inout,int *len,MPI_Datatype *dptr);
#endif
}


#endif // GMPI
