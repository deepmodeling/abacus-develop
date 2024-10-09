#ifndef PARALLEL_COMM_H
#define PARALLEL_COMM_H

#ifdef __MPI

#include "mpi.h"
extern MPI_Comm POOL_WORLD;
extern MPI_Comm INTER_POOL; // communicator among different pools
extern MPI_Comm STO_WORLD;
extern MPI_Comm PARAPW_WORLD;
extern MPI_Comm GRID_WORLD; // mohan add 2012-01-13
extern MPI_Comm DIAG_WORLD; // mohan add 2012-01-13

#endif

#endif // PARALLEL_COMM_H