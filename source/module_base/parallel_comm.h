#ifndef PARALLEL_COMM_H
#define PARALLEL_COMM_H

#ifdef __MPI
#include "mpi.h"

MPI_Comm POOL_WORLD;
MPI_Comm INTER_POOL = MPI_COMM_NULL; // communicator among different pools
MPI_Comm STO_WORLD;
MPI_Comm PARAPW_WORLD;
MPI_Comm GRID_WORLD; // mohan add 2012-01-13
MPI_Comm DIAG_WORLD; // mohan add 2012-01-13

#endif

#endif // PARALLEL_COMM_H