#ifndef PDIAG_COMMON_H
#define PDIAG_COMMON_H

#ifdef __MPI
#include "mpi.h"
#endif
#include "module_base/blas_connector.h"
#include "module_base/lapack_connector.h"			// Peize Lin add 2016-08-04

struct LocalMatrix
{
    int *row_set;
    int *col_set;
    
	int col_num;
    int row_num;
    
	int col_pos;
    int row_pos;
    
	int row_b;  //row block size
    int col_b;  //column block size
};

#ifdef __MPI
void indxg2l(int ia,int ja,int nb,int nprow,int npcol,int *iia,int *jja);
void indxg2p(MPI_Comm comm2D,int nb,int ia,int ja,int *iarow,int *iacol);
void mpi_sub_col(MPI_Comm comm2D, MPI_Comm *comm_col);
void mpi_sub_row(MPI_Comm comm2D, MPI_Comm *comm_row);
extern MPI_Comm DIAG_HPSEPS_WORLD; //mohan add 2012-01-13
#endif

#endif
