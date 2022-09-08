#pragma once
// blacs
    // Initialization
#include "mpi.h"
int Csys2blacs_handle(MPI_Comm SysCtxt);
void Cblacs_pinfo(int *myid, int *nprocs);
void Cblacs_get(int icontxt, int what, int *val);
void Cblacs_gridinit(int* icontxt, char *layout, int nprow, int npcol);
void Cblacs_gridmap(int* icontxt, int *usermap, int ldumap, int nprow, int npcol);
    // Destruction
void Cblacs_gridexit(int icontxt);
    // Informational and Miscellaneous
void Cblacs_gridinfo(int icontxt, int* nprow, int *npcol, int *myprow, int *mypcol);
int Cblacs_pnum(int icontxt, int prow, int pcol);
void Cblacs_pcoord(int icontxt, int pnum, int *prow, int *pcol);
void Cblacs_barrier(int icontxt, char *scope);
    // Point to  Point
void Cdgesd2d(int icontxt, int m, int n, double *a, int lda, int rdest, int cdest);
void Cdgerv2d(int icontxt, int m, int n, double *a, int lda, int rsrc, int csrc);
void Czgesd2d(int icontxt, int m, int n, double _Complex *a, int lda, int rdest, int cdest);
void Czgerv2d(int icontxt, int m, int n, double _Complex *a, int lda, int rsrc, int csrc);
    // Combine
//void Cdgamx2d(int icontxt, int scope, int top, int m, int n,
//              double *a, int lda, int *ra, int *ca, int rcflag, int rdest, int cdest);
