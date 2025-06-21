#pragma once
#include <complex>
#include <mpi.h>

static inline int globalIndex(int localIndex, int nblk, int nprocs, int myproc)
{
    int iblock, gIndex;
    iblock = localIndex / nblk;
    gIndex = (iblock * nprocs + myproc) * nblk + localIndex % nblk;
    return gIndex;
}

static inline int localIndex(int globalIndex, int nblk, int nprocs, int& lcoalProc)
{
    lcoalProc = int((globalIndex % (nblk * nprocs)) / nblk);
    return int(globalIndex / (nblk * nprocs)) * nblk + globalIndex % nblk;
}

#ifdef __MPI
void initBlacsGrid(int loglevel,
                   MPI_Comm comm,
                   int nFull,
                   int nblk,
                   int& blacs_ctxt,
                   int& narows,
                   int& nacols,
                   int desc[]);
#endif

// load matrix from the file
void loadMatrix(const char FileName[], int nFull, double* a, int* desca, int blacs_ctxt);

void saveLocalMatrix(const char filePrefix[], int narows, int nacols, double* a);

// use pdgemr2d to collect matrix from all processes to root process
// and save to one completed matrix file
void saveMatrix(const char FileName[], int nFull, double* a, int* desca, int blacs_ctxt);

void loadMatrix(const char FileName[], int nFull, std::complex<double>* a, int* desca, int blacs_ctxt);

void saveLocalMatrix(const char filePrefix[], int narows, int nacols, std::complex<double>* a);

// use pzgemr2d to collect matrix from all processes to root process
// and save to one completed matrix file
void saveMatrix(const char FileName[], int nFull, std::complex<double>* a, int* desca, int blacs_ctxt);
