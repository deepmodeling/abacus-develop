#pragma once
#include "mpi.h"

#include <complex>
#include <fstream>
#include <vector>
#include <cal.h>
#include <cusolverMp.h>
#include "module_base/macros.h"

template<typename inputT>
class Diag_CusolverMP_gvd
{
  private:
    using outputT = typename GetTypeReal<inputT>::type;

  public:
    Diag_CusolverMP_gvd(const MPI_Comm comm,
                          const int narows,
                          const int nacols,
                          const int* desc);


    // 20240530 zhanghaochong
    // Here eigen is the output rather than the input. Don't be confused by inputT.
    // Here inputT is the data type, the data type of EigenVector should be the same as
    // the data type of the input matrix.
    int generalized_eigenvector(inputT* A,
                                inputT* B,
                                outputT* EigenValue,
                                inputT* EigenVector);
    ~Diag_CusolverMP_gvd();
    void outputParameters();
  private:

    int nFull;
    int cblacs_ctxt;
    int lda;

    // num of processor rows and cols
    int nprows;
    int npcols;
    // my processor row and col
    int myprow;
    int mypcol;

    int m_local;
    int n_local;
    int comm_f;
    // my global mpi id
    int globalMpiRank;
    // global mpi size
    int globalMpiSize;
    cudaDataType_t datatype;

    cal_comm_t cusolverCalComm = NULL;
    cudaStream_t localStream = NULL;
    cusolverMpHandle_t cusolverMpHandle = NULL;
    cusolverMpGrid_t grid = NULL;

    /* cusolverMp matrix descriptors */
    cusolverMpMatrixDescriptor_t desc_for_cusolvermp = NULL;

    int64_t matrix_i;
    int64_t matrix_j;
};

// 实现模板类的成员函数

