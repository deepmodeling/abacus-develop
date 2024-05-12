#ifndef DIAG_CUSOLVER_CUH
#define DIAG_CUSOLVER_CUH

#include <complex>
#include "nvToolsExt.h"
#include <cuda_runtime.h>
#include <cusolverDn.h>


///
/// The diag_cusolver contains subroutines for generalized eigenvector problems of dense matrice in lcao on a single GPU.
/// Cusolver apis, i.e. cusolverDnDsygvd and cusolverDnZhegvd, are wrapped with respect to different data types for easy calling.
/// For details of low-level methods relevent, one can refer to [cuSOLVER documentation](https://docs.nvidia.com/cuda/cusolver/index.html#cuSolverDN-lt-t-gt-sygvd).
///
// Xu Shu add 2022-03-23
//

class Diag_Cusolver_gvd{

//-------------------
// private variables
//-------------------

    cusolverDnHandle_t cusolverH;

    cusolverEigType_t itype; //problem type: A*x = (lambda)*B*x
    cusolverEigMode_t jobz; // compute eigenvalues and eigenvectors.
    cublasFillMode_t uplo;

    int m;
    int lda;

    double *d_A;
    double *d_B;
    double *d_work;
    
    cuDoubleComplex *d_A2;
    cuDoubleComplex *d_B2;
    cuDoubleComplex *d_work2;

    double *d_W ;
    int *devInfo;

    int  lwork;
    int info_gpu;

//   subroutines that are related to initializing the class:
//  - init_double : initializing relevant double type data structures and gpu apis' handle and memory
//  - init_complex : initializing relevant complex type data structures and gpu apis' handle and memory
//      Input Parameters
//          N: the dimension of the matrix 
    void init_double(int N);
    void init_complex(int N);

    void finalize();  // for recycling the usage of the static class Diag_Cusolver_gvd
public:

    int is_init;    // For expensive gpu initialization only once when using cusolver for lcao

    Diag_Cusolver_gvd();
    ~Diag_Cusolver_gvd();

//   subroutines that are related to calculating generalized eigenvalues and eigenvectors for dense matrix pairs:
//  - Dngvd_double : dense double type matrix
//  - Dngvd_complex : dense complex type matrix
//      Input Parameters
//          N: the number of rows of the matrix 
//          M: the number of cols of the matrix  
//          A: the hermitian matrix A in A x=lambda B (column major) 
//          B: the SPD matrix B in A x=lambda B (column major) 
//      Output Parameter
//          W: generalized eigenvalues
//          V: generalized eigenvectors (column major)

    void Dngvd_double(int N, int M, double *A, double *B, double *W, double *V);
    void Dngvd_complex(int N, int M, std::complex<double> *A, std::complex<double> *B, double *W, std::complex<double> *V);
    
    void Dngvd(int N, int M, double *A, double *B, double *W, double *V)
    {
        return Dngvd_double(N, M, A, B, W, V);
    };

    void Dngvd(int N, int M, std::complex<double> *A, std::complex<double> *B, double *W, std::complex<double> *V)
    {
        return Dngvd_complex(N, M, A, B, W, V);
    };

};

#endif
