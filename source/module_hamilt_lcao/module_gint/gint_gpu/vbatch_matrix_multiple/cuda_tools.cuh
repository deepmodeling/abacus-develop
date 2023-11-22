#ifndef CUDA_TOOLS_CUH
#define CUDA_TOOLS_CUH
#include <cuda.h>    // for CUDA_VERSION
#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <fstream>
#include <sstream>
#include <iostream>

typedef cuDoubleComplex magmaDoubleComplex;
typedef cuFloatComplex  magmaFloatComplex;
#define C_MAKE(r,i)     make_cuFloatComplex(r, i)
#define Z_MAKE(r,i)     make_cuDoubleComplex(r, i)    ///< @return complex number r + i*sqrt(-1).

__device__ static __inline__ float 
satomic_add(float* address, float val)
{
    return atomicAdd(address, val);
}

/******************************************************************************/
__device__ static __inline__ double 
datomic_add(double* address, double val)
{
    return atomicAdd(address, val);
}
/******************************************************************************/
__device__ static __inline__ magmaFloatComplex 
catomic_add(magmaFloatComplex* address, magmaFloatComplex val)
{
    float re = satomic_add( (float*) (&(*address).x) ,val.x);
    float im = satomic_add( (float*) (&(*address).y) ,val.y);
    return C_MAKE(re, im);
}

/******************************************************************************/
__device__ static __inline__ magmaDoubleComplex 
zatomic_add(magmaDoubleComplex* address, magmaDoubleComplex val)
{
    double re = datomic_add( (double*) (&(*address).x) ,val.x);
    double im = datomic_add( (double*) (&(*address).y) ,val.y);
    return Z_MAKE(re, im);
}
#define fetch(A, m, n, bound)  offs_d##A[min(n*LD##A+m, bound)]

// =============================================================================
#if defined(PRECISION_z)
    #define add(A, B)        Z_ADD(A, B)
    #define mul(A, B)        Z_MUL(A, B)
    #define div(A, B)        Z_DIV(A, B)
    #define fma(A, B, C) C = cuCfma(A, B, C)
    #define make_FloatingPoint(x, y) Z_MAKE(x, y)
#elif defined(PRECISION_c)
    #define add(A, B)        C_ADD(A, B)
    #define mul(A, B)        C_MUL(A, B)
    #define div(A, B)        C_DIV(A, B)
    #define fma(A, B, C) C = cuCfmaf(A, B, C)
    #define make_FloatingPoint(x, y) C_MAKE(x, y)
#elif defined(PRECISION_h)
    #define add(A, B)         (A+B)
    #define mul(A, B)         (A*B)
    #define div(A, B)         (A/B)
    #define fma(A, B, C) C += (A*B)
    #define make_FloatingPoint(x, y) ((magmaHalf)x)
#else
    #define add(A, B)         (A+B)
    #define mul(A, B)         (A*B)
    #define div(A, B)         (A/B)
    #define fma(A, B, C) C += (A*B)
    #define make_FloatingPoint(x, y) (x)
#endif


#if defined(PRECISION_z)
    #define atomic_add zatomic_add
#elif defined(PRECISION_c)
    #define atomic_add catomic_add
#elif defined(PRECISION_d)
    #define atomic_add datomic_add
#else
    #define atomic_add satomic_add
#endif

void dump_cuda_array_to_file(double *cuda_array, int width, int hight, const std::string &filename);


#endif // CUDA_TOOLS_CUH#ifndef CUDA_TOOLS_CUH
