#ifndef BASE_THIRD_PARTY_BLAS_H_
#define BASE_THIRD_PARTY_BLAS_H_

#include <complex>
#include "source_base/module_external/blas_connector.h"

#if defined(__CUDA)
#include <base/third_party/cublas.h>
#elif defined(__ROCM)
#include <base/third_party/hipblas.h>
#endif


#endif // BASE_THIRD_PARTY_BLAS_H_
