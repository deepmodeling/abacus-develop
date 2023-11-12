#ifndef ATEN_KERNELS_LAPACK_H_
#define ATEN_KERNELS_LAPACK_H_

#include <ATen/core/tensor.h>
#include <ATen/core/tensor_types.h>

#include <base/third_party/lapack.h>

namespace container {
namespace kernels {


template <typename T, typename Device>
struct set_matrix {
    void operator() (
        const char& uplo,
        T* A,
        const int& dim);
};


template <typename T, typename Device>
struct lapack_trtri {
    void operator()(
        const char& uplo,
        const char& diag,
        const int& dim,
        T* Mat,
        const int& lda);
};


template <typename T, typename Device>
struct lapack_potrf {
    void operator()(
        const char& uplo,
        const int& dim,
        T* Mat, 
        const int& lda);
};


template <typename T, typename Device>
struct lapack_dnevd {
    using Real = typename GetTypeReal<T>::type;
    void operator()(
        const char& jobz,
        const char& uplo,
        T* Mat,
        const int& dim,
        Real* eigen_val);
};


template <typename T, typename Device>
struct lapack_dngvd {
    using Real = typename GetTypeReal<T>::type;
    void operator()(
        const int& itype,
        const char& jobz,
        const char& uplo,
        T* Mat_A,
        T* Mat_B,
        const int& dim,
        Real* eigen_val);
};

#if defined(__CUDA) || defined(__ROCM)
// TODO: Use C++ singleton to manage the GPU handles
void createGpuSolverHandle();  // create cusolver handle
void destroyGpuSolverHandle(); // destroy cusolver handle
#endif

} // namespace container
} // namespace kernels

#endif // ATEN_KERNELS_LAPACK_H_