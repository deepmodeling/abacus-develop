#include "blas_connector.h"
#include "../macros.h"

#ifdef __CUDA
#include <base/macros/macros.h>
#include <cuda_runtime.h>
#include "cublas_v2.h"
#include "source_base/kernels/math_kernel_op.h"
#include "source_base/module_device/memory_op.h"


namespace BlasUtils{

    cublasHandle_t cublas_handle = nullptr;

    void createGpuBlasHandle(){
        if (cublas_handle == nullptr) {
            CHECK_CUBLAS(cublasCreate(&cublas_handle));
        }
    }

    void destoryBLAShandle(){
        if (cublas_handle != nullptr) {
            CHECK_CUBLAS(cublasDestroy(cublas_handle));
            cublas_handle = nullptr;
        }
    }


    cublasOperation_t judge_trans(bool is_complex, const char& trans, const char*)
    {
        return GetCublasOperation(trans, is_complex);
    }

    cublasSideMode_t judge_side(const char& side)
    {
        if (side == 'L' || side == 'l')
        {
            return CUBLAS_SIDE_LEFT;
        }
        else if (side == 'R' || side == 'r')
        {
            return CUBLAS_SIDE_RIGHT;
        }
        throw std::invalid_argument(std::string("Invalid BLAS side: ") + side);
    }

    cublasFillMode_t judge_fill(const char& uplo)
    {
        return cublas_fill_mode(uplo);
    }

} // namespace BlasUtils

#endif
