#include "blas_connector.h"
#include "../macros.h"

#ifdef __CUDA
#include <base/macros/macros.h>
#include <cuda_runtime.h>
#include "cublas_v2.h"
#include "source_base/kernels/math_kernel_op.h"
#include "source_base/module_device/memory_op.h"


namespace BlasUtils{

	void createGpuBlasHandle(){
		if (cublas_handle == nullptr) {
			cublasErrcheck(cublasCreate(&cublas_handle));
		}
	}

	void destoryBLAShandle(){
		if (cublas_handle != nullptr) {
			cublasErrcheck(cublasDestroy(cublas_handle));
			cublas_handle = nullptr;
		}
	}


	cublasOperation_t judge_trans(bool is_complex, const char& trans, const char* name)
	{
		if (trans == 'N')
		{
			return CUBLAS_OP_N;
		}
		else if(trans == 'T')
		{
			return CUBLAS_OP_T;
		}
		else if(is_complex && trans == 'C')
		{
			return CUBLAS_OP_C;
		}
		return CUBLAS_OP_N;
	}

	cublasSideMode_t judge_side(const char& trans)
	{
		if (trans == 'L')
		{
			return CUBLAS_SIDE_LEFT;
		}
		else if (trans == 'R')
		{
			return CUBLAS_SIDE_RIGHT;
		}
		return CUBLAS_SIDE_LEFT;
	}

	cublasFillMode_t judge_fill(const char& trans)
	{
		if (trans == 'F')
		{
			return CUBLAS_FILL_MODE_FULL;
		}
		else if (trans == 'U')
		{
			return CUBLAS_FILL_MODE_UPPER;
		}
		else if (trans == 'D')
		{
			return CUBLAS_FILL_MODE_LOWER;
		}
		return CUBLAS_FILL_MODE_FULL;
	}

} // namespace BlasUtils

#endif