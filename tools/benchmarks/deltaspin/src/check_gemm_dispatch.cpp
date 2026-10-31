#include "module_hsolver/kernels/math_kernel_op.h"

#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <complex>
#include <cstdlib>
#include <iostream>

namespace
{
int three_real_calls = 0;

void check(cudaError_t status)
{
    if (status != cudaSuccess) { std::cerr << cudaGetErrorString(status) << '\n'; std::exit(1); }
}
} // namespace

extern "C" cublasStatus_t __real_cublasZgemm3m(
    cublasHandle_t, cublasOperation_t, cublasOperation_t, int, int, int,
    const cuDoubleComplex*, const cuDoubleComplex*, int, const cuDoubleComplex*, int,
    const cuDoubleComplex*, cuDoubleComplex*, int);

extern "C" cublasStatus_t __wrap_cublasZgemm3m(
    cublasHandle_t handle, cublasOperation_t ta, cublasOperation_t tb, int m, int n, int k,
    const cuDoubleComplex* alpha, const cuDoubleComplex* a, int lda, const cuDoubleComplex* b, int ldb,
    const cuDoubleComplex* beta, cuDoubleComplex* c, int ldc)
{
    ++three_real_calls;
    return __real_cublasZgemm3m(handle, ta, tb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}

int main()
{
    using Complex = std::complex<double>;
    constexpr int leading = 129, columns = 131, inner = 4097;
    Complex *a = nullptr, *b = nullptr, *c = nullptr;
    check(cudaMalloc(&a, size_t(leading) * inner * sizeof(Complex)));
    check(cudaMalloc(&b, size_t(inner) * columns * sizeof(Complex)));
    check(cudaMalloc(&c, size_t(leading) * columns * sizeof(Complex)));
    check(cudaMemset(a, 0, size_t(leading) * inner * sizeof(Complex)));
    check(cudaMemset(b, 0, size_t(inner) * columns * sizeof(Complex)));
    hsolver::createGpuBlasHandle();
    const Complex alpha(1, 0), beta(0, 0);
    // A wide product must reach the portable 3M API on every CUDA device.
    // The adjacent ineligible row count must retain the ordinary GEMM path.
    for (int rows : {129, 127})
    {
        hsolver::gemm_op<Complex, base_device::DEVICE_GPU>()(
            nullptr, 'N', 'N', rows, columns, inner, &alpha, a, leading, b, inner, &beta, c, leading);
        check(cudaDeviceSynchronize());
        if (three_real_calls != 1)
        {
            std::cerr << "rows=" << rows << " expected_3m_calls=1 actual_3m_calls=" << three_real_calls << '\n';
            return 2;
        }
    }
    hsolver::destoryBLAShandle();
    check(cudaFree(c));
    check(cudaFree(b));
    check(cudaFree(a));
    std::cout << "eligible_and_fallback_dispatch=PASS three_real_calls=" << three_real_calls << '\n';
}
