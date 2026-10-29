#include "module_hsolver/kernels/math_kernel_op.h"

#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <vector>

using Complex = std::complex<double>;

void require(bool condition, const char* message)
{
    if (!condition) { std::cerr << message << '\n'; std::exit(1); }
}
void cuda_check(cudaError_t status)
{
    require(status == cudaSuccess, cudaGetErrorString(status));
}
cublasOperation_t operation(char transpose)
{
    return transpose == 'N' ? CUBLAS_OP_N : transpose == 'T' ? CUBLAS_OP_T : CUBLAS_OP_C;
}
Complex entry(const std::vector<Complex>& data, int leading, int row, int column, char transpose)
{
    const Complex value = transpose == 'N' ? data[size_t(column) * leading + row]
                                          : data[size_t(row) * leading + column];
    return transpose == 'C' ? std::conj(value) : value;
}
void run(cublasHandle_t reference_handle, int m, int n, int k, char ta, char tb, bool cancellation)
{
    const int lda = (ta == 'N' ? m : k) + 3, ldb = (tb == 'N' ? k : n) + 5, ldc = m + 7;
    const Complex padding(51.25, -37.5);
    std::vector<Complex> a(size_t(lda) * (ta == 'N' ? k : m), padding);
    std::vector<Complex> b(size_t(ldb) * (tb == 'N' ? n : k), padding);
    std::vector<Complex> original(size_t(ldc) * n, padding), reference(original.size()), actual(original.size());
    for (int column = 0; column < (ta == 'N' ? k : m); ++column)
    for (int row = 0; row < (ta == 'N' ? m : k); ++row)
    {
        const size_t index = size_t(column) * lda + row;
        const double real = double((index * 37 + 13) % 997) / 997 - 0.5;
        const double imag = cancellation ? real * (1.0 + 1e-12) : double((index * 83 + 7) % 991) / 991 - 0.5;
        a[index] = Complex(real, imag) / std::sqrt(k);
    }
    for (int column = 0; column < (tb == 'N' ? n : k); ++column)
    for (int row = 0; row < (tb == 'N' ? k : n); ++row)
    {
        const size_t index = size_t(column) * ldb + row;
        const double real = double((index * 19 + 3) % 983) / 983 - 0.5;
        const double imag = cancellation ? real * (1.0 - 1e-12) : double((index * 47 + 11) % 977) / 977 - 0.5;
        b[index] = Complex(real, imag) / std::sqrt(k);
    }
    for (int column = 0; column < n; ++column)
    for (int row = 0; row < m; ++row) { original[size_t(column) * ldc + row] = Complex(row * 1e-4, -column * 2e-4); }
    Complex *da, *db, *dc, *dr;
    cuda_check(cudaMalloc(&da, a.size() * sizeof(Complex)));
    cuda_check(cudaMalloc(&db, b.size() * sizeof(Complex)));
    cuda_check(cudaMalloc(&dc, actual.size() * sizeof(Complex)));
    cuda_check(cudaMalloc(&dr, actual.size() * sizeof(Complex)));
    cuda_check(cudaMemcpy(da, a.data(), a.size() * sizeof(Complex), cudaMemcpyHostToDevice));
    cuda_check(cudaMemcpy(db, b.data(), b.size() * sizeof(Complex), cudaMemcpyHostToDevice));
    double reference_error = 0.0, extended_error = 0.0;
    const Complex alpha(0.7, -0.1);
    for (const Complex beta : {Complex{}, Complex(-0.3, 0.2)})
    {
        cuda_check(cudaMemcpy(dc, original.data(), actual.size() * sizeof(Complex), cudaMemcpyHostToDevice));
        cuda_check(cudaMemcpy(dr, original.data(), actual.size() * sizeof(Complex), cudaMemcpyHostToDevice));
        hsolver::gemm_op<Complex, base_device::DEVICE_GPU>()(
            nullptr, ta, tb, m, n, k, &alpha, da, lda, db, ldb, &beta, dc, ldc);
        require(cublasZgemm(reference_handle, operation(ta), operation(tb), m, n, k,
            reinterpret_cast<const cuDoubleComplex*>(&alpha), reinterpret_cast<const cuDoubleComplex*>(da), lda,
            reinterpret_cast<const cuDoubleComplex*>(db), ldb, reinterpret_cast<const cuDoubleComplex*>(&beta),
            reinterpret_cast<cuDoubleComplex*>(dr), ldc) == CUBLAS_STATUS_SUCCESS, "Reference GEMM failed");
        cuda_check(cudaMemcpy(actual.data(), dc, actual.size() * sizeof(Complex), cudaMemcpyDeviceToHost));
        cuda_check(cudaMemcpy(reference.data(), dr, reference.size() * sizeof(Complex), cudaMemcpyDeviceToHost));
        for (int column = 0; column < n; ++column)
        for (int row = 0; row < ldc; ++row)
        {
            const size_t index = size_t(column) * ldc + row;
            if (row >= m) { require(actual[index] == padding, "Output padding changed"); }
            require(std::isfinite(actual[index].real()) && std::isfinite(actual[index].imag()), "Nonfinite output");
            reference_error = std::max(reference_error, std::abs(actual[index] - reference[index]) / (1 + std::abs(reference[index])));
        }
        for (int j = 0; j < 9; ++j)
        for (int i = 0; i < 9; ++i)
        {
            const int row = i * (m - 1) / 8, column = j * (n - 1) / 8;
            std::complex<long double> sum{};
            for (int inner = 0; inner < k; ++inner)
            {
                sum += std::complex<long double>(entry(a, lda, row, inner, ta))
                       * std::complex<long double>(entry(b, ldb, inner, column, tb));
            }
            const auto expected = std::complex<long double>(alpha) * sum
                                  + std::complex<long double>(beta) * std::complex<long double>(original[size_t(column) * ldc + row]);
            extended_error = std::max(extended_error, double(std::abs(std::complex<long double>(actual[size_t(column) * ldc + row]) - expected)
                                                          / (1 + std::abs(expected))));
        }
    }
    std::vector<Complex> preserved(a.size());
    cuda_check(cudaMemcpy(preserved.data(), da, a.size() * sizeof(Complex), cudaMemcpyDeviceToHost));
    require(preserved == a, "First input changed");
    preserved.resize(b.size());
    cuda_check(cudaMemcpy(preserved.data(), db, b.size() * sizeof(Complex), cudaMemcpyDeviceToHost));
    require(preserved == b, "Second input changed");
    cuda_check(cudaFree(da)); cuda_check(cudaFree(db)); cuda_check(cudaFree(dc)); cuda_check(cudaFree(dr));
    require(reference_error < 1e-12 && extended_error < 1e-12, "GEMM error exceeds bound");
    std::cout << "PASS " << m << ' ' << n << ' ' << k << ' ' << ta << tb << " cancellation=" << cancellation
              << " reference_error=" << reference_error << " extended_error=" << extended_error << std::endl;
}
int main()
{
    hsolver::createGpuBlasHandle();
    cublasHandle_t reference_handle;
    require(cublasCreate(&reference_handle) == CUBLAS_STATUS_SUCCESS, "Reference handle creation failed");
    for (char ta : {'N', 'T', 'C'})
    for (char tb : {'N', 'T', 'C'}) { run(reference_handle, 129, 131, 4097, ta, tb, false); }
    run(reference_handle, 4097, 129, 65, 'N', 'N', false);
    run(reference_handle, 4097, 129, 65, 'N', 'T', true);
    run(reference_handle, 129, 4097, 65, 'C', 'N', false);
    run(reference_handle, 128, 128, 4096, 'C', 'N', true);
    run(reference_handle, 127, 128, 4096, 'C', 'N', false);
    run(reference_handle, 128, 127, 4096, 'C', 'N', false);
    run(reference_handle, 4096, 128, 63, 'N', 'N', false);
    run(reference_handle, 128, 128, 4095, 'C', 'N', false);
    for (char ta : {'N', 'T', 'C'})
    for (char tb : {'N', 'T', 'C'}) { run(reference_handle, 129, 64, 4097, ta, tb, true); }
    run(reference_handle, 128, 64, 4096, 'C', 'N', false);
    run(reference_handle, 127, 64, 4096, 'C', 'N', false);
    run(reference_handle, 128, 63, 4096, 'C', 'N', false);
    run(reference_handle, 128, 64, 4095, 'C', 'N', false);
    run(reference_handle, 68, 104, 4097, 'C', 'N', false);
    run(reference_handle, 4097, 104, 68, 'N', 'T', false);
    require(cublasDestroy(reference_handle) == CUBLAS_STATUS_SUCCESS, "Reference handle destruction failed");
    hsolver::destoryBLAShandle();
}
