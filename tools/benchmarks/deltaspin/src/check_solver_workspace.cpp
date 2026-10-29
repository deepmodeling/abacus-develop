#include "module_hsolver/kernels/dngvd_op.h"
#include <cuda_runtime_api.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <vector>

using Complex = std::complex<double>;
using GPU = base_device::DEVICE_GPU;

void check(cudaError_t status)
{
    if (status != cudaSuccess) { std::cerr << cudaGetErrorString(status) << '\n'; std::exit(1); }
}

void run(int n, bool invalid_overlap)
{
    const int leading = n + 7;
    const size_t count = size_t(leading) * n, bytes = count * sizeof(Complex);
    std::vector<Complex> a(count, Complex(99, -17)), b(a), vectors(a), preserved(a);
    std::vector<double> expected(n), values(n);
    const Complex alpha(0.15, 0.08), beta(0.025, -0.013);
    for (int j = 0; j < n; ++j)
    {
        for (int i = 0; i < n; ++i) { a[j * leading + i] = b[j * leading + i] = 0; }
        a[j * leading + j] = 2;
        b[j * leading + j] = 1;
        a[j * leading + (j + 1) % n] = alpha;
        a[j * leading + (j + n - 1) % n] = std::conj(alpha);
        b[j * leading + (j + 1) % n] = beta;
        b[j * leading + (j + n - 1) % n] = std::conj(beta);
        const Complex phase = std::polar(1.0, -2 * std::acos(-1.0) * j / n);
        expected[j] = (2 + 2 * (alpha * phase).real()) / (1 + 2 * (beta * phase).real());
    }
    if (invalid_overlap) { b[0] = -1; }
    std::sort(expected.begin(), expected.end());
    Complex *da = nullptr, *db = nullptr, *dv = nullptr;
    double* dw = nullptr;
    for (auto** pointer : {&da, &db, &dv}) { check(cudaMalloc(reinterpret_cast<void**>(pointer), bytes)); }
    check(cudaMalloc(reinterpret_cast<void**>(&dw), n * sizeof(double)));
    check(cudaMemcpy(da, a.data(), bytes, cudaMemcpyHostToDevice));
    check(cudaMemcpy(db, b.data(), bytes, cudaMemcpyHostToDevice));
    int status = -1;
    hsolver::dngvd_op<Complex, GPU>()(nullptr, n, leading, da, db, dw, dv, &status);
    if (invalid_overlap)
    {
        if (status <= n) { std::cerr << "Missing indefinite-overlap status: " << status << '\n'; std::exit(2); }
        std::cout << "n=" << n << " invalid_overlap_status=" << status << std::endl;
    }
    else
    {
        if (status != 0) { std::cerr << "Solver status " << status << '\n'; std::exit(3); }
        check(cudaMemcpy(values.data(), dw, n * sizeof(double), cudaMemcpyDeviceToHost));
        check(cudaMemcpy(vectors.data(), dv, bytes, cudaMemcpyDeviceToHost));
        check(cudaMemcpy(preserved.data(), da, bytes, cudaMemcpyDeviceToHost));
        if (preserved != a) { std::cerr << "Input matrix changed\n"; std::exit(4); }
        double eigen_error = 0, residual = 0;
        for (int column = 0; column < n; ++column)
        {
            if (!std::isfinite(values[column])) { std::exit(5); }
            eigen_error = std::max(eigen_error, std::abs(values[column] - expected[column]));
            const auto* v = vectors.data() + column * leading;
            double numerator = 0, denominator = 0;
            for (int i = 0; i < n; ++i)
            {
                const Complex av = 2.0 * v[i] + alpha * v[(i + n - 1) % n] + std::conj(alpha) * v[(i + 1) % n];
                const Complex bv = v[i] + beta * v[(i + n - 1) % n] + std::conj(beta) * v[(i + 1) % n];
                numerator += std::norm(av - values[column] * bv);
                denominator += std::norm(av) + values[column] * values[column] * std::norm(bv);
            }
            const double relative = std::sqrt(numerator / denominator);
            if (!std::isfinite(relative)) { std::exit(6); }
            residual = std::max(residual, relative);
            for (int i = n; i < leading; ++i)
            {
                if (v[i] != a[column * leading + i]) { std::cerr << "Padding changed\n"; std::exit(7); }
            }
        }
        std::cout << "n=" << n << " eigen_error=" << eigen_error << " relative_residual=" << residual << std::endl;
        if (eigen_error > 1e-12 || residual > 1e-12) { std::exit(8); }
    }
    for (auto* pointer : {da, db, dv}) { check(cudaFree(pointer)); }
    check(cudaFree(dw));
}

int main()
{
    for (int lifecycle = 0; lifecycle < 2; ++lifecycle)
    {
        hsolver::createGpuSolverHandle();
        for (int n : {31, 308, 615, 1230, 7, 615}) { run(n, false); }
        run(31, true);
        run(31, false);
        hsolver::destroyGpuSolverHandle();
    }
}
