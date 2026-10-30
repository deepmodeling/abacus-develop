#include "module_hsolver/kernels/fixed_overlap.h"
#include "module_hsolver/kernels/dngvd_op.h"
#include "module_hsolver/kernels/math_kernel_op.h"
#include <cuda_runtime_api.h>
#include <mpi.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <vector>

using Complex = std::complex<double>;
using Wide = std::complex<long double>;
using GPU = base_device::DEVICE_GPU;

void check(cudaError_t status)
{
    if (status != cudaSuccess) { std::cerr << cudaGetErrorString(status) << '\n'; std::exit(1); }
}

void require(bool condition, const char* message)
{
    if (!condition) { std::cerr << message << '\n'; std::exit(2); }
}

class DeviceBuffer
{
  public:
    explicit DeviceBuffer(size_t count) : values_(count + 4, Complex(91, -37))
    {
        check(cudaMalloc(&storage_, values_.size() * sizeof(Complex)));
        check(cudaMemcpy(storage_, values_.data(), values_.size() * sizeof(Complex), cudaMemcpyHostToDevice));
    }
    ~DeviceBuffer() { check(cudaFree(storage_)); }
    DeviceBuffer(const DeviceBuffer&) = delete;
    DeviceBuffer& operator=(const DeviceBuffer&) = delete;
    Complex* data() { return storage_ + 2; }
    void set(const std::vector<Complex>& values)
    {
        require(values.size() + 4 == values_.size(), "Device input length mismatch");
        check(cudaMemcpy(data(), values.data(), values.size() * sizeof(Complex), cudaMemcpyHostToDevice));
    }
    std::vector<Complex> get()
    {
        check(cudaMemcpy(values_.data(), storage_, values_.size() * sizeof(Complex), cudaMemcpyDeviceToHost));
        for (size_t index : {size_t(0), size_t(1), values_.size() - 2, values_.size() - 1})
        { require(values_[index] == Complex(91, -37), "Device buffer guard changed"); }
        return {values_.begin() + 2, values_.end() - 2};
    }
  private:
    Complex* storage_ = nullptr;
    std::vector<Complex> values_;
};

std::vector<Complex> multiply(const std::vector<Complex>& a, const std::vector<Complex>& b, int n)
{
    std::vector<Complex> product(size_t(n) * n);
    for (int column = 0; column < n; ++column)
    for (int row = 0; row < n; ++row)
    {
        Wide sum = 0;
        for (int inner = 0; inner < n; ++inner)
        { sum += Wide(a[size_t(inner) * n + row]) * Wide(b[size_t(column) * n + inner]); }
        product[size_t(column) * n + row] = Complex(sum);
    }
    return product;
}

void run(int n, int points)
{
    const int rows = n + 3;
    const size_t matrix = size_t(n) * n;
    const long double pi = std::acos(-1.0L);
    std::vector<std::vector<Complex>> h(points, std::vector<Complex>(matrix));
    auto s = h;
    std::vector<Complex> projectors(size_t(rows) * n);
    for (int column = 0; column < n; ++column)
    {
        projectors[size_t(column) * rows + column] = 1;
        for (int row = n; row < rows; ++row)
        { projectors[size_t(column) * rows + row] = std::polar(0.1 / std::sqrt(n), 0.17 * column * (row - n + 1)); }
    }
    for (int point = 0; point < points; ++point)
    for (int column = 0; column < n; ++column)
    for (int row = 0; row < n; ++row)
    {
        Wide hs = 0, ss = 0;
        for (int state = 0; state < n; ++state)
        {
            const long double overlap = std::exp(-std::log(1000.0L) * state / n);
            const long double energy = 0.3L * point + 0.01L * state - 0.4L;
            const Wide phase = std::polar(1.0L / n, 2 * pi * (row - column) * state / n);
            hs += phase * overlap * energy;
            ss += phase * overlap;
        }
        h[point][size_t(column) * n + row] = Complex(hs);
        s[point][size_t(column) * n + row] = Complex(ss);
    }
    DeviceBuffer input_h(matrix), input_s(matrix), input_r(projectors.size());
    DeviceBuffer shifted(matrix), rotated(projectors.size()), reference_vectors(matrix);
    DeviceBuffer reference_h(matrix), reference_s(matrix);
    double* device_eigenvalues = nullptr;
    check(cudaMalloc(&device_eigenvalues, n * sizeof(double)));
    input_r.set(projectors);
    hsolver::FixedOverlapSubspace cache(points, n, rows);
    require(cache.available(), "Small fixed-overlap fixture could not allocate");
    for (int point = 0; point < points; ++point)
    {
        input_h.set(h[point]); input_s.set(s[point]);
        cache.prepare(point, input_h.data(), input_s.data(), input_r.data());
        require(input_h.get() == h[point] && input_s.get() == s[point]
                && input_r.get() == projectors, "Preparation changed an input");
    }
    double maximum_residual = 0, maximum_orthogonality = 0, maximum_spectrum = 0, maximum_projection = 0;
    for (int point = points - 1; point >= 0; --point)
    for (double field : {-0.003, 0.0, 0.0027})
    {
        std::vector<Complex> transformed_h(matrix), transformed_r(projectors.size());
        check(cudaMemcpy(transformed_h.data(), cache.hamiltonian(point), matrix * sizeof(Complex), cudaMemcpyDeviceToHost));
        check(cudaMemcpy(transformed_r.data(), cache.projectors(point), transformed_r.size() * sizeof(Complex), cudaMemcpyDeviceToHost));
        auto candidate_h = transformed_h;
        auto physical_h = h[point];
        for (int column = 0; column < n; ++column)
        for (int row = 0; row < n; ++row)
        {
            Wide physical_shift = 0, transformed_shift = 0;
            for (int inner = 0; inner < rows; ++inner)
            {
                const long double weight = field * (inner % 3 - 1);
                physical_shift += weight * std::conj(Wide(projectors[size_t(row) * rows + inner]))
                                  * Wide(projectors[size_t(column) * rows + inner]);
                transformed_shift += weight * std::conj(Wide(transformed_r[size_t(row) * rows + inner]))
                                     * Wide(transformed_r[size_t(column) * rows + inner]);
            }
            physical_h[size_t(column) * n + row] += Complex(physical_shift);
            candidate_h[size_t(column) * n + row] += Complex(transformed_shift);
        }
        shifted.set(candidate_h);
        std::vector<double> values(n), reference_values(n);
        cache.diagonalize(shifted.data(), cache.projectors(point), rotated.data(), values.data());
        require(shifted.get() == candidate_h, "Diagonalization changed the projected input");
        const auto projections = rotated.get();
        std::vector<Complex> vectors(matrix);
        for (int column = 0; column < n; ++column)
        { std::copy_n(projections.data() + size_t(column) * rows, n, vectors.data() + size_t(column) * n); }
        for (int column = 0; column < n; ++column)
        for (int row = n; row < rows; ++row)
        {
            Wide expected = 0;
            for (int inner = 0; inner < n; ++inner)
            { expected += Wide(projectors[size_t(inner) * rows + row]) * Wide(vectors[size_t(column) * n + inner]); }
            const double error = double(std::abs(expected - Wide(projections[size_t(column) * rows + row])));
            require(std::isfinite(error), "Non-finite projected eigenvector");
            maximum_projection = std::max(maximum_projection, error);
        }
        const auto hv = multiply(physical_h, vectors, n);
        const auto sv = multiply(s[point], vectors, n);
        long double norm_h = 0, norm_s = 0;
        for (size_t index = 0; index < matrix; ++index)
        {
            norm_h += std::norm(Wide(physical_h[index]));
            norm_s += std::norm(Wide(s[point][index]));
        }
        norm_h = std::sqrt(norm_h);
        norm_s = std::sqrt(norm_s);
        for (int column = 0; column < n; ++column)
        {
            long double numerator = 0, vector_norm = 0;
            for (int row = 0; row < n; ++row)
            {
                const size_t index = size_t(column) * n + row;
                numerator += std::norm(hv[index] - values[column] * sv[index]);
                vector_norm += std::norm(Wide(vectors[index]));
                Wide dot = 0;
                for (int inner = 0; inner < n; ++inner)
                { dot += std::conj(Wide(vectors[size_t(row) * n + inner])) * Wide(sv[size_t(column) * n + inner]); }
                maximum_orthogonality = std::max(maximum_orthogonality, double(std::abs(dot - Wide(row == column))));
            }
            // Normwise backward error remains defined for zero eigenvalues.
            const long double denominator = (norm_h + std::abs(values[column]) * norm_s) * std::sqrt(vector_norm);
            const double residual = double(std::sqrt(numerator) / denominator);
            require(std::isfinite(values[column]) && std::isfinite(residual), "Non-finite eigenpair");
            maximum_residual = std::max(maximum_residual, residual);
        }
        reference_h.set(physical_h); reference_s.set(s[point]);
        hsolver::dngvd_op<Complex, GPU>()(nullptr, n, n, reference_h.data(), reference_s.data(),
                                         device_eigenvalues, reference_vectors.data());
        check(cudaMemcpy(reference_values.data(), device_eigenvalues, n * sizeof(double), cudaMemcpyDeviceToHost));
        for (int band = 0; band < n; ++band)
        { maximum_spectrum = std::max(maximum_spectrum, std::abs(values[band] - reference_values[band])); }
        std::vector<Complex> preserved(matrix), preserved_r(projectors.size());
        check(cudaMemcpy(preserved.data(), cache.hamiltonian(point), matrix * sizeof(Complex), cudaMemcpyDeviceToHost));
        check(cudaMemcpy(preserved_r.data(), cache.projectors(point), preserved_r.size() * sizeof(Complex), cudaMemcpyDeviceToHost));
        require(preserved == transformed_h && preserved_r == transformed_r, "Cached matrices changed during response");
        reference_vectors.get(); reference_h.get(); reference_s.get();
    }
    std::cout << "bands=" << n << " points=" << points << " backward_error=" << maximum_residual
              << " orthogonality=" << maximum_orthogonality << " spectrum=" << maximum_spectrum
              << " projection=" << maximum_projection << std::endl;
    require(maximum_residual < 1e-12 && maximum_orthogonality < 1e-10 && maximum_spectrum < 1e-10
            && maximum_projection < 1e-10,
            "Fixed-overlap eigenpair checks failed");
    check(cudaFree(device_eigenvalues));
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    hsolver::createGpuBlasHandle();
    hsolver::createGpuSolverHandle();
    size_t available = 0, total = 0;
    check(cudaMemGetInfo(&available, &total));
    {
        const int points = int(total / (16 * (16 + 19) * sizeof(Complex))) + 1024;
        hsolver::FixedOverlapSubspace unavailable(points, 16, 19);
        require(!unavailable.available(), "Oversized optional cache did not report unavailable");
    }
    for (int generation = 0; generation < 2; ++generation)
    {
        for (int bands : {7, 31, 65, 257}) { run(bands, generation + 1); }
    }
    hsolver::destroyGpuSolverHandle();
    hsolver::destoryBLAShandle();
    MPI_Finalize();
}
