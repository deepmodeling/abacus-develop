#include "module_basis/module_pw/pw_basis_k.h"
#include "module_elecstate/kernels/elecstate_op.h"

#include <mpi.h>
#include <cuda_runtime.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <random>
#include <vector>

using Complex = std::complex<double>;

void check_cuda(cudaError_t status)
{
    if (status != cudaSuccess) { std::cerr << cudaGetErrorString(status) << '\n'; std::exit(1); }
}

void require(bool condition, const char* message)
{
    if (!condition) { std::cerr << message << '\n'; std::exit(2); }
}

void run(ModulePW::PW_Basis_K& basis, double cutoff, double& maximum, double& relative, double& charge_error)
{
    basis.initmpi(1, 0, MPI_COMM_WORLD);
    basis.initgrids(30.0, ModuleBase::Matrix3(1, 0.1, 0, 0, 1, 0.15, 0.05, 0, 1), 200, 200, 200);
    ModuleBase::Vector3<double> points[2] = {{0.25, -0.125, 0.375}, {-0.2, 0.15, -0.1}};
    basis.initparameters(false, cutoff, 2, points);
    basis.setuptransform();
    basis.collect_local_pw();
    basis.get_ig2ixyz_k();
    const int stride = basis.npwk_max + 7, size = basis.nxyz;
    const int compact_size = basis.ft.fused_veff_grid_size();
    require(basis.ft.fused_veff_stick_size() > 0, "Stick transform not available");
    std::mt19937 generator(1234);
    std::uniform_real_distribution<double> random(-1, 1);
    std::vector<Complex> input(2 * stride), dense(size), preserved(input.size());
    for (auto& value : input) { value = Complex(random(generator), random(generator)) / std::sqrt(stride); }
    std::vector<double> expected(4 * size), actual(4 * size);
    Complex *device_input = nullptr, *grid = nullptr, *real[2] = {}, *compact[2] = {};
    double *reference_density = nullptr, *density = nullptr;
    check_cuda(cudaMalloc(&device_input, input.size() * sizeof(Complex)));
    check_cuda(cudaMalloc(&grid, size * sizeof(Complex)));
    check_cuda(cudaMalloc(&reference_density, expected.size() * sizeof(double)));
    check_cuda(cudaMalloc(&density, expected.size() * sizeof(double)));
    for (int spin = 0; spin < 2; ++spin)
    {
        check_cuda(cudaMalloc(&real[spin], size * sizeof(Complex)));
        check_cuda(cudaMalloc(&compact[spin], compact_size * sizeof(Complex)));
    }
    check_cuda(cudaMemcpy(device_input, input.data(), input.size() * sizeof(Complex), cudaMemcpyHostToDevice));
    double* channels[4] = {reference_density, reference_density + size,
                           reference_density + 2 * size, reference_density + 3 * size};
    for (int ik = 0; ik < 2; ++ik)
    {
        long double norm = 0;
        for (int spin = 0; spin < 2; ++spin)
        {
            std::fill(dense.begin(), dense.end(), Complex{});
            for (int ig = 0; ig < basis.npwk[ik]; ++ig)
            {
                const auto value = input[spin * stride + ig];
                dense[basis.ig2ixyz_k_[ik * basis.npwk_max + ig]] = value;
                norm += std::norm(value);
            }
            check_cuda(cudaMemcpy(grid, dense.data(), size * sizeof(Complex), cudaMemcpyHostToDevice));
            basis.ft.fft3D_backward(static_cast<const base_device::DEVICE_GPU*>(nullptr), grid, real[spin]);
        }
        check_cuda(cudaMemset(reference_density, 0, expected.size() * sizeof(double)));
        check_cuda(cudaMemset(density, 0, expected.size() * sizeof(double)));
        double weight_sum = 0;
        for (double weight : {0.2, -0.03, 0.0, 1e-20})
        {
            elecstate::elecstate_pw_op<double, base_device::DEVICE_GPU>()(
                nullptr, true, false, size, weight, channels, real[0], real[1]);
            basis.add_spinor_density(device_input, weight, density, ik, stride, compact[0], compact[1]);
            weight_sum += weight;
        }
        check_cuda(cudaMemcpy(expected.data(), reference_density, expected.size() * sizeof(double), cudaMemcpyDeviceToHost));
        check_cuda(cudaMemcpy(actual.data(), density, actual.size() * sizeof(double), cudaMemcpyDeviceToHost));
        long double delta_norm = 0, expected_norm = 0, charge = 0;
        for (size_t i = 0; i < expected.size(); ++i)
        {
            require(std::isfinite(actual[i]), "Non-finite density");
            const double delta = actual[i] - expected[i];
            maximum = std::max(maximum, std::abs(delta));
            delta_norm += delta * delta;
            expected_norm += expected[i] * expected[i];
            if (i < size) { charge += actual[i]; }
        }
        relative = std::max(relative, std::sqrt(double(delta_norm / expected_norm)));
        charge_error = std::max(charge_error, std::abs(double(charge / (weight_sum * size * norm)) - 1.0));
    }
    check_cuda(cudaMemcpy(preserved.data(), device_input, input.size() * sizeof(Complex), cudaMemcpyDeviceToHost));
    require(preserved == input, "Input wavefunctions changed");
    for (int spin = 0; spin < 2; ++spin) { check_cuda(cudaFree(real[spin])); check_cuda(cudaFree(compact[spin])); }
    check_cuda(cudaFree(device_input)); check_cuda(cudaFree(grid));
    check_cuda(cudaFree(reference_density)); check_cuda(cudaFree(density));
    std::cout << "cutoff=" << cutoff << " max_error=" << maximum << " relative_error=" << relative
              << " charge_error=" << charge_error << std::endl;
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    double maximum = 0, relative = 0, charge_error = 0;
    size_t free_bytes[2] = {}, total = 0;
    for (int repeat = 0; repeat < 2; ++repeat)
    {
        {
            ModulePW::PW_Basis_K basis("gpu", "double");
            run(basis, 50.0, maximum, relative, charge_error);
            run(basis, 100.0, maximum, relative, charge_error);
        }
        check_cuda(cudaDeviceSynchronize());
        check_cuda(cudaMemGetInfo(&free_bytes[repeat], &total));
    }
    require(maximum < 1e-12 && relative < 1e-12 && charge_error < 1e-12, "Density comparison failed");
    require(free_bytes[1] >= free_bytes[0], "Device memory retained after repeated basis lifetimes");
    std::cout << "retained_bytes=" << static_cast<long long>(free_bytes[0]) - static_cast<long long>(free_bytes[1]) << std::endl;
    MPI_Finalize();
}
