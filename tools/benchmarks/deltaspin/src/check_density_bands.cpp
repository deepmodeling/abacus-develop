#include "module_basis/module_pw/kernels/fused_veff.h"
#include "module_basis/module_pw/pw_basis_k.h"

#include <cuda_runtime.h>
#include <mpi.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

using Complex = std::complex<double>;

void check_cuda(cudaError_t status)
{
    if (status != cudaSuccess) { std::cerr << cudaGetErrorString(status) << '\n'; std::exit(1); }
}

void run(int length, int nz, int bands, bool zeros)
{
    ModulePW::PW_Basis_K basis("gpu", "double");
    basis.initmpi(1, 0, MPI_COMM_WORLD);
    basis.initgrids((length / 72.0) * 2 * 2.866 / 0.529177210903,
        ModuleBase::Matrix3(1, 0, 0, 0, 1, 0, 0, 0, double(nz) / length), length, length, nz);
    ModuleBase::Vector3<double> points[2] = {{0.25, -0.125, 0.375}, {-0.2, 0.15, -0.1}};
    basis.initparameters(false, 100.0, 2, points, 1, true);
    basis.setuptransform();
    basis.collect_local_pw();
    basis.get_ig2ixyz_k();
    if (basis.ft.fused_veff_stick_size() <= 0) { std::exit(2); }
    bands = std::min(bands, basis.ft.spinor_batch_size());
    const int compact = basis.ft.fused_veff_grid_size();
    const int leading = basis.npwk_max + 7;
    const size_t density_size = size_t(4) * basis.nxyz;
    const size_t input_size = size_t(2) * bands * leading;
    double weights[4] = {0.3, -0.07, 0.13, 0.04};
    if (zeros) { weights[1] = 0; weights[2] = 1e-30; }
    Complex *input = nullptr, *slabs = nullptr;
    double* density_allocation = nullptr;
    constexpr int guard = 31;
    std::vector<double> guarded(density_size + 2 * guard, 47.25);
    check_cuda(cudaMalloc(&input, input_size * sizeof(Complex)));
    check_cuda(cudaMalloc(&slabs, size_t(2) * bands * compact * sizeof(Complex)));
    check_cuda(cudaMalloc(&density_allocation, guarded.size() * sizeof(double)));
    check_cuda(cudaMemcpy(density_allocation, guarded.data(), guarded.size() * sizeof(double), cudaMemcpyHostToDevice));
    auto* density = density_allocation + guard;
    auto* up = slabs;
    auto* down = slabs + bands * compact;
    std::vector<Complex> host_input(input_size), preserved(input_size);
    std::vector<double> initial(density_size), reference(density_size), actual(density_size);
    for (size_t i = 0; i < density_size; ++i) { initial[i] = 0.01 * std::sin(0.01 * i); }
    std::mt19937 generator(7311);
    std::uniform_real_distribution<double> uniform(-0.5, 0.5);
    for (int ik = 0; ik < 2; ++ik)
    {
        long double expected_charge = 0;
        const int rows = basis.npwk[ik];
        for (int band = 0; band < bands; ++band)
        for (int spin = 0; spin < 2; ++spin)
        for (int row = 0; row < leading; ++row)
        {
            const Complex value = row < rows ? Complex(uniform(generator), uniform(generator)) / std::sqrt(double(rows))
                                             : Complex(17.25, -23.75);
            host_input[(2 * band + spin) * leading + row] = value;
            if (row < rows) { expected_charge += weights[band] * std::norm(value); }
        }
        check_cuda(cudaMemcpy(input, host_input.data(), input_size * sizeof(Complex), cudaMemcpyHostToDevice));
        auto execute = [&](bool batch) {
            if (!batch)
            {
                for (int band = 0; band < bands; ++band)
                {
                    if (weights[band] == 0) { continue; }
                    basis.add_spinor_density(input + size_t(2 * band) * leading, weights[band],
                                             density, ik, leading, up, down);
                }
                return;
            }
            basis.add_spinor_density(input, weights, bands, density, ik, leading, up, down);
        };
        for (int method = 0; method < 2; ++method)
        {
            check_cuda(cudaMemcpy(density, initial.data(), density_size * sizeof(double), cudaMemcpyHostToDevice));
            execute(method);
            check_cuda(cudaGetLastError());
            check_cuda(cudaMemcpy(method ? actual.data() : reference.data(), density,
                                  density_size * sizeof(double), cudaMemcpyDeviceToHost));
        }
        double error = 0;
        long double charge = 0;
        for (size_t i = 0; i < density_size; ++i)
        {
            const double difference = std::abs(actual[i] - reference[i]);
            if (!std::isfinite(difference)) { std::exit(3); }
            error = std::max(error, difference);
            if (i < size_t(basis.nxyz)) { charge += actual[i] - initial[i]; }
        }
        const double charge_error = std::abs(double(charge / basis.nxyz - expected_charge));
        if (std::memcmp(actual.data(), reference.data(), density_size * sizeof(double)) != 0 || charge_error > 2e-12)
        {
            std::cerr << "Density error " << error << ", charge error " << charge_error << '\n';
            std::exit(4);
        }
        cudaEvent_t start, stop;
        check_cuda(cudaEventCreate(&start));
        check_cuda(cudaEventCreate(&stop));
        double samples[2][5] = {};
        for (int round = 0; round < 5; ++round)
        for (int order = 0; order < 2; ++order)
        {
            const int method = (order + round) % 2;
            execute(method);
            check_cuda(cudaEventRecord(start));
            for (int repetition = 0; repetition < 15; ++repetition) { execute(method); }
            check_cuda(cudaEventRecord(stop));
            check_cuda(cudaEventSynchronize(stop));
            float elapsed = 0;
            check_cuda(cudaEventElapsedTime(&elapsed, start, stop));
            samples[method][round] = elapsed / 15;
        }
        for (auto& method : samples) { std::sort(method, method + 5); }
        check_cuda(cudaMemcpy(preserved.data(), input, input_size * sizeof(Complex), cudaMemcpyDeviceToHost));
        if (preserved != host_input) { std::cerr << "Input changed" << '\n'; std::exit(5); }
        check_cuda(cudaMemcpy(guarded.data(), density_allocation, guarded.size() * sizeof(double), cudaMemcpyDeviceToHost));
        for (int i = 0; i < guard; ++i)
        {
            if (guarded[i] != 47.25 || guarded[guard + density_size + i] != 47.25)
            {
                std::cerr << "Density guard overwritten" << '\n'; std::exit(6);
            }
        }
        std::cout << std::setprecision(12) << "length=" << length << " nz=" << nz << " bands=" << bands << " zeros=" << zeros << " k=" << ik
                  << " npw=" << rows << " max_error=" << error << " charge_error=" << charge_error
                  << " sequential_ms=" << samples[0][2] << " batched_ms=" << samples[1][2]
                  << " speedup=" << samples[0][2] / samples[1][2]
                  << " additional_slab_bytes=" << size_t(2) * (bands - 1) * compact * sizeof(Complex) << std::endl;
        check_cuda(cudaEventDestroy(start));
        check_cuda(cudaEventDestroy(stop));
    }
    check_cuda(cudaFree(input));
    check_cuda(cudaFree(slabs));
    check_cuda(cudaFree(density_allocation));
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    for (int bands = 1; bands <= 4; ++bands)
    {
        run(36, 36, bands, false);
        run(72, 144, bands, false);
    }
    for (int length : {48, 64, 81, 108, 125, 144, 200})
    {
        run(length, length == 125 ? 250 : length, 4, false);
        run(length, length == 125 ? 250 : length, 4, true);
    }
    if (ModulePW::fused_fft(17) != nullptr || ModulePW::fused_fft(257) != nullptr)
    {
        std::cerr << "Unsupported FFT dispatch" << '\n'; return 7;
    }
    MPI_Finalize();
}
