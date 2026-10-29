#include "module_hamilt_lcao/module_deltaspin/spin_constrain.h"
#include "module_hsolver/kernels/math_kernel_op.h"
#include "module_parameter/parameter.h"
#include <ATen/core/tensor.h>
#include <mpi.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <iostream>
#include <numeric>
#include <vector>

using Complex = std::complex<double>;
using Cpu = base_device::DEVICE_CPU;
using Gpu = base_device::DEVICE_GPU;

double run(int atoms, int bands, int npol, int sign, int field, bool measure)
{
    auto& sc = SpinConstrain<Complex>::getScInstance();
    sc.set_atomCounts({{0, atoms}});
    sc.set_nspin(npol == 2 ? 4 : 2);
    sc.set_npol(npol);
    sc.zero_Mi();
    std::vector<int> projectors(atoms);
    std::vector<ModuleBase::Vector3<double>> lambda(atoms);
    for (int atom = 0; atom < atoms; ++atom)
    {
        projectors[atom] = measure ? 27 : (atom * 3) % 8;
        lambda[atom] = {field == 0 || field == 1 ? 0.0 : 0.03 * (atom + 1),
                        field == 0 || field == 2 ? 0.0 : -0.05 * (atom + 1),
                        field == 0 || field == 3 ? 0.0 : 0.07 * (atom + 1)};
    }
    const int nkb = std::accumulate(projectors.begin(), projectors.end(), 0);
    const int rows = nkb * npol;
    std::vector<Complex> input(rows * bands), initial(bands * bands);
    for (int i = 0; i < int(input.size()); ++i)
    {
        input[i] = Complex(std::sin(0.13 * i), std::cos(0.17 * i)) / std::sqrt(rows);
    }
    for (int j = 0; j < bands; ++j)
    {
        for (int i = 0; i < bands; ++i)
        {
            initial[j * bands + i] = {0.01 * (i + j + 1), 0.02 * (i - j)};
        }
    }
    auto expected = initial;
    const_cast<Input_para&>(PARAM.inp).device = "cpu";
    sc.calculate_delta_hcc(expected.data(), input.data(), lambda.data(), bands, nkb, projectors.data(), sign);
    double analytic_error = 0;
    if (!measure)
    {
        for (int j = 0; j < bands; ++j)
        {
            for (int i = 0; i < bands; ++i)
            {
                std::complex<long double> value = initial[j * bands + i];
                int offset = 0;
                for (int atom = 0; atom < atoms; ++atom)
                {
                    const std::complex<long double> cross(lambda[atom].x, -lambda[atom].y);
                    for (int ip = 0; ip < projectors[atom]; ++ip)
                    {
                        const int index = offset + ip;
                        const std::complex<long double> a = input[i * rows + index], b = input[j * rows + index];
                        if (npol == 1) { value += std::conj(a) * b * static_cast<long double>(lambda[atom].z * sign); }
                        else
                        {
                            const std::complex<long double> c = input[i * rows + index + nkb];
                            const std::complex<long double> d = input[j * rows + index + nkb];
                            value += std::conj(a) * (static_cast<long double>(lambda[atom].z) * b + cross * d)
                                     + std::conj(c) * (std::conj(cross) * b - static_cast<long double>(lambda[atom].z) * d);
                        }
                    }
                    offset += projectors[atom];
                }
                analytic_error = std::max(analytic_error, double(std::abs(value - std::complex<long double>(expected[j * bands + i]))));
            }
        }
    }
    ct::Tensor coeff(ct::DataType::DT_COMPLEX_DOUBLE, ct::DeviceType::GpuDevice, {bands, rows});
    ct::Tensor matrix(ct::DataType::DT_COMPLEX_DOUBLE, ct::DeviceType::GpuDevice, {bands, bands});
    base_device::memory::synchronize_memory_op<Complex, Gpu, Cpu>()(nullptr, nullptr, coeff.data<Complex>(), input.data(), input.size());
    const_cast<Input_para&>(PARAM.inp).device = "gpu";
    std::vector<Complex> actual(initial.size()), preserved(input.size());
    auto execute = [&]()
    {
        base_device::memory::synchronize_memory_op<Complex, Gpu, Cpu>()(nullptr, nullptr, matrix.data<Complex>(), initial.data(), initial.size());
        sc.calculate_delta_hcc(matrix.data<Complex>(), coeff.data<Complex>(), lambda.data(), bands, nkb, projectors.data(), sign);
        base_device::memory::synchronize_memory_op<Complex, Cpu, Gpu>()(nullptr, nullptr, actual.data(), matrix.data<Complex>(), actual.size());
    };
    execute();
    double error = 0, hermitian = 0;
    for (int j = 0; j < bands; ++j)
    {
        for (int i = 0; i < bands; ++i)
        {
            const int k = j * bands + i;
            if (!std::isfinite(std::abs(actual[k]))) { MPI_Abort(MPI_COMM_WORLD, 2); }
            error = std::max(error, std::abs(actual[k] - expected[k]) / (1 + std::abs(expected[k])));
            hermitian = std::max(hermitian, std::abs(actual[k] - std::conj(actual[i * bands + j])));
        }
    }
    base_device::memory::synchronize_memory_op<Complex, Cpu, Gpu>()(nullptr, nullptr, preserved.data(), coeff.data<Complex>(), preserved.size());
    if (input != preserved || error > 1e-12 || analytic_error > 1e-12 || hermitian > 1e-10)
    {
        std::cerr << "Delta Hamiltonian mismatch " << error << ' ' << analytic_error << ' ' << hermitian << '\n';
        MPI_Abort(MPI_COMM_WORLD, 3);
    }
    double elapsed_ms = 0;
    if (measure)
    {
        std::vector<double> times;
        for (int repeat = 0; repeat < 3; ++repeat)
        {
            const auto start = std::chrono::steady_clock::now();
            for (int i = 0; i < 10; ++i) { execute(); }
            times.push_back(std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count() / 10);
        }
        std::sort(times.begin(), times.end());
        elapsed_ms = times[1];
    }
    std::cout << "atoms=" << atoms << " bands=" << bands << " npol=" << npol << " sign=" << sign
              << " field=" << field << " relative_error=" << error << " analytic_error=" << analytic_error
              << " hermitian=" << hermitian << " elapsed_ms=" << elapsed_ms << std::endl;
    return error;
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    hsolver::createGpuBlasHandle();
    for (int npol : {1, 2})
    {
        for (int sign : {-1, 1})
        {
            for (int field = 0; field < 5; ++field) { run(5, 7, npol, sign, field, false); }
        }
    }
    run(16, 308, 2, 1, 4, true);
    run(54, 1037, 2, 1, 4, true);
    hsolver::destoryBLAShandle();
    MPI_Finalize();
}
