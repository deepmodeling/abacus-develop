#include "module_basis/module_pw/pw_basis_k.h"
#include "module_basis/module_pw/kernels/fused_veff.h"
#include "module_hamilt_pw/hamilt_pwdft/operator_pw/veff_pw.h"
#include "module_hamilt_pw/hamilt_pwdft/kernels/veff_op.h"

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

double compare(const std::vector<Complex>& actual, const std::vector<Complex>& expected)
{
    double error = 0;
    for (size_t i = 0; i < actual.size(); ++i)
    {
        require(std::isfinite(actual[i].real()) && std::isfinite(actual[i].imag()), "Non-finite output");
        error = std::max(error, std::abs(actual[i] - expected[i]) / std::max(1.0, std::abs(expected[i])));
    }
    return error;
}

std::vector<Complex> reference(const ModulePW::PW_Basis_K& basis, int ik, int stride, int npol, int bands,
                               const std::vector<Complex>& input, const double* potential,
                               Complex* grid, Complex* up, Complex* down)
{
    const base_device::DEVICE_GPU* context = nullptr;
    const size_t bytes = basis.nxyz * sizeof(Complex);
    std::vector<Complex> dense(basis.nxyz), result(input.size(), 0.0);
    Complex* real[2] = {up, down};
    for (int band = 0; band < bands; ++band)
    {
        for (int spin = 0; spin < npol; ++spin)
        {
            std::fill(dense.begin(), dense.end(), Complex(0, 0));
            for (int ig = 0; ig < basis.npwk[ik]; ++ig)
            {
                dense[basis.ig2ixyz_k_[ik * basis.npwk_max + ig]] = input[(band * npol + spin) * stride + ig];
            }
            check_cuda(cudaMemcpy(grid, dense.data(), bytes, cudaMemcpyHostToDevice));
            basis.ft.fft3D_backward(context, grid, real[spin]);
        }
        if (npol == 2)
        {
            const double* channels[4] = {potential, potential + basis.nxyz,
                                          potential + 2 * basis.nxyz, potential + 3 * basis.nxyz};
            hamilt::veff_pw_op<double, base_device::DEVICE_GPU>()(context, basis.nxyz, up, down, channels);
        }
        else { hamilt::veff_pw_op<double, base_device::DEVICE_GPU>()(context, basis.nxyz, up, potential); }
        for (int spin = 0; spin < npol; ++spin)
        {
            basis.ft.fft3D_forward(context, real[spin], grid);
            check_cuda(cudaMemcpy(dense.data(), grid, bytes, cudaMemcpyDeviceToHost));
            for (int ig = 0; ig < basis.npwk[ik]; ++ig)
            {
                result[(band * npol + spin) * stride + ig]
                    = dense[basis.ig2ixyz_k_[ik * basis.npwk_max + ig]] / double(basis.nxyz);
            }
        }
    }
    return result;
}

void run(ModulePW::PW_Basis_K& basis, int nx, int ny, int nz, bool xprime,
         double& max_error, double& max_hermitian, double cutoff = 50.0, int bands = 2)
{
    basis.initmpi(1, 0, MPI_COMM_WORLD);
    const double length = nx * 0.15;
    basis.initgrids(length, ModuleBase::Matrix3(1, 0.1, 0, 0, double(ny) / nx, 0.15, 0.05, 0, double(nz) / nx), nx, ny, nz);
    ModuleBase::Vector3<double> points[2] = {{0.25, -0.125, 0.375}, {-0.2, 0.15, -0.1}};
    basis.initparameters(false, cutoff, 2, points, 1, xprime);
    basis.setuptransform();
    basis.collect_local_pw();
    basis.get_ig2ixyz_k();
    const int stride = basis.npwk_max + 7;
    const int count = 2 * bands * stride;
    const int compact_size = basis.ft.fused_veff_grid_size();
    require((compact_size > 0) == (xprime && ModulePW::fused_fft(nx) != nullptr), "Unexpected fused grid eligibility");
    const int stick_size = basis.ft.fused_veff_stick_size();
    std::mt19937 generator(1234);
    std::uniform_real_distribution<double> random(-1, 1);
    std::vector<Complex> input(count), initial(count), expected(count), actual(count);
    const double scale = 1 / std::sqrt(basis.npwk_max);
    for (int i = 0; i < count; ++i)
    {
        input[i] = Complex(random(generator), random(generator)) * scale;
        initial[i] = Complex(random(generator), random(generator)) * scale;
    }
    std::vector<double> potential(4 * basis.nxyz);
    for (size_t i = 0; i < potential.size(); ++i) { potential[i] = 0.1 * std::sin(0.017 * i) + 0.03; }
    Complex *device_input = nullptr, *device_output = nullptr, *grid = nullptr, *up = nullptr, *down = nullptr;
    Complex *compact_up = nullptr, *compact_down = nullptr;
    double* device_potential = nullptr;
    const size_t wave_bytes = count * sizeof(Complex), grid_bytes = basis.nxyz * sizeof(Complex);
    check_cuda(cudaMalloc(&device_input, wave_bytes));
    check_cuda(cudaMalloc(&device_output, wave_bytes));
    check_cuda(cudaMalloc(&grid, grid_bytes));
    check_cuda(cudaMalloc(&up, grid_bytes));
    check_cuda(cudaMalloc(&down, grid_bytes));
    check_cuda(cudaMalloc(&device_potential, potential.size() * sizeof(double)));
    check_cuda(cudaMemcpy(device_input, input.data(), wave_bytes, cudaMemcpyHostToDevice));
    check_cuda(cudaMemcpy(device_potential, potential.data(), potential.size() * sizeof(double), cudaMemcpyHostToDevice));
    if (compact_size > 0)
    {
        check_cuda(cudaMalloc(&compact_up, basis.ft.spinor_batch_size() * compact_size * sizeof(Complex)));
        check_cuda(cudaMalloc(&compact_down, basis.ft.spinor_batch_size() * compact_size * sizeof(Complex)));
    }
    const int spins[2] = {0, 0};
    {
        hamilt::Veff<hamilt::OperatorPW<Complex, base_device::DEVICE_GPU>> op(spins, device_potential, 4, basis.nxyz, &basis);
        for (int ik = 0; ik < 2; ++ik)
        {
            op.init(ik);
            for (int npol : {2, 1, 2})
            {
                const auto pure = reference(basis, ik, stride, npol, bands, input, device_potential, grid, up, down);
                check_cuda(cudaMemcpy(device_output, initial.data(), wave_bytes, cudaMemcpyHostToDevice));
                for (int repeat = 1; repeat <= 2; ++repeat)
                {
                    op.act(bands * npol, npol * stride, npol, device_input, device_output, basis.npwk[ik]);
                    check_cuda(cudaGetLastError());
                    check_cuda(cudaMemcpy(actual.data(), device_output, wave_bytes, cudaMemcpyDeviceToHost));
                    for (int i = 0; i < count; ++i) { expected[i] = initial[i] + double(repeat) * pure[i]; }
                    max_error = std::max(max_error, compare(actual, expected));
                }
                if (npol == 2 && compact_size > 0)
                {
                    check_cuda(cudaMemcpy(device_output, initial.data(), wave_bytes, cudaMemcpyHostToDevice));
                    for (int band = 0; band < bands;)
                    {
                        const int remaining = bands - band;
                        const int batch = std::min(remaining, basis.ft.spinor_batch_size());
                        basis.add_spinor_potential(device_input + band * 2 * stride, device_potential,
                                                   device_output + band * 2 * stride, ik, stride, compact_up, compact_down, batch);
                        band += batch;
                    }
                    check_cuda(cudaMemcpy(actual.data(), device_output, wave_bytes, cudaMemcpyDeviceToHost));
                    for (int i = 0; i < count; ++i) { expected[i] = initial[i] + pure[i]; }
                    max_error = std::max(max_error, compare(actual, expected));
                    Complex a_vb = 0, va_b = 0;
                    for (int spin = 0; spin < 2; ++spin)
                    {
                        for (int ig = 0; ig < basis.npwk[ik]; ++ig)
                        {
                            const int a = spin * stride + ig, b = (2 + spin) * stride + ig;
                            a_vb += std::conj(input[a]) * (actual[b] - initial[b]);
                            va_b += std::conj(actual[a] - initial[a]) * input[b];
                        }
                    }
                    max_hermitian = std::max(max_hermitian, std::abs(a_vb - va_b) / std::max(1.0, std::abs(a_vb)));
                }
            }
        }
    }
    check_cuda(cudaMemcpy(actual.data(), device_input, wave_bytes, cudaMemcpyDeviceToHost));
    require(actual == input, "Input wavefunctions changed");
    check_cuda(cudaFree(device_input));
    check_cuda(cudaFree(device_output));
    check_cuda(cudaFree(grid));
    check_cuda(cudaFree(up));
    check_cuda(cudaFree(down));
    check_cuda(cudaFree(compact_up));
    check_cuda(cudaFree(compact_down));
    check_cuda(cudaFree(device_potential));
    std::cout << "grid=" << nx << ',' << ny << ',' << nz << " xprime=" << xprime
              << " bands=" << bands << " cutoff=" << cutoff << " compact=" << compact_size << " sticks=" << stick_size
              << " error=" << max_error << " hermitian=" << max_hermitian << std::endl;
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    double max_error = 0, max_hermitian = 0;
    size_t free_bytes[2] = {}, total_bytes = 0;
    for (int repeat = 0; repeat < 2; ++repeat)
    {
        {
            ModulePW::PW_Basis_K basis("gpu", "double");
            for (int length = 16; length <= 256; ++length)
            {
                if (ModulePW::fused_fft(length) != nullptr)
                {
                    run(basis, length, length, 31, true, max_error, max_hermitian, 50.0, 3);
                }
            }
            run(basis, 79, 83, 89, true, max_error, max_hermitian, 50.0, 3);
            run(basis, 50, 72, 80, true, max_error, max_hermitian, 50.0, 3);
            run(basis, 72, 108, 125, true, max_error, max_hermitian, 100.0, 19);
            run(basis, 125, 81, 90, true, max_error, max_hermitian, 100.0, 3);
            run(basis, 72, 72, 72, true, max_error, max_hermitian, 100.0, 19);
            run(basis, 72, 72, 144, true, max_error, max_hermitian, 100.0, 33);
            run(basis, 108, 108, 108, true, max_error, max_hermitian, 100.0, 17);
            run(basis, 72, 73, 74, true, max_error, max_hermitian, 50.0, 3);
            run(basis, 72, 72, 72, false, max_error, max_hermitian, 100.0, 3);
            run(basis, 200, 81, 83, true, max_error, max_hermitian);
            run(basis, 200, 200, 200, true, max_error, max_hermitian, 100.0, 7);
            run(basis, 200, 81, 83, false, max_error, max_hermitian);
            run(basis, 64, 65, 66, true, max_error, max_hermitian);
            run(basis, 200, 200, 200, true, max_error, max_hermitian, 50.0, 5);
            run(basis, 200, 200, 200, true, max_error, max_hermitian, 160.0, 2);
        }
        check_cuda(cudaDeviceSynchronize());
        check_cuda(cudaMemGetInfo(&free_bytes[repeat], &total_bytes));
    }
    require(free_bytes[1] >= free_bytes[0], "GPU memory retained after repeated basis lifetimes");
    std::cout << "retained_bytes_after_repeat=" << static_cast<long long>(free_bytes[0])
                 - static_cast<long long>(free_bytes[1]) << std::endl;
    require(max_error < 1e-12 && max_hermitian < 1e-12, "Spinor potential comparison failed");
    MPI_Finalize();
}
