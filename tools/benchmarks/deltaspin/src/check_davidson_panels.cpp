#include "module_hsolver/diago_dav_subspace.h"
#include "module_hsolver/kernels/dngvd_op.h"
#include "module_hsolver/kernels/math_kernel_op.h"
#include "module_base/module_device/memory_op.h"
#include <mpi.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <vector>

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    hsolver::createGpuBlasHandle();
    hsolver::createGpuSolverHandle();
    using Complex = std::complex<double>;
    using GPU = base_device::DEVICE_GPU;
    using CPU = base_device::DEVICE_CPU;
    const int rows = 1024, bands = 257, leading = rows + 7;
    const double pi = std::acos(-1.0);
    std::vector<Complex> basis(rows * rows), hamiltonian(rows * rows);
    for (int column = 0; column < rows; ++column)
    {
        for (int row = 0; row < rows; ++row)
        {
            basis[column * rows + row] = std::polar(1.0 / std::sqrt(double(rows)),
                                                  2 * pi * row * column / rows);
        }
    }
    for (int row = 0; row < rows; ++row)
    {
        for (int column = 0; column < rows; ++column)
        {
            const long double angle = 2 * std::acos(-1.0L) * (row - column) / rows;
            const std::complex<long double> phase = std::polar(1.0L, angle);
            const std::complex<long double> value = row == column
                ? std::complex<long double>((rows + 1.0L) / (2 * rows), 0)
                : -1.0L / (static_cast<long double>(rows) * (1.0L - phase));
            hamiltonian[row * rows + column] = Complex(value.real(), value.imag());
        }
    }
    Complex* device = nullptr;
    base_device::memory::resize_memory_op<Complex, GPU>()(nullptr, device, leading * bands);
    Complex* device_hamiltonian = nullptr;
    base_device::memory::resize_memory_op<Complex, GPU>()(nullptr, device_hamiltonian, rows * rows);
    base_device::memory::synchronize_memory_op<Complex, GPU, CPU>()(
        nullptr, nullptr, device_hamiltonian, hamiltonian.data(), hamiltonian.size());
    const Complex one(1, 0), zero(0, 0);
    std::vector<double> precondition(rows, 1.0), eigenvalues(bands), thresholds(bands, 1e-13);
    hsolver::diag_comm_info comm(MPI_COMM_SELF, 0, 1);
    double eigen_error = 0, residual_error = 0, orthogonality_error = 0;
    {
        hsolver::Diago_DavSubspace<Complex, GPU> solver(precondition, bands, rows, 4, 1e-13, 100, false, comm);
        for (int repetition = 0; repetition < 2; ++repetition)
        {
            std::vector<Complex> psi(leading * bands, Complex(7, -3));
            for (int band = 0; band < bands; ++band)
            {
                for (int row = 0; row < rows; ++row)
                {
                    psi[band * leading + row] = (basis[band * rows + row]
                        + 0.2 * basis[(band + bands) * rows + row]
                        + 0.05 * basis[(band + 2 * bands) * rows + row]) / std::sqrt(1.0425);
                }
            }
            for (int row = 0; row < rows; ++row) { precondition[row] = 1.0 + repetition; }
            base_device::memory::synchronize_memory_op<Complex, GPU, CPU>()(nullptr, nullptr, device,
                                                                           psi.data(), psi.size());
            auto apply = [&](Complex* output, Complex* input, int, int input_rows, int first, int last) {
                hsolver::gemm_op<Complex, GPU>()(nullptr, 'T', 'N', rows, last - first + 1, rows,
                    &one, device_hamiltonian, rows, input + first * input_rows, input_rows, &zero, output, rows);
            };
            const int iterations = solver.diag(apply, device, leading, eigenvalues.data(), thresholds.data(), false);
            if (iterations < 0) { return 2; }
            base_device::memory::synchronize_memory_op<Complex, CPU, GPU>()(nullptr, nullptr, psi.data(),
                                                                           device, psi.size());
            for (int band = 0; band < bands; ++band)
            {
                eigen_error = std::max(eigen_error, std::abs(eigenvalues[band] - double(band + 1) / rows));
                double residual = 0;
                for (int row = 0; row < rows; ++row)
                {
                    Complex value = -eigenvalues[band] * psi[band * leading + row];
                    for (int column = 0; column < rows; ++column)
                    {
                        value += hamiltonian[row * rows + column] * psi[band * leading + column];
                    }
                    residual += std::norm(value);
                }
                residual_error = std::max(residual_error, std::sqrt(residual));
                for (int other = 0; other < bands; ++other)
                {
                    Complex overlap = 0;
                    for (int row = 0; row < rows; ++row)
                    {
                        overlap += std::conj(psi[band * leading + row]) * psi[other * leading + row];
                    }
                    orthogonality_error = std::max(orthogonality_error, std::abs(overlap - double(band == other)));
                }
            }
            std::cout << "repetition=" << repetition << " iterations=" << iterations << '\n';
        }
    }
    base_device::memory::delete_memory_op<Complex, GPU>()(nullptr, device);
    base_device::memory::delete_memory_op<Complex, GPU>()(nullptr, device_hamiltonian);
    hsolver::destroyGpuSolverHandle();
    hsolver::destoryBLAShandle();
    MPI_Finalize();
    std::cout << "eigenvalue_error=" << eigen_error << " residual=" << residual_error
              << " orthogonality_error=" << orthogonality_error << '\n';
    return eigen_error < 1e-10 && residual_error < 1e-6 && orthogonality_error < 1e-10 ? 0 : 3;
}
