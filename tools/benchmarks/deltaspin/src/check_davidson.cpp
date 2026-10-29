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
    const int rows = 96, bands = 6, leading = rows + 7;
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
            for (int state = 0; state < rows; ++state)
            {
                hamiltonian[row * rows + column] += basis[state * rows + row] * double(state + 1)
                                                   * std::conj(basis[state * rows + column]);
            }
        }
    }
    Complex* device = nullptr;
    base_device::memory::resize_memory_op<Complex, GPU>()(nullptr, device, leading * bands);
    std::vector<double> precondition(rows, 1.0), eigenvalues(bands), thresholds(bands, 1e-13);
    hsolver::diag_comm_info comm(MPI_COMM_SELF, 0, 1);
    double eigen_error = 0, residual_error = 0, orthogonality_error = 0;
    double projector_error = 0, projected_input_error = 0;
    std::vector<Complex> reference_basis;
    {
        hsolver::Diago_DavSubspace<Complex, GPU> solver(precondition, bands, rows, 4, 1e-13, 100, false, comm);
        for (int repetition = 0; repetition < 5; ++repetition)
        {
            std::vector<Complex> psi(leading * bands, Complex(7, -3));
            for (int band = 0; band < bands; ++band)
            {
                for (int row = 0; row < rows; ++row)
                {
                    psi[band * leading + row] = (basis[band * rows + row]
                        + 0.2 * basis[(band + bands) * rows + row]) / std::sqrt(1.04);
                }
            }
            if (repetition == 2 || repetition == 3)
            {
                const auto original = psi;
                for (int band = 0; band < bands; ++band)
                {
                    for (int row = 0; row < rows; ++row)
                    {
                        psi[band * leading + row] = 0;
                        for (int column = 0; column < bands; ++column)
                        {
                            const Complex coefficient = std::polar((1 + 0.03 * band) / std::sqrt(double(bands)),
                                                                   2 * pi * column * band / bands);
                            psi[band * leading + row] += original[column * leading + row] * coefficient;
                        }
                    }
                }
            }
            for (int row = 0; row < rows; ++row)
            {
                precondition[row] = 1.0 + std::min(repetition, 1) * (1.0 + row / 96.0);
            }
            base_device::memory::synchronize_memory_op<Complex, GPU, CPU>()(nullptr, nullptr, device,
                                                                           psi.data(), psi.size());
            Complex* initial_h = nullptr;
            Complex* initial_s = nullptr;
            std::vector<Complex> projected_h(bands * bands), projected_s(bands * bands);
            if (repetition == 3)
            {
                std::vector<Complex> hpsi(rows * bands);
                for (int band = 0; band < bands; ++band)
                {
                    for (int row = 0; row < rows; ++row)
                    {
                        for (int column = 0; column < rows; ++column)
                        {
                            hpsi[band * rows + row] += hamiltonian[row * rows + column]
                                                       * psi[band * leading + column];
                        }
                    }
                }
                for (int column = 0; column < bands; ++column)
                {
                    for (int row = 0; row < bands; ++row)
                    {
                        for (int basis_row = 0; basis_row < rows; ++basis_row)
                        {
                            const Complex left = std::conj(psi[row * leading + basis_row]);
                            projected_h[column * bands + row] += left * hpsi[column * rows + basis_row];
                            projected_s[column * bands + row] += left * psi[column * leading + basis_row];
                        }
                    }
                }
                base_device::memory::resize_memory_op<Complex, GPU>()(nullptr, initial_h, bands * bands);
                base_device::memory::resize_memory_op<Complex, GPU>()(nullptr, initial_s, bands * bands);
                base_device::memory::synchronize_memory_op<Complex, GPU, CPU>()(
                    nullptr, nullptr, initial_h, projected_h.data(), projected_h.size());
                base_device::memory::synchronize_memory_op<Complex, GPU, CPU>()(
                    nullptr, nullptr, initial_s, projected_s.data(), projected_s.size());
            }
            auto apply = [&](Complex* output, Complex* input, int input_bands, int input_rows, int first, int last) {
                std::vector<Complex> host_input(input_bands * input_rows), host_output((last - first + 1) * rows);
                base_device::memory::synchronize_memory_op<Complex, CPU, GPU>()(nullptr, nullptr,
                    host_input.data(), input, host_input.size());
                for (int band = first; band <= last; ++band)
                {
                    for (int row = 0; row < rows; ++row)
                    {
                        for (int column = 0; column < rows; ++column)
                        {
                            host_output[(band - first) * rows + row] += hamiltonian[row * rows + column]
                                                                    * host_input[band * input_rows + column];
                        }
                    }
                }
                base_device::memory::synchronize_memory_op<Complex, GPU, CPU>()(nullptr, nullptr,
                    output, host_output.data(), host_output.size());
            };
            const int iterations = solver.diag(apply, device, leading, eigenvalues.data(), thresholds.data(),
                                                false, initial_h, initial_s);
            if (iterations < 0) { return 2; }
            if (repetition == 3)
            {
                std::vector<Complex> retained_h(bands * bands), retained_s(bands * bands);
                base_device::memory::synchronize_memory_op<Complex, CPU, GPU>()(
                    nullptr, nullptr, retained_h.data(), initial_h, retained_h.size());
                base_device::memory::synchronize_memory_op<Complex, CPU, GPU>()(
                    nullptr, nullptr, retained_s.data(), initial_s, retained_s.size());
                for (int index = 0; index < bands * bands; ++index)
                {
                    projected_input_error = std::max(projected_input_error, std::abs(retained_h[index] - projected_h[index]));
                    projected_input_error = std::max(projected_input_error, std::abs(retained_s[index] - projected_s[index]));
                }
                base_device::memory::delete_memory_op<Complex, GPU>()(nullptr, initial_h);
                base_device::memory::delete_memory_op<Complex, GPU>()(nullptr, initial_s);
            }
            base_device::memory::synchronize_memory_op<Complex, CPU, GPU>()(nullptr, nullptr, psi.data(),
                                                                           device, psi.size());
            for (int band = 0; band < bands; ++band)
            {
                eigen_error = std::max(eigen_error, std::abs(eigenvalues[band] - (band + 1)));
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
            if (repetition == 1) { reference_basis = psi; }
            else if (repetition >= 2)
            {
                // The occupied-space projector is independent of eigenvector phases
                // and of invertible changes to the initial subspace basis.
                for (int row = 0; row < rows; ++row)
                {
                    for (int column = 0; column < rows; ++column)
                    {
                        Complex delta = 0;
                        for (int band = 0; band < bands; ++band)
                        {
                            delta += psi[band * leading + row] * std::conj(psi[band * leading + column])
                                   - reference_basis[band * leading + row] * std::conj(reference_basis[band * leading + column]);
                        }
                        projector_error = std::max(projector_error, std::abs(delta));
                    }
                }
            }
            std::cout << "repetition=" << repetition << " iterations=" << iterations << '\n';
        }
    }
    base_device::memory::delete_memory_op<Complex, GPU>()(nullptr, device);
    hsolver::destroyGpuSolverHandle();
    hsolver::destoryBLAShandle();
    MPI_Finalize();
    std::cout << "eigenvalue_error=" << eigen_error << " residual=" << residual_error
              << " orthogonality_error=" << orthogonality_error << '\n';
    std::cout << "projector_error=" << projector_error << " projected_input_error=" << projected_input_error << std::endl;
    return projected_input_error == 0 && projector_error < 1e-6 && eigen_error < 1e-10
           && residual_error < 2e-6 && orthogonality_error < 1e-10 ? 0 : 3;
}
