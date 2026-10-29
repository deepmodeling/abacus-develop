#include "module_hsolver/kernels/math_kernel_op.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <vector>

void check_cuda(cudaError_t status)
{
    if (status != cudaSuccess) { std::cerr << cudaGetErrorString(status) << '\n'; std::exit(1); }
}

void run(int rows, int basis, int bands)
{
    using Complex = std::complex<double>;
    const int leading = rows + 7;
    const size_t input_size = size_t(leading) * basis;
    const size_t output_size = size_t(leading) * bands;
    std::vector<Complex> psi(input_size + output_size), hpsi(input_size);
    std::vector<int> selected(bands);
    std::vector<double> energies(bands);
    for (size_t i = 0; i < psi.size(); ++i)
    {
        psi[i] = Complex(double((37 * i + 13) % 997) / 997 - 0.5,
                         double((83 * i + 17) % 991) / 991 - 0.5);
    }
    for (size_t i = 0; i < hpsi.size(); ++i) { hpsi[i] = Complex(0.7, 0.3) * psi[i]; }
    std::vector<Complex> expected(psi.begin() + input_size, psi.end());
    for (int band = 0; band < bands; ++band)
    {
        selected[band] = band * basis / bands;
        energies[band] = (band % 13 - 6) * 0.137;
        for (int row = 0; row < rows; ++row)
        {
            const size_t source = size_t(selected[band]) * leading + row;
            expected[size_t(band) * leading + row] = hpsi[source] - energies[band] * psi[source];
        }
    }
    Complex *d_psi = nullptr, *d_hpsi = nullptr;
    int* d_selected = nullptr;
    double* d_energies = nullptr;
    check_cuda(cudaMalloc(&d_psi, psi.size() * sizeof(Complex)));
    check_cuda(cudaMalloc(&d_hpsi, hpsi.size() * sizeof(Complex)));
    check_cuda(cudaMalloc(&d_selected, bands * sizeof(int)));
    check_cuda(cudaMalloc(&d_energies, bands * sizeof(double)));
    check_cuda(cudaMemcpy(d_psi, psi.data(), psi.size() * sizeof(Complex), cudaMemcpyHostToDevice));
    check_cuda(cudaMemcpy(d_hpsi, hpsi.data(), hpsi.size() * sizeof(Complex), cudaMemcpyHostToDevice));
    check_cuda(cudaMemcpy(d_selected, selected.data(), bands * sizeof(int), cudaMemcpyHostToDevice));
    check_cuda(cudaMemcpy(d_energies, energies.data(), bands * sizeof(double), cudaMemcpyHostToDevice));
    hsolver::selected_basis_residual_gpu(rows, bands, leading, d_psi, d_hpsi,
                                        d_selected, d_energies, d_psi + input_size);
    std::vector<Complex> actual(psi.size()), actual_hpsi(hpsi.size());
    check_cuda(cudaMemcpy(actual.data(), d_psi, psi.size() * sizeof(Complex), cudaMemcpyDeviceToHost));
    check_cuda(cudaMemcpy(actual_hpsi.data(), d_hpsi, hpsi.size() * sizeof(Complex), cudaMemcpyDeviceToHost));
    if (!std::equal(psi.begin(), psi.begin() + input_size, actual.begin()) || actual_hpsi != hpsi)
    {
        std::cerr << "Input modified\n"; std::exit(2);
    }
    double maximum = 0;
    for (size_t i = 0; i < output_size; ++i)
    {
        const Complex value = actual[input_size + i];
        if (!std::isfinite(value.real()) || !std::isfinite(value.imag())) { std::exit(3); }
        if (i % leading >= size_t(rows) && value != expected[i]) { std::exit(4); }
        maximum = std::max(maximum, std::abs(value - expected[i]));
    }
    check_cuda(cudaFree(d_psi)); check_cuda(cudaFree(d_hpsi));
    check_cuda(cudaFree(d_selected)); check_cuda(cudaFree(d_energies));
    if (maximum > 1e-15) { std::cerr << "Residual mismatch " << maximum << '\n'; std::exit(5); }
    std::cout << "rows=" << rows << " basis=" << basis << " bands=" << bands << " max_error=" << maximum << '\n';
}

int main()
{
    hsolver::selected_basis_residual_gpu(0, 7, 7, nullptr, nullptr, nullptr, nullptr, nullptr);
    hsolver::selected_basis_residual_gpu(7, 0, 7, nullptr, nullptr, nullptr, nullptr, nullptr);
    run(19, 11, 7);
    run(257, 37, 37);
    run(1025, 37, 13);
    run(1025, 37, 1);
}
