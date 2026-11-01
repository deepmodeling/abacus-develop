#include "module_hsolver/davidson_action.h"
#include "module_psi/kpoint_storage.h"
#include "module_base/global_variable.h"

#include <mpi.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <vector>

using Complex = std::complex<double>;
using GPU = base_device::DEVICE_GPU;

void require(bool condition, const char* message)
{
    if (!condition) { throw std::runtime_error(message); }
}

void check_actions(int rows, int columns, bool host_backed, bool host_failure)
{
    const int outputs = 7, leading = rows + 11, coefficient_leading = columns + 5;
    std::vector<Complex> basis(size_t(rows) * columns), coefficients(coefficient_leading * outputs);
    for (size_t i = 0; i < basis.size(); ++i)
    {
        basis[i] = Complex(std::sin(0.13 * i), std::cos(0.17 * i));
    }
    for (size_t i = 0; i < coefficients.size(); ++i)
    {
        coefficients[i] = Complex(std::sin(0.23 * i), std::cos(0.07 * i)) / double(columns);
    }
    Complex *device_basis = nullptr, *device_coefficients = nullptr, *device_output = nullptr;
    cudaErrcheck(cudaMalloc(reinterpret_cast<void**>(&device_basis), basis.size() * sizeof(Complex)));
    cudaErrcheck(cudaMalloc(reinterpret_cast<void**>(&device_coefficients), coefficients.size() * sizeof(Complex)));
    cudaErrcheck(cudaMalloc(reinterpret_cast<void**>(&device_output), size_t(leading) * outputs * sizeof(Complex)));
    cudaErrcheck(cudaMemcpy(device_basis, basis.data(), basis.size() * sizeof(Complex), cudaMemcpyHostToDevice));
    cudaErrcheck(cudaMemcpy(device_coefficients, coefficients.data(), coefficients.size() * sizeof(Complex), cudaMemcpyHostToDevice));
    double maximum_error = 0;
    {
        hsolver::DavidsonAction<Complex, GPU> action(rows, columns, host_backed);
        require(action.host_backed() == (host_backed && !host_failure), "Requested action storage is unavailable");
        int applied = 0;
        auto identity = [&](Complex* output, Complex* input, int input_columns, int input_rows, int first, int last) {
            require(input == device_basis && input_columns == columns && input_rows == rows,
                    "Hamiltonian input layout changed");
            require(first == applied && last >= first && last < columns, "Hamiltonian band coverage changed");
            cudaErrcheck(cudaMemcpy(output, input + size_t(first) * rows,
                                   size_t(last - first + 1) * rows * sizeof(Complex), cudaMemcpyDeviceToDevice));
            applied = last + 1;
        };
        const int split = columns / 3;
        action.apply(identity, device_basis, columns, 0, split);
        action.apply(identity, device_basis, columns, split, columns - split);
        require(applied == columns, "Hamiltonian did not cover every band");
        for (int first = 0; first < columns; first += action.panel_columns())
        {
            const int width = std::min(action.panel_columns(), columns - first);
            std::vector<Complex> actual(size_t(rows) * width);
            cudaErrcheck(cudaMemcpy(actual.data(), action.panel(first, width),
                                   actual.size() * sizeof(Complex), cudaMemcpyDeviceToHost));
            require(std::memcmp(actual.data(), basis.data() + size_t(first) * rows,
                                actual.size() * sizeof(Complex)) == 0, "Action panel changed stored values");
        }
        action.store(3, 2, device_basis);
        std::copy_n(basis.begin(), 2 * rows, basis.begin() + 3 * rows);
        for (int inner : {1, columns / 2, columns})
        {
            const Complex padding(19.0, -7.0);
            std::vector<Complex> actual(leading * outputs, padding);
            cudaErrcheck(cudaMemcpy(device_output, actual.data(), actual.size() * sizeof(Complex), cudaMemcpyHostToDevice));
            action.multiply(inner, outputs, device_coefficients, coefficient_leading, device_output, leading);
            cudaErrcheck(cudaMemcpy(actual.data(), device_output, actual.size() * sizeof(Complex), cudaMemcpyDeviceToHost));
            for (int column = 0; column < outputs; ++column)
            {
                for (int row = 0; row < rows; ++row)
                {
                    Complex expected = 0;
                    for (int k = 0; k < inner; ++k)
                    {
                        expected += basis[size_t(k) * rows + row] * coefficients[column * coefficient_leading + k];
                    }
                    const double error = std::abs(actual[column * leading + row] - expected);
                    require(std::isfinite(error), "Action product is nonfinite");
                    maximum_error = std::max(maximum_error, error);
                }
                for (int row = rows; row < leading; ++row)
                {
                    require(actual[column * leading + row] == padding, "Action product overwrote output padding");
                }
            }
        }
    }
    cudaErrcheck(cudaFree(device_basis));
    cudaErrcheck(cudaFree(device_coefficients));
    cudaErrcheck(cudaFree(device_output));
    require(maximum_error < 1e-12, "Action product differs from independent CPU contraction");
    std::cout << "ACTION_STORAGE_OK rows=" << rows << " columns=" << columns
              << " host_backed=" << host_backed << " error=" << maximum_error << '\n';
}

void check_kpoints(int points, int polarizations, bool enabled, bool host_failure)
{
    GlobalV::NPOL = polarizations;
    const int bands = 5, leading = 37;
    const size_t block = size_t(bands) * leading;
    std::vector<int> ngk(points);
    for (int k = 0; k < points; ++k) { ngk[k] = leading - 2 - k % 5; }
    psi::Psi<Complex, GPU> wavefunctions(points, bands, leading, ngk.data());
    std::vector<Complex> expected(points * block);
    for (size_t i = 0; i < expected.size(); ++i) { expected[i] = Complex(i * 0.125, -0.25 * i); }
    cudaErrcheck(cudaMemcpy(wavefunctions.get_pointer(), expected.data(), expected.size() * sizeof(Complex), cudaMemcpyHostToDevice));
    wavefunctions.fix_k(points - 1);
    wavefunctions.fix_b(2);
    {
        psi::KPointStorage<Complex, GPU> unused(wavefunctions, points, enabled);
    }
    require(wavefunctions.get_current_k() == points - 1 && wavefunctions.get_current_b() == 2,
            "Unused storage lifetime changed selection");
    {
        psi::KPointStorage<Complex, GPU> storage(wavefunctions, points, enabled);
        require(storage.active() == (enabled && points > 1 && !host_failure), "Unexpected k-point storage selection");
        for (int point : {points - 1, 0, points / 2, points - 1})
        {
            storage.load(point);
            wavefunctions.fix_k(point);
            require(wavefunctions.get_current_k() == point, "Global k-point index was lost");
            require(wavefunctions.get_current_nbas() == (polarizations == 2 ? leading : ngk[point]),
                    "Active k-point basis size changed");
            std::vector<Complex> actual(block);
            cudaErrcheck(cudaMemcpy(actual.data(), wavefunctions.get_pointer(), block * sizeof(Complex), cudaMemcpyDeviceToHost));
            require(std::memcmp(actual.data(), expected.data() + point * block, block * sizeof(Complex)) == 0,
                    "K-point load changed wavefunctions or padding");
            for (int band = 0; band < bands; ++band)
            {
                actual[band * leading + 2] += Complex(1.0, -2.0);
                expected[point * block + band * leading + 2] += Complex(1.0, -2.0);
            }
            cudaErrcheck(cudaMemcpy(wavefunctions.get_pointer(), actual.data(), block * sizeof(Complex), cudaMemcpyHostToDevice));
            storage.store(point);
        }
        storage.restore();
        storage.restore();
    }
    require(wavefunctions.get_nk() == points && wavefunctions.get_current_k() == points - 1,
            "K-point restore changed shape or selection");
    std::vector<Complex> actual(expected.size());
    Complex* base = wavefunctions.get_pointer() - wavefunctions.get_psi_bias();
    cudaErrcheck(cudaMemcpy(actual.data(), base, actual.size() * sizeof(Complex), cudaMemcpyDeviceToHost));
    require(std::memcmp(actual.data(), expected.data(), actual.size() * sizeof(Complex)) == 0,
            "K-point restore changed the full wavefunction");
    psi::Psi<Complex, GPU> wrapped(base, points, bands, leading, ngk.data());
    psi::KPointStorage<Complex, GPU> borrowed(wrapped, points, true);
    require(!borrowed.active() && wrapped.get_pointer() == base, "Borrowed storage was reallocated");
    psi::Psi<Complex, GPU> band_first(points, bands, leading, ngk.data(), false);
    psi::KPointStorage<Complex, GPU> interleaved(band_first, points, true);
    require(!interleaved.active(), "Band-first storage was staged");
    psi::KPointStorage<Complex, GPU> mismatched(wavefunctions, points + 1, true);
    require(!mismatched.active(), "Incomplete k-point storage was staged");
    std::cout << "KPOINT_STORAGE_OK points=" << points << " npol=" << polarizations
              << " enabled=" << enabled << '\n';
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    hsolver::createGpuBlasHandle();
    const bool host_failure = argc == 2 && std::string(argv[1]) == "--expect-host-failure";
    for (bool host_backed : {false, true})
    {
        for (int rows : {17, 131})
        {
            for (int columns : {9, 145}) { check_actions(rows, columns, host_backed, host_failure); }
        }
    }
    for (int points : {1, 2, 3, 8})
    {
        for (int polarizations : {1, 2})
        {
            for (bool enabled : {false, true}) { check_kpoints(points, polarizations, enabled, host_failure); }
        }
    }
    hsolver::destoryBLAShandle();
    MPI_Finalize();
}
