#include "module_hsolver/diago_iter_assist.h"
#include "module_hsolver/kernels/math_kernel_op.h"
#include "module_hamilt_pw/hamilt_pwdft/kernels/ekinetic_op.h"
#include "module_base/global_variable.h"
#include "module_base/module_device/memory_op.h"
#include <cuda_runtime_api.h>
#include <mpi.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <type_traits>
#include <vector>

using Complex = std::complex<double>;
using CPU = base_device::DEVICE_CPU;
using GPU = base_device::DEVICE_GPU;

void cuda_check(cudaError_t status)
{
    if (status != cudaSuccess) { std::cerr << cudaGetErrorString(status) << '\n'; std::exit(1); }
}

template <typename Device>
class DiagonalOperator : public hamilt::Operator<Complex, Device>
{
public:
    using Info = typename hamilt::Operator<Complex, Device>::hpsi_info;
    DiagonalOperator(int rows, double offset, double slope)
    {
        std::vector<double> diagonal(rows);
        for (int row = 0; row < rows; ++row) { diagonal[row] = offset + slope * row / rows; }
        base_device::memory::resize_memory_op<double, Device>()(nullptr, diagonal_, rows);
        base_device::memory::synchronize_memory_op<double, Device, CPU>()(
            nullptr, nullptr, diagonal_, diagonal.data(), rows);
    }
    ~DiagonalOperator() override
    {
        base_device::memory::delete_memory_op<double, Device>()(nullptr, diagonal_);
    }
    void apply(const Complex* input, Complex* output, int leading, int rows, int bands) const
    {
        base_device::memory::set_memory_op<Complex, Device>()(nullptr, output, 0, size_t(leading) * bands);
        hamilt::ekinetic_pw_op<double, Device>()(nullptr, bands, rows, leading, 1.0,
                                                diagonal_, output, input);
    }
    Info hPsi(Info& input) const override
    {
        const auto& psi = *std::get<0>(input);
        const int leading = psi.get_nbasis(), rows = psi.get_current_nbas(), bands = psi.get_nbands();
        apply(psi.get_pointer(), std::get<2>(input), leading, rows, bands);
        return input;
    }
private:
    double* diagonal_ = nullptr;
};

template <typename Device>
class DiagonalHamiltonian : public hamilt::Hamilt<Complex, Device>
{
public:
    DiagonalHamiltonian(int rows, bool identity) : action_(rows, 2.0, 1.0),
        overlap_(rows, 1.0, identity ? 0.0 : 0.25) { this->ops = &action_; }
    void sPsi(const Complex* input, Complex* output, int leading, int rows, int bands) const override
    {
        ++overlap_calls;
        overlap_.apply(input, output, leading, rows, bands);
    }
    mutable int overlap_calls = 0;
private:
    DiagonalOperator<Device> action_;
    DiagonalOperator<Device> overlap_;
};

template <typename Device>
void run(int rows, int leading, int bands, bool identity)
{
    psi::Psi<Complex, Device> psi(1, bands, leading, &rows);
    psi.fix_k(0);
    const size_t count = size_t(leading) * bands, matrix_count = size_t(bands) * bands;
    std::vector<Complex> input(count), expected_action(count), host_h(matrix_count), host_s(matrix_count);
    for (size_t i = 0; i < count; ++i)
    {
        input[i] = Complex(double((i * 37 + 3) % 997) / 997.0,
                           double((i * 83 + 17) % 991) / 991.0) / std::sqrt(double(rows));
        const int row = i % leading;
        expected_action[i] = input[i] * (row < rows ? 2.0 + double(row) / rows : 0.0);
    }
    base_device::memory::synchronize_memory_op<Complex, Device, CPU>()(
        nullptr, nullptr, psi.get_pointer(), input.data(), count);
    Complex *h = nullptr, *s = nullptr;
    base_device::memory::resize_memory_op<Complex, Device>()(nullptr, h, matrix_count);
    base_device::memory::resize_memory_op<Complex, Device>()(nullptr, s, matrix_count);
    DiagonalHamiltonian<Device> model(rows, identity);
    hsolver::DiagoIterAssist<Complex, Device>::cal_hs_subspace(&model, psi, h, s);
    base_device::memory::synchronize_memory_op<Complex, CPU, Device>()(
        nullptr, nullptr, host_h.data(), h, matrix_count);
    base_device::memory::synchronize_memory_op<Complex, CPU, Device>()(
        nullptr, nullptr, host_s.data(), s, matrix_count);
    if (rows < 100)
    {
        double error = 0;
        for (int column = 0; column < bands; ++column)
        {
            for (int row = 0; row < bands; ++row)
            {
                Complex expected_h = 0, expected_s = 0;
                for (int inner = 0; inner < rows; ++inner)
                {
                    const auto product = std::conj(input[size_t(row) * leading + inner])
                                         * input[size_t(column) * leading + inner];
                    expected_h += product * (2.0 + double(inner) / rows);
                    expected_s += product * (identity ? 1.0 : 1.0 + 0.25 * inner / rows);
                }
                error = std::max(error, std::abs(expected_h - host_h[size_t(column) * bands + row]));
                error = std::max(error, std::abs(expected_s - host_s[size_t(column) * bands + row]));
            }
        }
        if (error > 1e-12) { std::cerr << "Analytic projection error " << error << '\n'; std::exit(2); }
    }
    std::vector<Complex> pageable(count), actual(matrix_count), retained_input(count);
    Complex* pinned = nullptr;
    cuda_check(cudaHostAlloc(&pinned, count * sizeof(Complex), cudaHostAllocDefault));
    for (int repetition = 0; repetition < 4; ++repetition)
    {
        Complex* cache = repetition % 2 == 0 ? pinned : pageable.data();
        std::fill(cache, cache + count, Complex(-9, 7));
        const int calls = model.overlap_calls;
        const bool direct_overlap = identity && repetition >= 2;
        hsolver::DiagoIterAssist<Complex, Device>::cal_hs_subspace(
            &model, psi, h, s, cache, direct_overlap);
        if (model.overlap_calls - calls != int(!direct_overlap)) { std::exit(3); }
        base_device::memory::synchronize_memory_op<Complex, CPU, Device>()(
            nullptr, nullptr, actual.data(), h, matrix_count);
        if (std::memcmp(actual.data(), host_h.data(), matrix_count * sizeof(Complex))) { std::exit(4); }
        base_device::memory::synchronize_memory_op<Complex, CPU, Device>()(
            nullptr, nullptr, actual.data(), s, matrix_count);
        if (std::memcmp(actual.data(), host_s.data(), matrix_count * sizeof(Complex))) { std::exit(5); }
        if (std::memcmp(cache, expected_action.data(), count * sizeof(Complex))) { std::exit(6); }
    }
    base_device::memory::synchronize_memory_op<Complex, CPU, Device>()(
        nullptr, nullptr, retained_input.data(), psi.get_pointer(), count);
    if (std::memcmp(retained_input.data(), input.data(), count * sizeof(Complex))) { std::exit(7); }
    cuda_check(cudaFreeHost(pinned));
    base_device::memory::delete_memory_op<Complex, Device>()(nullptr, h);
    base_device::memory::delete_memory_op<Complex, Device>()(nullptr, s);
    std::cout << "device=" << (std::is_same<Device, GPU>::value ? "gpu" : "cpu")
              << " rows=" << rows << " bands=" << bands << " identity=" << identity
              << " bit_errors=0\n";
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    GlobalV::NPOL = 1;
    GlobalV::NPROC_IN_POOL = 1;
    hsolver::createGpuBlasHandle();
    for (bool identity : {false, true})
    {
        run<CPU>(19, 23, 7, identity);
        run<GPU>(19, 23, 7, identity);
        run<GPU>(8192, 8201, 128, identity);
    }
    hsolver::destoryBLAShandle();
    MPI_Finalize();
}
