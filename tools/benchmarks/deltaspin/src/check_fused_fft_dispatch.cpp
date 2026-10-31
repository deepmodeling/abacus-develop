#include "module_basis/module_pw/kernels/fused_veff.h"

#include <cuda_runtime.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace
{
using Complex = std::complex<double>;
using Reference = std::complex<long double>;

void require(bool condition, const char* message)
{
    if (!condition) { std::cerr << message << '\n'; std::exit(2); }
}

void check(cudaError_t status)
{
    if (status != cudaSuccess) { std::cerr << cudaGetErrorString(status) << '\n'; std::exit(1); }
}

bool smooth(int length)
{
    for (int factor : {2, 3, 5})
    {
        while (length % factor == 0) { length /= factor; }
    }
    return length == 1;
}

double check_transform(int length, const ModulePW::FusedFFT& fft)
{
    constexpr int sticks = 3, bands = 3;
    const int rows = sticks * length;
    std::vector<int> reverse(rows, -1);
    int npw = 0;
    for (int i = 0; i < rows; ++i)
    {
        if (i % 4 != 1) { reverse[i] = npw++; }
    }
    const int stride = npw + 7;
    std::vector<Complex> input(2 * bands * stride), initial(input.size()), actual(input.size());
    for (size_t i = 0; i < input.size(); ++i)
    {
        input[i] = Complex(std::sin(0.27 * i), std::cos(0.19 * i));
        initial[i] = Complex(std::cos(0.13 * i), std::sin(0.11 * i));
    }
    const Complex guard(17.0, -23.0);
    std::vector<Complex> fields(2 * bands * rows + 2, guard);
    Complex *device_input = nullptr, *device_output = nullptr, *device_fields = nullptr;
    int* device_reverse = nullptr;
    check(cudaMalloc(&device_input, input.size() * sizeof(Complex)));
    check(cudaMalloc(&device_output, input.size() * sizeof(Complex)));
    check(cudaMalloc(&device_fields, fields.size() * sizeof(Complex)));
    check(cudaMalloc(&device_reverse, reverse.size() * sizeof(int)));
    check(cudaMemcpy(device_input, input.data(), input.size() * sizeof(Complex), cudaMemcpyHostToDevice));
    check(cudaMemcpy(device_output, initial.data(), input.size() * sizeof(Complex), cudaMemcpyHostToDevice));
    check(cudaMemcpy(device_fields, fields.data(), fields.size() * sizeof(Complex), cudaMemcpyHostToDevice));
    check(cudaMemcpy(device_reverse, reverse.data(), reverse.size() * sizeof(int), cudaMemcpyHostToDevice));
    fft.load_sticks(device_input, device_fields + 1, device_reverse, sticks, bands, stride);
    check(cudaMemcpy(fields.data(), device_fields, fields.size() * sizeof(Complex), cudaMemcpyDeviceToHost));
    require(fields.front() == guard && fields.back() == guard, "Sparse inverse FFT overwrote a guard");
    std::vector<Reference> roots(length);
    for (int i = 0; i < length; ++i)
    {
        roots[i] = std::polar(1.0L, 2.0L * std::acos(-1.0L) * i / length);
    }
    double error = 0;
    for (int spin = 0; spin < 2; ++spin)
    for (int band = 0; band < bands; ++band)
    for (int stick = 0; stick < sticks; ++stick)
    for (int z = 0; z < length; ++z)
    {
        Reference expected = 0;
        long double norm = 0;
        for (int k = 0; k < length; ++k)
        {
            const int ig = reverse[stick * length + k];
            if (ig < 0) { continue; }
            const Reference value = input[(2 * band + spin) * stride + ig];
            expected += value * roots[(z * k) % length];
            norm += std::abs(value);
        }
        const Reference value = fields[1 + (spin * bands + band) * rows + stick * length + z];
        require(std::isfinite(value.real()) && std::isfinite(value.imag()), "Non-finite inverse FFT output");
        error = std::max(error, double(std::abs(value - expected) / std::max(1.0L, norm)));
    }
    fft.add_sticks(device_fields + 1, device_output, device_reverse, sticks, bands, stride, length);
    check(cudaMemcpy(actual.data(), device_output, actual.size() * sizeof(Complex), cudaMemcpyDeviceToHost));
    for (size_t i = 0; i < actual.size(); ++i)
    {
        const Complex expected = initial[i] + (i % stride < size_t(npw) ? input[i] : Complex(0));
        require(std::isfinite(actual[i].real()) && std::isfinite(actual[i].imag()), "Non-finite roundtrip output");
        error = std::max(error, std::abs(actual[i] - expected));
    }
    check(cudaMemcpy(actual.data(), device_input, actual.size() * sizeof(Complex), cudaMemcpyDeviceToHost));
    require(actual == input, "Sparse FFT modified its input");
    check(cudaFree(device_reverse));
    check(cudaFree(device_fields));
    check(cudaFree(device_output));
    check(cudaFree(device_input));
    require(error < 2e-12, "Sparse FFT disagrees with independent DFT or normalized roundtrip");
    return error;
}
} // namespace

int main(int argc, char** argv)
{
    const bool unavailable = argc == 2 && std::string(argv[1]) == "--expect-unavailable";
    require(argc == 1 || unavailable, "Expected no arguments or --expect-unavailable");
    int devices = 0;
    check(cudaGetDeviceCount(&devices));
    require(devices > 0, "No CUDA device available");
    // Revisit the first device to exercise invalidation after a device switch.
    for (int visit = 0; visit <= std::min(devices, 2); ++visit)
    {
        const int device = visit % std::min(devices, 2);
        check(cudaSetDevice(device));
        require(ModulePW::fused_veff_supported() != unavailable, "Unexpected architecture eligibility");
        double error = 0;
        int count = 0;
        for (int length = 16; length <= 256; ++length)
        {
            const auto* fft = ModulePW::fused_fft(length);
            const bool power_of_two = (length & (length - 1)) == 0;
            const bool valid = unavailable ? fft == nullptr
                                           : (smooth(length) || fft == nullptr) && (!power_of_two || fft != nullptr);
            if (!valid)
            {
                std::cerr << "device=" << device << " length=" << length
                          << " available=" << (fft != nullptr) << '\n';
            }
            require(valid, "Unexpected transform registry entry");
            require(ModulePW::fused_fft(length) == fft, "Unstable cached transform lookup");
            if (fft == nullptr) { continue; }
            if (visit < std::min(devices, 2)) { error = std::max(error, check_transform(length, *fft)); }
            ++count;
        }
        for (int length : {-1, 0, 15, 257})
        {
            require(ModulePW::fused_fft(length) == nullptr, "Unsupported length returned kernels");
        }
        std::cout << "device=" << device << " visit=" << visit << " lengths=" << count
                  << " error=" << error << " checks=PASS\n";
    }
}
