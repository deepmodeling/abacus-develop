#include "module_hamilt_lcao/module_deltaspin/hpsi_transfer.h"

#include <algorithm>
#include <cstring>
#include <iostream>
#include <vector>

using Complex = std::complex<double>;

void run(int points, int source_leading, int leading, int bands)
{
    const size_t count = size_t(points) * bands * source_leading;
    Complex* host = nullptr;
    Complex* output = nullptr;
    const std::vector<int> indices = points == 1 ? std::vector<int>{0, 0, 0}
        : std::vector<int>{0, 1, 2, 3, 4, 5, 6, 7, 7, 1, 3, 2, 0, 7, 2};
    const size_t stride = size_t(bands) * leading + 2;
    std::vector<Complex> actual(stride * indices.size(), Complex(999, -777));
    auto expected = actual;
    cudaErrcheck(cudaHostAlloc(&host, count * sizeof(Complex), cudaHostAllocDefault));
    cudaErrcheck(cudaMalloc(&output, actual.size() * sizeof(Complex)));
    cudaErrcheck(cudaMemcpy(output, actual.data(), actual.size() * sizeof(Complex), cudaMemcpyHostToDevice));
    for (size_t i = 0; i < count; ++i) { host[i] = Complex(double(i % 1789) / 2048, -double(i % 7919) / 16384); }
    const std::vector<Complex> original(host, host + count);
    {
        spinconstrain::HPsiTransfer unused(host, points, bands, source_leading);
    }
    {
        spinconstrain::HPsiTransfer transfer(host, points, bands, source_leading);
        for (size_t access = 0; access < indices.size(); ++access)
        {
            const int rows = leading - access % 3;
            const int ik = indices[access];
            const size_t offset = access * stride + 1;
            transfer.restore(ik, output + offset, rows);
            for (int band = 0; band < bands; ++band)
            {
                const Complex* source = host + (size_t(ik) * bands + band) * source_leading;
                std::copy(source, source + rows, expected.data() + offset + size_t(band) * rows);
            }
        }
        // Destruction drains the final restore and an unconsumed prefetch.
        // No host synchronization separates consecutive requests above.
    }
    cudaErrcheck(cudaMemcpy(actual.data(), output, actual.size() * sizeof(Complex), cudaMemcpyDeviceToHost));
    if (std::memcmp(actual.data(), expected.data(), actual.size() * sizeof(Complex)) != 0)
    {
        std::cerr << "Transfer output or guard changed: points=" << points << " source_leading=" << source_leading << '\n';
        std::exit(2);
    }
    if (std::memcmp(original.data(), host, count * sizeof(Complex)) != 0) { std::exit(3); }
    cudaErrcheck(cudaFree(output));
    cudaErrcheck(cudaFreeHost(host));
    std::cout << "points=" << points << " source_leading=" << source_leading << " leading=" << leading
              << " bands=" << bands << " requests=" << indices.size()
              << " bit_errors=0 guards_unchanged=true input_unchanged=true\n";
}

int main()
{
    for (int lifecycle = 0; lifecycle < 2; ++lifecycle)
    {
        run(8, 17, 11, 7);
        run(8, 513, 513, 19);
        run(8, 4099, 4093, 73);
        run(1, 17, 11, 7);
    }
}
