#include "module_hamilt_pw/hamilt_pwdft/kernels/force_op.h"
#include <base/macros/macros.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <string>
#include <vector>

using Complex = std::complex<double>;

int check(int projectors, int count, int direction)
{
    const int rows = count + 7, stride = rows + 3;
    const size_t elements = size_t(projectors) * rows;
    std::vector<int> indices(3 * stride, -1);
    indices[direction * stride] = count;
    for (int i = 0; i < count; ++i) { indices[direction * stride + 1 + i] = rows - 1 - i; }
    std::vector<Complex> values(elements);
    for (size_t i = 0; i < elements; ++i) { values[i] = Complex(0.125 * (i % 101), -0.25 * (i % 59)); }
    int* device_indices = nullptr;
    Complex *device_values = nullptr, *device_saved = nullptr;
    cudaErrcheck(cudaMalloc(reinterpret_cast<void**>(&device_indices), indices.size() * sizeof(int)));
    cudaErrcheck(cudaMalloc(reinterpret_cast<void**>(&device_values), std::max(size_t(1), elements) * sizeof(Complex)));
    cudaErrcheck(cudaMalloc(reinterpret_cast<void**>(&device_saved), std::max(size_t(1), size_t(projectors) * count) * sizeof(Complex)));
    cudaErrcheck(cudaMemcpy(device_indices, indices.data(), indices.size() * sizeof(int), cudaMemcpyHostToDevice));
    cudaErrcheck(cudaMemcpy(device_values, values.data(), elements * sizeof(Complex), cudaMemcpyHostToDevice));
    hamilt::saveVkbValues<double>(device_indices, device_values, device_saved, projectors, count,
                                 rows, direction, stride);
    const auto saved_status = cudaDeviceSynchronize();
    if (saved_status != cudaSuccess)
    {
        std::cerr << "SAVE_FAILED projectors=" << projectors << " count=" << count
                  << " status=" << int(saved_status) << '\n';
        return 2;
    }
    std::vector<Complex> saved(size_t(projectors) * count);
    cudaErrcheck(cudaMemcpy(saved.data(), device_saved, saved.size() * sizeof(Complex), cudaMemcpyDeviceToHost));
    for (int projector = 0; projector < projectors; ++projector)
    for (int i = 0; i < count; ++i)
    {
        if (saved[size_t(projector) * count + i] != values[size_t(projector) * rows + rows - 1 - i]) { return 3; }
    }
    const Complex sentinel(23, -17), coefficient(-0.4, 0.7);
    std::vector<Complex> actual(elements, sentinel);
    cudaErrcheck(cudaMemcpy(device_values, actual.data(), elements * sizeof(Complex), cudaMemcpyHostToDevice));
    hamilt::revertVkbValues<double>(device_indices, device_values, device_saved, projectors, count,
                                   rows, direction, stride, coefficient);
    const auto restored_status = cudaDeviceSynchronize();
    if (restored_status != cudaSuccess)
    {
        std::cerr << "RESTORE_FAILED status=" << int(restored_status) << '\n';
        return 4;
    }
    cudaErrcheck(cudaMemcpy(actual.data(), device_values, elements * sizeof(Complex), cudaMemcpyDeviceToHost));
    double error = 0;
    for (size_t i = 0; i < elements; ++i)
    {
        const Complex expected = int(i % rows) >= rows - count ? values[i] * coefficient : sentinel;
        const double difference = std::abs(actual[i] - expected);
        if (!std::isfinite(difference)) { return 5; }
        error = std::max(error, difference);
    }
    cudaErrcheck(cudaFree(device_indices));
    cudaErrcheck(cudaFree(device_values));
    cudaErrcheck(cudaFree(device_saved));
    std::cout << "PROJECTOR_STORAGE_OK projectors=" << projectors << " count=" << count
              << " direction=" << direction << " error=" << error << '\n';
    return error < 1e-12 ? 0 : 6;
}

int main(int argc, char** argv)
{
    if (argc == 2 && std::string(argv[1]) == "--large") { return check(2049, 4096, 1); }
    for (int direction = 0; direction < 3; ++direction)
    {
        for (int count : {0, 1, 9, 257})
        {
            const int status = check(7, count, direction);
            if (status != 0) { return status; }
        }
    }
    return 0;
}
