#include "module_hamilt_pw/hamilt_pwdft/kernels/onsite_op.h"
#include <cuda_runtime.h>
#include <mpi.h>
#include <vector>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <numeric>

void check(cudaError_t status)
{
    if (status != cudaSuccess)
    {
        std::cerr << cudaGetErrorString(status) << '\n';
        std::exit(1);
    }
}

class Moments
{
    int atoms_, bands_, points_, npol_, projectors_;
    int *offsets_ = nullptr, *spins_ = nullptr;
    double *moments_ = nullptr, *weights_ = nullptr, *output_ = nullptr;

  public:
    Moments(int atoms, int bands, int points, int npol, const int* counts, const int* signs)
        : atoms_(atoms), bands_(bands), points_(points), npol_(npol)
    {
        std::vector<int> offsets(atoms + 1, 0), spins(points);
        for (int atom = 0; atom < atoms; ++atom) { offsets[atom + 1] = offsets[atom] + counts[atom]; }
        for (int point = 0; point < points; ++point) { spins[point] = signs[point] > 0 ? 0 : 1; }
        projectors_ = offsets.back();
        check(cudaMalloc(&offsets_, offsets.size() * sizeof(int)));
        check(cudaMalloc(&spins_, points * sizeof(int)));
        check(cudaMalloc(&moments_, (size_t(points) * bands * atoms * 3 + 7) * sizeof(double)));
        check(cudaMalloc(&weights_, size_t(points) * bands * sizeof(double)));
        check(cudaMalloc(&output_, (atoms * 3 + 7) * sizeof(double)));
        check(cudaMemcpy(offsets_, offsets.data(), offsets.size() * sizeof(int), cudaMemcpyHostToDevice));
        check(cudaMemcpy(spins_, spins.data(), points * sizeof(int), cudaMemcpyHostToDevice));
        const double guard[] = {971, 971, 971, 971, 971, 971, 971};
        check(cudaMemcpy(output_ + atoms * 3, guard, sizeof(guard), cudaMemcpyHostToDevice));
        check(cudaMemcpy(moments_ + size_t(points) * bands * atoms * 3, guard, sizeof(guard), cudaMemcpyHostToDevice));
    }

    Moments(const Moments&) = delete;
    Moments& operator=(const Moments&) = delete;

    ~Moments()
    {
        check(cudaFree(offsets_));
        check(cudaFree(spins_));
        check(cudaFree(moments_));
        check(cudaFree(weights_));
        check(cudaFree(output_));
    }

    void capture(int point, const std::complex<double>* projections)
    {
        hamilt::projector_band_moments_gpu(bands_, npol_, atoms_, projectors_, offsets_, projections,
            moments_ + size_t(point) * bands_ * atoms_ * 3);
    }

    void finish(const double* weights, double* output)
    {
        check(cudaMemcpy(weights_, weights, size_t(points_) * bands_ * sizeof(double), cudaMemcpyHostToDevice));
        hamilt::sum_projector_moments_gpu(points_, bands_, atoms_, npol_, moments_, weights_, spins_, output_);
        check(cudaMemcpy(output, output_, atoms_ * 3 * sizeof(double), cudaMemcpyDeviceToHost));
        double guard[7];
        for (const double* end : {output_ + atoms_ * 3, moments_ + size_t(points_) * bands_ * atoms_ * 3})
        {
            check(cudaMemcpy(guard, end, sizeof(guard), cudaMemcpyDeviceToHost));
            if (!std::all_of(guard, guard + 7, [](double value) { return value == 971; }))
            { std::cerr << "Device output guard modified\n"; std::exit(5); }
        }
    }
};
using Complex = std::complex<double>;

void reference(const Complex* projections, const double* weights, const std::vector<int>& counts,
               const std::vector<int>& signs, int points, int bands, int polarizations, double* output)
{
    const int atoms = counts.size();
    const int projectors = std::accumulate(counts.begin(), counts.end(), 0);
    std::fill(output, output + atoms * 3, 0.0);
    for (int point = 0; point < points; ++point)
    {
        for (int band = 0; band < bands; ++band)
        {
            const double weight = weights[point * bands + band];
            int offset = 0;
            for (int atom = 0; atom < atoms; ++atom)
            {
                Complex occupations[4] = {};
                for (int projector = 0; projector < counts[atom]; ++projector)
                {
                    const size_t row = size_t(point * bands + band) * polarizations * projectors + offset + projector;
                    const auto up = projections[row];
                    occupations[0] += std::conj(up) * up;
                    if (polarizations == 2)
                    {
                        const auto down = projections[row + projectors];
                        occupations[1] += std::conj(up) * down;
                        occupations[2] += std::conj(down) * up;
                        occupations[3] += std::conj(down) * down;
                    }
                }
                if (polarizations == 2)
                {
                    output[3 * atom] += weight * (occupations[1] + occupations[2]).real();
                    output[3 * atom + 1] += weight * (occupations[1] - occupations[2]).imag();
                    output[3 * atom + 2] += weight * (occupations[0] - occupations[3]).real();
                }
                else { output[3 * atom + 2] += weight * occupations[0].real() * signs[point]; }
                offset += counts[atom];
            }
        }
    }
}

void run(int atoms, int bands, int points, int polarizations, bool irregular)
{
    std::vector<int> counts(atoms, 9), signs(points, 1);
    if (irregular) { for (int i = 0; i < atoms; ++i) { counts[i] = (i * 7) % 14; } }
    for (int i = 0; i < points; ++i) { signs[i] = i % 2 == 0 ? 1 : -1; }
    const int projectors = std::accumulate(counts.begin(), counts.end(), 0);
    const size_t per_point = size_t(bands) * polarizations * projectors;
    const size_t count = per_point * points, bytes = count * sizeof(Complex);
    std::vector<Complex> input(count), downloaded(count);
    std::vector<double> weights(points * bands), expected(atoms * 3), actual(atoms * 3 + 7, 971.0);
    for (size_t i = 0; i < count; ++i)
    {
        const double scale = i % 19 == 0 ? 1e-140 : 0.01;
        input[i] = scale * Complex((int(i % 509) - 257) / 509.0, (int(i % 337) - 171) / 337.0);
    }
    for (int i = 0; i < points * bands; ++i)
    { weights[i] = i % 7 == 0 ? 0.0 : i % 7 == 1 ? 1e-140 : (i % 37 - 16) / 37.0; }
    Complex* device = nullptr;
    check(cudaMalloc(&device, std::max(size_t(1), bytes)));
    if (bytes > 0) { check(cudaMemcpy(device, input.data(), bytes, cudaMemcpyHostToDevice)); }
    reference(input.data(), weights.data(), counts, signs, points, bands, polarizations, expected.data());
    auto original = [&] {
        for (int point = 0; point < points; ++point)
        {
            if (per_point > 0)
            { check(cudaMemcpy(downloaded.data() + point * per_point, device + point * per_point,
                               per_point * sizeof(Complex), cudaMemcpyDeviceToHost)); }
        }
        reference(downloaded.data(), weights.data(), counts, signs, points, bands, polarizations, actual.data());
    };
    auto candidate = [&] {
        Moments moments(atoms, bands, points, polarizations, counts.data(), signs.data());
        for (int point = 0; point < points; ++point)
        { moments.capture(point, device + point * per_point); }
        moments.finish(weights.data(), actual.data());
    };
    original();
    candidate();
    double error = 0;
    int unequal = 0;
    for (int i = 0; i < atoms * 3; ++i)
    {
        if (!std::isfinite(actual[i])) { std::exit(2); }
        error = std::max(error, std::abs(actual[i] - expected[i]));
        unequal += std::memcmp(&actual[i], &expected[i], sizeof(double)) != 0;
    }
    if (error > 2e-14 || !std::all_of(actual.begin() + atoms * 3, actual.end(), [](double x) { return x == 971; }))
    { std::cerr << "Moment or output guard mismatch: " << error << '\n'; std::exit(3); }
    if (bytes > 0) { check(cudaMemcpy(downloaded.data(), device, bytes, cudaMemcpyDeviceToHost)); }
    if (bytes > 0 && std::memcmp(input.data(), downloaded.data(), bytes) != 0)
    { std::cerr << "Projection input modified\n"; std::exit(4); }
    double milliseconds[2] = {};
    constexpr int rounds = 4;
    for (int round = 0; round < rounds; ++round)
    {
        for (int order = 0; order < 2; ++order)
        {
            const int method = (order + round) % 2;
            check(cudaDeviceSynchronize());
            const auto start = std::chrono::steady_clock::now();
            if (method == 0) { original(); } else { candidate(); }
            check(cudaDeviceSynchronize());
            milliseconds[method] += std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now() - start).count() / rounds;
        }
    }
    std::cout << "atoms=" << atoms << " bands=" << bands << " points=" << points
              << " polarizations=" << polarizations << " irregular=" << irregular
              << " unequal=" << unequal << " error=" << error << " original_ms=" << milliseconds[0]
              << " candidate_ms=" << milliseconds[1] << " speedup=" << milliseconds[0] / milliseconds[1]
              << " scratch_bytes=" << size_t(points) * bands * atoms * 3 * sizeof(double) << std::endl;
    check(cudaFree(device));
}

void check_pauli_states()
{
    const double half = 1.0 / std::sqrt(2.0);
    std::vector<Complex> input = {{1, 0}, {0, 0}, {0, 0}, {1, 0},
                                  {half, 0}, {half, 0}, {half, 0}, {0, half}};
    Complex* device = nullptr;
    check(cudaMalloc(&device, input.size() * sizeof(Complex)));
    const int counts[] = {1}, signs[] = {1};
    Moments moments(1, 4, 1, 2, counts, signs);
    const double expected[4][3] = {{0, 0, 1}, {0, 0, -1}, {1, 0, 0}, {0, 1, 0}};
    double error = 0;
    for (double phase : {0.0, 0.37})
    {
        auto rotated = input;
        for (auto& value : rotated) { value *= std::polar(1.0, phase); }
        check(cudaMemcpy(device, rotated.data(), rotated.size() * sizeof(Complex), cudaMemcpyHostToDevice));
        moments.capture(0, device);
        for (int state = 0; state < 5; ++state)
        {
            double weights[4] = {}, actual[3];
            if (state < 4) { weights[state] = 1; }
            moments.finish(weights, actual);
            for (int axis = 0; axis < 3; ++axis)
            {
                const double reference = state < 4 ? expected[state][axis] : 0;
                error = std::max(error, std::abs(actual[axis] - reference));
                if (!std::isfinite(actual[axis]) || std::abs(actual[axis] - reference) > 1e-14)
                { std::cerr << "Pauli moment or phase-invariance mismatch\n"; std::exit(6); }
            }
        }
    }
    check(cudaFree(device));
    std::cout << "Pauli states, phase invariance and repeated occupations: error=" << error << std::endl;
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    std::cout << std::setprecision(12);
    check_pauli_states();
    for (int polarizations : {1, 2})
    {
        run(1, 1, 1, polarizations, true);
        run(5, 7, 3, polarizations, true);
        run(17, 31, 2, polarizations, true);
        run(16, 308, 8, polarizations, false);
        run(32, 615, 8, polarizations, false);
    }
    MPI_Finalize();
}
