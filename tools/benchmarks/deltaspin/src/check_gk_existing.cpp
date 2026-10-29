#include "module_base/vector3.h"
#include "module_hamilt_pw/hamilt_pwdft/kernels/vnl_op.h"
#include <cstring>
#include <cuda_runtime.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

void check(cudaError_t status)
{
    if (status != cudaSuccess) { std::cerr << cudaGetErrorString(status) << '\n'; std::exit(1); }
}

int main()
{
    constexpr int nk = 8, npw = 492057, npwx = npw + 17;
    using Vector = ModuleBase::Vector3<double>;
    std::vector<Vector> g(nk * npwx), k(nk);
    for (int ik = 0; ik < nk; ++ik)
    {
        k[ik] = {0.125 * ik - 0.5, -0.03125 * ik, 0.17 * ik};
        for (int i = 0; i < npw; ++i)
        {
            g[ik * npwx + i] = {std::sin(0.013 * i), std::cos(0.017 * i), std::sin(0.11 * i + 0.3)};
        }
    }
    double *dg, *dk, *output;
    check(cudaMalloc(&dg, g.size() * sizeof(Vector)));
    check(cudaMalloc(&dk, k.size() * sizeof(Vector)));
    check(cudaMalloc(&output, (npw + 2) * sizeof(Vector)));
    check(cudaMemcpy(dg, g.data(), g.size() * sizeof(Vector), cudaMemcpyHostToDevice));
    check(cudaMemcpy(dk, k.data(), k.size() * sizeof(Vector), cudaMemcpyHostToDevice));
    auto execute = [&](int method, int ik)
    {
        if (method == 0)
        {
            std::unique_ptr<Vector[]> host(new Vector[npw]);
            #pragma omp parallel for schedule(static, 512)
            for (int i = 0; i < npw; ++i) { host[i] = g[ik * npwx + i] + k[ik]; }
            check(cudaMemcpy(output, host.get(), npw * sizeof(Vector), cudaMemcpyHostToDevice));
        }
        else
        {
            hamilt::cal_gk_gpu(npw, dg + ik * npwx * 3, dk + ik * 3, output);
            check(cudaGetLastError());
        }
    };
    std::vector<Vector> reference(npw), actual(npw);
    double error = 0;
    for (int ik = 0; ik < nk; ++ik)
    {
        execute(0, ik);
        check(cudaMemcpy(reference.data(), output, npw * sizeof(Vector), cudaMemcpyDeviceToHost));
        execute(1, ik);
        check(cudaMemcpy(actual.data(), output, npw * sizeof(Vector), cudaMemcpyDeviceToHost));
        for (int i = 0; i < npw; ++i)
        {
            error = std::max(error, (actual[i] - reference[i]).norm());
        }
    }
    if (error != 0) { std::cerr << "G+k mismatch " << error << '\n'; return 2; }
    for (int count : {0, 1, 85, 86, 257, npw})
    {
        std::vector<double> guarded(3 * (npw + 2), -731.25);
        check(cudaMemcpy(output, guarded.data(), guarded.size() * sizeof(double), cudaMemcpyHostToDevice));
        hamilt::cal_gk_gpu(count, dg + 3 * npwx, dk + 3, output + 3);
        check(cudaGetLastError());
        check(cudaMemcpy(guarded.data(), output, guarded.size() * sizeof(double), cudaMemcpyDeviceToHost));
        for (int i = 0; i < 3; ++i)
        {
            if (guarded[i] != -731.25) { return 3; }
        }
        for (int i = 0; i < count; ++i)
        {
            Vector expected = g[npwx + i] + k[1];
            if (std::memcmp(guarded.data() + 3 * (i + 1), &expected, sizeof(Vector)) != 0) { return 4; }
        }
        for (int i = 3 * (count + 1); i < static_cast<int>(guarded.size()); ++i)
        {
            if (guarded[i] != -731.25) { return 5; }
        }
    }
    std::vector<Vector> copied_g(g.size()), copied_k(k.size());
    check(cudaMemcpy(copied_g.data(), dg, g.size() * sizeof(Vector), cudaMemcpyDeviceToHost));
    check(cudaMemcpy(copied_k.data(), dk, k.size() * sizeof(Vector), cudaMemcpyDeviceToHost));
    if (std::memcmp(g.data(), copied_g.data(), g.size() * sizeof(Vector)) != 0
        || std::memcmp(k.data(), copied_k.data(), k.size() * sizeof(Vector)) != 0) { return 6; }
    double times[2][5] = {};
    for (int repeat = 0; repeat < 5; ++repeat)
    {
        for (int order = 0; order < 2; ++order)
        {
            const int method = (order + repeat) % 2;
            for (int i = 0; i < 8; ++i) { execute(method, i); }
            check(cudaDeviceSynchronize());
            const auto start = std::chrono::steady_clock::now();
            for (int i = 0; i < 40; ++i) { execute(method, i % nk); }
            check(cudaDeviceSynchronize());
            times[method][repeat] = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count() / 40;
        }
    }
    for (auto& row : times) { std::sort(row, row + 5); }
    std::cout << "npw=" << npw << " max_error=" << error
              << " baseline_ms=" << times[0][2] * 1000 << " device_ms=" << times[1][2] * 1000
              << " speedup=" << times[0][2] / times[1][2] << std::endl;
    check(cudaFree(dg)); check(cudaFree(dk)); check(cudaFree(output));
}
