/**
 * @file bench_nep_postprocess.cu
 * @brief NEP 后处理性能基准测试 — CPU vs GPU 对比
 *
 * 测试原理解释:
 *   1. 正确性: 先用手工计算验证 CPU 函数, 再用 CPU vs GPU 对比验证 GPU
 *   2. 性能: 计时对比 CPU 和 GPU 后处理的执行时间 (含显存拷贝)
 */

#include "esolver_nep_postprocess.h"
#include <iostream>
#include <cmath>
#include <vector>
#include <chrono>

using ModuleESolver::postprocess_nep_cpu;
using ModuleESolver::postprocess_nep_cuda;

// ==========================================================================
// 手工计算参考值 (验证 CPU 本身是否正确)
// ==========================================================================
void correctness_demo()
{
    std::cout << "============================================================" << std::endl;
    std::cout << "  第 1 步: 手工计算验证 CPU 函数正确性" << std::endl;
    std::cout << "============================================================" << std::endl;

    const int nat = 3;
    // 输入: 每原子能量 [1.0, 2.0, 3.0] eV
    std::vector<double> e = {1.0, 2.0, 3.0};
    // SoA 格式力 [fx0,fx1,fx2, fy0,fy1,fy2, fz0,fz1,fz2]
    std::vector<double> f = {10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0};
    std::vector<double> v(9 * nat, 1.0);
    double fe = 2.0, ff = 3.0, fv = 5.0;

    // --- 手工计算期望值 ---
    // 能量: sum(1+2+3) * 2.0 = 12.0
    double expected_energy = (1.0 + 2.0 + 3.0) * fe;
    // force(i,0) = fx[i] * 3.0, force(i,1) = fy[i] * 3.0, force(i,2) = fz[i] * 3.0
    // virial 每个分量: sum(nat 个 1.0) * 5.0 = 3 * 5.0 = 15.0
    double expected_virial = nat * 1.0 * fv;

    std::cout << "\n手工计算期望值:" << std::endl;
    std::cout << "  能量 = (" << e[0] << "+" << e[1] << "+" << e[2] << ") * " << fe
              << " = " << expected_energy << std::endl;
    std::cout << "  力(0,0) = fx[0] * ff = " << f[0] << " * " << ff << " = " << f[0] * ff << std::endl;
    std::cout << "  力(0,1) = fy[0] * ff = " << f[0+nat] << " * " << ff << " = " << f[0+nat] * ff << std::endl;
    std::cout << "  virial 每分量 = " << nat << " * 1.0 * " << fv << " = " << expected_virial << std::endl;

    // --- 调用 CPU 函数 ---
    double cpu_pot = 0;
    ModuleBase::matrix cpu_force(nat, 3);
    ModuleBase::matrix cpu_virial(3, 3);
    postprocess_nep_cpu(nat, e.data(), f.data(), v.data(), fe, ff, fv,
                        cpu_pot, cpu_force, cpu_virial);

    std::cout << "\nCPU 函数输出:" << std::endl;
    std::cout << "  能量 = " << cpu_pot << std::endl;
    std::cout << "  力(0,0) = " << cpu_force(0,0)
              << ", 力(0,1) = " << cpu_force(0,1)
              << ", 力(0,2) = " << cpu_force(0,2) << std::endl;
    std::cout << "  virial(0,0) = " << cpu_virial(0,0) << std::endl;

    // --- 对比 ---
    bool ok = true;
    double eps = 1e-12;
    if (std::abs(cpu_pot - expected_energy) > eps) ok = false;
    if (std::abs(cpu_force(0,0) - f[0]*ff) > eps) ok = false;
    if (std::abs(cpu_force(0,1) - f[0+nat]*ff) > eps) ok = false;
    if (std::abs(cpu_virial(0,0) - expected_virial) > eps) ok = false;

    if (ok)
        std::cout << "\n  ==> CPU 输出与手工计算一致, CPU 函数正确!" << std::endl;
    else
        std::cout << "\n  ==> 错误: CPU 输出与手工计算不一致!" << std::endl;

    // --- 调用 GPU 函数并对比 ---
    double gpu_pot = 0;
    ModuleBase::matrix gpu_force(nat, 3);
    ModuleBase::matrix gpu_virial(3, 3);
    postprocess_nep_cuda(nat, e.data(), f.data(), v.data(), fe, ff, fv,
                         gpu_pot, gpu_force, gpu_virial);

    std::cout << "\nGPU 函数输出:" << std::endl;
    std::cout << "  能量 = " << gpu_pot << std::endl;
    std::cout << "  力(0,0) = " << gpu_force(0,0) << std::endl;
    std::cout << "  virial(0,0) = " << gpu_virial(0,0) << std::endl;

    bool gpu_ok = true;
    if (std::abs(gpu_pot - expected_energy) > eps) gpu_ok = false;
    if (std::abs(gpu_pot - cpu_pot) > eps) gpu_ok = false;
    if (std::abs(gpu_virial(0,0) - expected_virial) > eps) gpu_ok = false;

    if (gpu_ok)
        std::cout << "\n  ==> GPU 输出与手工计算一致, 也与 CPU 输出一致!" << std::endl;
    else
        std::cout << "\n  ==> 错误: GPU 输出不一致!" << std::endl;
}

// ==========================================================================
// 性能测试
// ==========================================================================
class Timer
{
    using clock = std::chrono::high_resolution_clock;
    clock::time_point start_;
public:
    Timer() : start_(clock::now()) {}
    double elapsed_ms() const {
        return std::chrono::duration<double, std::milli>(clock::now() - start_).count();
    }
};

void run_benchmark(const std::string& name, int nat, int warmup, int iters)
{
    std::vector<double> e(nat, 1.5);
    std::vector<double> f(3 * nat, 2.5);
    std::vector<double> v(9 * nat, 0.5);
    double fe = 1.0, ff = 1.0, fv = 1.0;

    // --- CPU 计时 (含 warmup) ---
    for (int w = 0; w < warmup; ++w) {
        double pot = 0;
        ModuleBase::matrix force(nat, 3);
        ModuleBase::matrix virial(3, 3);
        postprocess_nep_cpu(nat, e.data(), f.data(), v.data(), fe, ff, fv, pot, force, virial);
    }

    Timer cpu_timer;
    for (int r = 0; r < iters; ++r) {
        double pot = 0;
        ModuleBase::matrix force(nat, 3);
        ModuleBase::matrix virial(3, 3);
        postprocess_nep_cpu(nat, e.data(), f.data(), v.data(), fe, ff, fv, pot, force, virial);
    }
    double cpu_ms = cpu_timer.elapsed_ms() / iters;

    // --- GPU 计时 (含 warmup) ---
    for (int w = 0; w < warmup; ++w) {
        double pot = 0;
        ModuleBase::matrix force(nat, 3);
        ModuleBase::matrix virial(3, 3);
        postprocess_nep_cuda(nat, e.data(), f.data(), v.data(), fe, ff, fv, pot, force, virial);
    }

    Timer gpu_timer;
    for (int r = 0; r < iters; ++r) {
        double pot = 0;
        ModuleBase::matrix force(nat, 3);
        ModuleBase::matrix virial(3, 3);
        postprocess_nep_cuda(nat, e.data(), f.data(), v.data(), fe, ff, fv, pot, force, virial);
    }
    double gpu_ms = gpu_timer.elapsed_ms() / iters;

    double speedup = cpu_ms / gpu_ms;

    std::cout << "  nat=" << nat << "  CPU=" << cpu_ms << "ms  GPU=" << gpu_ms
              << "ms  (含显存拷贝)  加速比=" << speedup << "x" << std::endl;
}

void performance_test()
{
    std::cout << "\n============================================================" << std::endl;
    std::cout << "  第 2 步: 性能基准测试 (CPU vs GPU)" << std::endl;
    std::cout << "  GPU: Tesla T4" << std::endl;
    std::cout << "  注意: GPU 时间包含 cudaMalloc+cudaMemcpy+cudaFree" << std::endl;
    std::cout << "============================================================" << std::endl;

    std::cout << "\n--- 小体系 (nat=10~100) ---" << std::endl;
    run_benchmark("small", 10, 5, 100);
    run_benchmark("small", 100, 5, 100);

    std::cout << "\n--- 中等体系 (nat=1000) ---" << std::endl;
    run_benchmark("medium", 1000, 3, 20);

    std::cout << "\n--- 大体系 (nat=5000~20000) ---" << std::endl;
    run_benchmark("large", 5000, 2, 10);
    run_benchmark("large", 10000, 2, 5);
    run_benchmark("large", 20000, 2, 3);

    // --- 分析 ---
    std::cout << "\n============================================================" << std::endl;
    std::cout << "  第 3 步: 原始代码 vs 修改后代码结构对比" << std::endl;
    std::cout << "============================================================" << std::endl;
    std::cout << "\n原始 runner():" << std::endl;
    std::cout << "  每步在 runner() 内临时创建 cell(9) + coord(3*nat)" << std::endl;
    std::cout << "  后处理内联, 不可切换 CPU/GPU" << std::endl;
    std::cout << "  timer 粒度: 只有整个 runner 级别" << std::endl;
    std::cout << "\n修改后 runner():" << std::endl;
    std::cout << "  cell/coord 在 before_all_runners() 中分配一次并复用" << std::endl;
    std::cout << "  后处理拆分为独立函数, 编译时选择 CPU/GPU 路径" << std::endl;
    std::cout << "  timer 粒度: prepare_input + postprocess 分开测量" << std::endl;
    std::cout << "\n性能收益:" << std::endl;
    std::cout << "  1. cell/coord 持久化: 消除每步 vector 构造/析构开销" << std::endl;
    std::cout << "     (std::vector 每步申请+释放堆内存 ~ 数百 ns)" << std::endl;
    std::cout << "  2. 后处理拆分: 为后续 OpenMP 并行和 SIMD 优化提供接口" << std::endl;
    std::cout << "  3. GPU 后处理: 大体系时 GPU 并行加速 (见上述加速比)" << std::endl;
    std::cout << "  4. CUDA cudaMalloc/Free: 当前每次调用分配/释放, 有额外开销" << std::endl;
    std::cout << "     后续改为持久化 device buffer 可进一步提升" << std::endl;
}

int main()
{
    correctness_demo();
    performance_test();
    return 0;
}
