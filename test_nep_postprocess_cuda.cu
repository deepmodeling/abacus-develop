/**
 * @file test_nep_postprocess_cuda.cu
 * @brief GPU 单元测试: 对比 postprocess_nep_cpu 与 postprocess_nep_cuda 的输出一致性
 *
 * 测试方法:
 *   对每个测试用例, 分别调用 CPU 和 GPU 后处理函数,
 *   逐项对比两者的能量、力、virial 输出, 验证 GPU 实现与 CPU 参考一致。
 *
 * 编译:
 *   nvcc -std=c++11 -D__CUDA -I source -I source/source_base -I source/source_esolver \
 *        test_nep_postprocess_cuda.cu \
 *        source/source_esolver/esolver_nep_postprocess.cpp \
 *        source/source_esolver/esolver_nep_postprocess.cu \
 *        source/source_base/matrix.cpp \
 *        source/source_base/module_external/blas_connector_base.cpp \
 *        source/source_base/module_external/blas_connector_vector.cpp \
 *        source/source_base/module_external/blas_connector_matrix.cpp \
 *        -L/usr/lib/x86_64-linux-gnu -lblas \
 *        -o test_nep_postprocess_cuda
 *
 * 运行:
 *   ./test_nep_postprocess_cuda
 */

#include "esolver_nep_postprocess.h"
#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <cstring>

using ModuleESolver::postprocess_nep_cpu;
using ModuleESolver::postprocess_nep_cuda;

int tests_passed = 0;
int tests_failed = 0;

const double epsilon = 1e-10;

void assert_double_eq(const std::string& label, double cpu_val, double gpu_val)
{
    if (std::abs(cpu_val - gpu_val) < epsilon)
    {
        std::cout << "  [PASS] " << label << ": CPU=" << cpu_val << " GPU=" << gpu_val << std::endl;
        tests_passed++;
    }
    else
    {
        std::cout << "  [FAIL] " << label << ": CPU=" << cpu_val << " GPU=" << gpu_val
                  << " (diff=" << std::abs(cpu_val - gpu_val) << ")" << std::endl;
        tests_failed++;
    }
}

/**
 * @brief 单次 CPU vs GPU 对比
 *
 * @param nat 原子数
 * @param atomic_energy 每原子能量 (nat 个)
 * @param raw_force 力 (3*nat 个, SoA 布局)
 * @param raw_virial virial (9*nat 个, SoA 布局)
 * @param fact_e/ fact_f/ fact_v 换算因子
 * @param test_name 测试名称
 */
void compare_cpu_gpu(int nat,
                     const std::vector<double>& atomic_energy,
                     const std::vector<double>& raw_force,
                     const std::vector<double>& raw_virial,
                     double fact_e, double fact_f, double fact_v,
                     const std::string& test_name)
{
    std::cout << "\n=== " << test_name << " (nat=" << nat << ") ===" << std::endl;

    // === CPU 路径 ===
    double cpu_potential = 0.0;
    ModuleBase::matrix cpu_force(nat, 3);
    ModuleBase::matrix cpu_virial(3, 3);

    postprocess_nep_cpu(nat,
                        atomic_energy.data(),
                        raw_force.data(),
                        raw_virial.data(),
                        fact_e, fact_f, fact_v,
                        cpu_potential,
                        cpu_force,
                        cpu_virial);

    // === GPU 路径 ===
    double gpu_potential = 0.0;
    ModuleBase::matrix gpu_force(nat, 3);
    ModuleBase::matrix gpu_virial(3, 3);

    postprocess_nep_cuda(nat,
                         atomic_energy.data(),
                         raw_force.data(),
                         raw_virial.data(),
                         fact_e, fact_f, fact_v,
                         gpu_potential,
                         gpu_force,
                         gpu_virial);

    // === 对比 ===
    assert_double_eq("energy", cpu_potential, gpu_potential);

    // 对比力 (抽样对比前 10 个和后 10 个, 避免输出过长)
    for (int i = 0; i < nat; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            std::string label = "force(" + std::to_string(i) + "," + std::to_string(j) + ")";
            if (std::abs(cpu_force(i, j) - gpu_force(i, j)) >= epsilon)
            {
                assert_double_eq(label, cpu_force(i, j), gpu_force(i, j));
            }
            else
            {
                tests_passed++; // 快速通过
            }
        }
    }

    // 对比 virial (9 个分量)
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            std::string label = "virial(" + std::to_string(i) + "," + std::to_string(j) + ")";
            assert_double_eq(label, cpu_virial(i, j), gpu_virial(i, j));
        }
    }
}

// ============================================================================
// Test 1: 单原子 — 基础 GPU 正确性
// ============================================================================
void test_single_atom()
{
    const int nat = 1;
    std::vector<double> e = {2.0};
    std::vector<double> f = {3.0, 4.0, 5.0};
    std::vector<double> v = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
    double fe = 1.0, ff = 1.0, fv = 1.0;

    compare_cpu_gpu(nat, e, f, v, fe, ff, fv, "单原子基础测试");
}

// ============================================================================
// Test 2: 多原子 (nat=4) — SoA 布局验证
// ============================================================================
void test_multi_atom()
{
    const int nat = 4;
    std::vector<double> e = {1.0, 2.0, 3.0, 4.0};
    // SoA 布局: [fx0,fx1,fx2,fx3, fy0,fy1,fy2,fy3, fz0,fz1,fz2,fz3]
    std::vector<double> f = {
        1.0, 2.0, 3.0, 4.0,      // fx
        5.0, 6.0, 7.0, 8.0,      // fy
        9.0, 10.0, 11.0, 12.0    // fz
    };
    std::vector<double> v(9 * nat, 1.0);
    double fe = 2.0, ff = 0.5, fv = 3.0;

    compare_cpu_gpu(nat, e, f, v, fe, ff, fv, "多原子 SoA 测试");
}

// ============================================================================
// Test 3: 中等大小 (nat=100) — 数值稳定性
// ============================================================================
void test_medium()
{
    const int nat = 100;
    std::vector<double> e(nat, 1.0);
    std::vector<double> f(3 * nat, 0.5);
    std::vector<double> v(9 * nat, 0.1);
    double fe = 2.0, ff = 3.0, fv = 4.0;

    compare_cpu_gpu(nat, e, f, v, fe, ff, fv, "中等体系 (nat=100)");
}

// ============================================================================
// Test 4: 大体系 (nat=5000) — GPU 并行压力测试
// ============================================================================
void test_large()
{
    const int nat = 5000;
    std::vector<double> e(nat);
    std::vector<double> f(3 * nat);
    std::vector<double> v(9 * nat);

    // 用随机模式填充 (避免全相同值掩盖 bug)
    for (int i = 0; i < nat; ++i)
        e[i] = (i % 17) * 0.1 + 1.0;
    for (int i = 0; i < 3 * nat; ++i)
        f[i] = (i % 23) * 0.05 - 0.5;
    for (int i = 0; i < 9 * nat; ++i)
        v[i] = (i % 31) * 0.02 - 0.3;

    double fe = 1.5, ff = 0.8, fv = 2.0;

    compare_cpu_gpu(nat, e, f, v, fe, ff, fv, "大体系 (nat=5000)");
}

// ============================================================================
// Test 5: 非均匀换算因子
// ============================================================================
void test_unit_conversion()
{
    const int nat = 10;
    // 使用真实物理单位换算因子
    std::vector<double> e(nat, 1.5);
    std::vector<double> f(3 * nat, 2.5);
    std::vector<double> v(9 * nat, 0.5);

    // 模拟 esolver_nep.cpp 中的实际换算因子
    const double Ry_to_eV = 13.605703976;
    const double ANGSTROM_AU = 1.88972612546;
    double fact_e = 1.0 / Ry_to_eV;
    double fact_f = 1.0 / (Ry_to_eV * ANGSTROM_AU);
    double fact_v = 1.0 / (100.0 * Ry_to_eV); // 假设 omega=100

    compare_cpu_gpu(nat, e, f, v, fact_e, fact_f, fact_v, "真实物理单位换算");
}

// ============================================================================
// Test 6: atomicAdd 原子操作正确性 (重复运行验证确定性)
// ============================================================================
void test_atomic_add_consistency()
{
    std::cout << "\n=== 原子操作一致性 (nat=2000, 运行 3 次) ===" << std::endl;
    const int nat = 2000;
    std::vector<double> e(nat, 0.1);
    std::vector<double> f(3 * nat, 0.2);
    std::vector<double> v(9 * nat, 0.05);

    double fe = 1.0, ff = 1.0, fv = 1.0;

    // 运行 3 次 GPU, 验证每次结果一致 (atomicAdd 确定性)
    double prev_potential = -1.0;
    std::vector<double> prev_virial(9, -1.0);

    for (int run = 0; run < 3; ++run)
    {
        double potential = 0.0;
        ModuleBase::matrix force(nat, 3);
        ModuleBase::matrix virial(3, 3);

        postprocess_nep_cuda(nat, e.data(), f.data(), v.data(),
                             fe, ff, fv, potential, force, virial);

        if (run > 0)
        {
            assert_double_eq("run" + std::to_string(run) + " energy consistent",
                             prev_potential, potential);
        }
        prev_potential = potential;

        std::cout << "  Run " << run << ": potential=" << potential
                  << ", virial(0,0)=" << virial(0, 0) << std::endl;
    }
}

// ============================================================================
// main
// ============================================================================
int main()
{
    std::cout << "============================================================" << std::endl;
    std::cout << "  NEP CUDA Postprocess Test — CPU vs GPU 对比验证" << std::endl;
    std::cout << "  GPU: Tesla T4, CUDA Driver 12.2, nvcc 11.5" << std::endl;
    std::cout << "============================================================" << std::endl;

    test_single_atom();
    test_multi_atom();
    test_medium();
    test_large();
    test_unit_conversion();
    test_atomic_add_consistency();

    std::cout << "\n============================================================" << std::endl;
    std::cout << "  Results: " << tests_passed << " passed, "
              << tests_failed << " failed" << std::endl;
    std::cout << "============================================================" << std::endl;

    if (tests_failed > 0)
    {
        std::cerr << "\n[FAIL] 存在 " << tests_failed << " 项 GPU 对比测试失败!" << std::endl;
        return 1;
    }

    std::cout << "\n[PASS] CPU 与 GPU 输出完全一致, CUDA 后处理正确性验证通过!" << std::endl;
    return 0;
}
