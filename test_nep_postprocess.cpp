/**
 * @file test_nep_postprocess.cpp
 * @brief 独立单元测试: 验证 postprocess_nep_cpu 能量、力、应力后处理的正确性
 *
 * 本测试文件不依赖 ABACUS 完整构建系统, 仅需 matrix.h/matrix.cpp 和
 * esolver_nep_postprocess.h/esolver_nep_postprocess.cpp 即可编译运行.
 *
 * 编译方法:
 *   g++ -std=c++11 -I source -I source/source_base \
 *       test_nep_postprocess.cpp \
 *       source/source_esolver/esolver_nep_postprocess.cpp \
 *       source/source_base/matrix.cpp \
 *       -o test_nep_postprocess
 *
 * 运行:
 *   ./test_nep_postprocess
 */

#include "esolver_nep_postprocess.h"
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>

using ModuleESolver::postprocess_nep_cpu;

// 测试结果全局变量
int tests_passed = 0;
int tests_failed = 0;

// 浮点数比较精度 (考虑单位换算后的舍入误差)
// 对于大体系 (nat=1000), 累加累积误差约 5e-12, 使用 1e-10 作为阈值
const double epsilon = 1e-10;

void assert_double_eq(const std::string& label, double val, double expected)
{
    if (std::abs(val - expected) < epsilon)
    {
        std::cout << "  [PASS] " << label << ": " << val << " == " << expected << std::endl;
        tests_passed++;
    }
    else
    {
        std::cout << "  [FAIL] " << label << ": got " << val << ", expected " << expected
                  << " (diff=" << std::abs(val - expected) << ")" << std::endl;
        tests_failed++;
    }
}

// ============================================================================
// Test 1: 单原子 (nat=1) — 基础功能验证
// ============================================================================
void test_single_atom()
{
    std::cout << "\n=== Test 1: 单原子 (nat=1) ===" << std::endl;

    const int nat = 1;

    // 输入: 模拟 NEP 外部库返回的原始数据 (eV 单位, 假想值)
    std::vector<double> atomic_energy = {2.0};           // 每原子能量 = 2 eV
    // raw_force: SoA 布局 [fx0, fy0, fz0]
    std::vector<double> raw_force = {3.0, 4.0, 5.0};    // 力 (eV/A)
    // raw_virial: SoA 布局 [v0_0, v1_0, ..., v8_0]  每个分量 1 个原子值
    std::vector<double> raw_virial = {
        1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0
    };

    // 单位换算因子 (与 esolver_nep.cpp 中一致)
    const double fact_e = 1.0 / 13.605703976;  // Ry_to_eV 的倒数 -> eV -> Ry
    const double fact_f = 1.0 / (13.605703976 * 1.88972612546); // eV/A -> Ry/Bohr
    const double fact_v = 1.0 / 1.0;           // 简化: 体积归一化因子在此设为1

    double potential = 0.0;
    ModuleBase::matrix force(nat, 3);
    ModuleBase::matrix virial(3, 3);

    postprocess_nep_cpu(nat, atomic_energy.data(), raw_force.data(),
                        raw_virial.data(), fact_e, fact_f, fact_v,
                        potential, force, virial);

    // 验证能量 (Ry)
    double expected_energy = 2.0 * fact_e;
    assert_double_eq("energy", potential, expected_energy);

    // 验证力 (Ry/Bohr)
    assert_double_eq("force(0,0)=fx", force(0, 0), 3.0 * fact_f);
    assert_double_eq("force(0,1)=fy", force(0, 1), 4.0 * fact_f);
    assert_double_eq("force(0,2)=fz", force(0, 2), 5.0 * fact_f);

    // 验证 virial (单原子时求和等于自身)
    // virial_sum[j] = sum_over_atoms(raw_virial[j*nat + i])
    // virial(i,j)  = virial_sum[3*i + j] * fact_v
    assert_double_eq("virial(0,0)=v0", virial(0, 0), 1.0 * fact_v);
    assert_double_eq("virial(0,1)=v1", virial(0, 1), 2.0 * fact_v);
    assert_double_eq("virial(0,2)=v2", virial(0, 2), 3.0 * fact_v);
    assert_double_eq("virial(1,0)=v3", virial(1, 0), 4.0 * fact_v);
    assert_double_eq("virial(1,1)=v4", virial(1, 1), 5.0 * fact_v);
    assert_double_eq("virial(1,2)=v5", virial(1, 2), 6.0 * fact_v);
    assert_double_eq("virial(2,0)=v6", virial(2, 0), 7.0 * fact_v);
    assert_double_eq("virial(2,1)=v7", virial(2, 1), 8.0 * fact_v);
    assert_double_eq("virial(2,2)=v8", virial(2, 2), 9.0 * fact_v);
}

// ============================================================================
// Test 2: 多原子 (nat=4) — 验证能量求和和力的 SoA→行主序转换正确性
// ============================================================================
void test_multi_atom()
{
    std::cout << "\n=== Test 2: 多原子 (nat=4) — 能量求和与力格式转换 ===" << std::endl;

    const int nat = 4;

    // 每原子能量: 随意赋值
    std::vector<double> atomic_energy = {1.0, 2.0, 3.0, 4.0};
    // raw_force: SoA 布局 [fx0,fx1,fx2,fx3, fy0,fy1,fy2,fy3, fz0,fz1,fz2,fz3]
    std::vector<double> raw_force = {
        /* fx */ 1.0, 2.0, 3.0, 4.0,
        /* fy */ 5.0, 6.0, 7.0, 8.0,
        /* fz */ 9.0, 10.0, 11.0, 12.0
    };
    std::vector<double> raw_virial(9 * nat, 1.0); // 全部给 1.0

    const double fact_e = 2.0;
    const double fact_f = 0.5;
    const double fact_v = 3.0;

    double potential = 0.0;
    ModuleBase::matrix force(nat, 3);
    ModuleBase::matrix virial(3, 3);

    postprocess_nep_cpu(nat, atomic_energy.data(), raw_force.data(),
                        raw_virial.data(), fact_e, fact_f, fact_v,
                        potential, force, virial);

    // 验证能量求和: sum(1+2+3+4)*2.0 = 20.0
    double expected_potential = (1.0 + 2.0 + 3.0 + 4.0) * fact_e;
    assert_double_eq("energy sum nat=4", potential, expected_potential);

    // 验证力格式转换: SoA -> 行主序
    // force(i,0) = raw_force[i + 0*nat] * fact_f = raw_force[i] * fact_f
    // force(i,1) = raw_force[i + 1*nat] * fact_f
    // force(i,2) = raw_force[i + 2*nat] * fact_f
    for (int i = 0; i < nat; ++i)
    {
        assert_double_eq("force(" + std::to_string(i) + ",0)=fx" + std::to_string(i),
                         force(i, 0), raw_force[i] * fact_f);
        assert_double_eq("force(" + std::to_string(i) + ",1)=fy" + std::to_string(i),
                         force(i, 1), raw_force[i + nat] * fact_f);
        assert_double_eq("force(" + std::to_string(i) + ",2)=fz" + std::to_string(i),
                         force(i, 2), raw_force[i + 2 * nat] * fact_f);
    }

    // 验证 virial: 每个分量为 nat 个 1.0 求和 * fact_v = nat * fact_v
    double expected_v = nat * fact_v; // 4 * 3.0 = 12.0
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            assert_double_eq("virial(" + std::to_string(i) + "," + std::to_string(j) + ")",
                             virial(i, j), expected_v);
        }
    }
}

// ============================================================================
// Test 3: 零值输入 — 验证边界条件
// ============================================================================
void test_zero_input()
{
    std::cout << "\n=== Test 3: 零值输入 — 边界条件 ===" << std::endl;

    const int nat = 3;

    std::vector<double> atomic_energy(nat, 0.0);
    std::vector<double> raw_force(3 * nat, 0.0);
    std::vector<double> raw_virial(9 * nat, 0.0);

    const double fact_e = 1.0;
    const double fact_f = 1.0;
    const double fact_v = 1.0;

    double potential = -999.0; // 故意给非零值
    ModuleBase::matrix force(nat, 3);
    ModuleBase::matrix virial(3, 3);

    postprocess_nep_cpu(nat, atomic_energy.data(), raw_force.data(),
                        raw_virial.data(), fact_e, fact_f, fact_v,
                        potential, force, virial);

    assert_double_eq("zero energy", potential, 0.0);
    for (int i = 0; i < nat; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            assert_double_eq("zero force(" + std::to_string(i) + "," + std::to_string(j) + ")",
                             force(i, j), 0.0);
        }
    }
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            assert_double_eq("zero virial(" + std::to_string(i) + "," + std::to_string(j) + ")",
                             virial(i, j), 0.0);
        }
    }
}

// ============================================================================
// Test 4: 大体系 (nat=1000) — 验证数值稳定性和性能
// ============================================================================
void test_large_system()
{
    std::cout << "\n=== Test 4: 大体系 (nat=1000) — 数值稳定性 ===" << std::endl;

    const int nat = 1000;

    std::vector<double> atomic_energy(nat, 1.0);
    std::vector<double> raw_force(3 * nat, 0.5);
    std::vector<double> raw_virial(9 * nat, 0.1);

    const double fact_e = 2.0;
    const double fact_f = 3.0;
    const double fact_v = 4.0;

    double potential = 0.0;
    ModuleBase::matrix force(nat, 3);
    ModuleBase::matrix virial(3, 3);

    postprocess_nep_cpu(nat, atomic_energy.data(), raw_force.data(),
                        raw_virial.data(), fact_e, fact_f, fact_v,
                        potential, force, virial);

    // 能量: nat * 1.0 * 2.0 = 2000.0
    assert_double_eq("energy nat=1000", potential, nat * 1.0 * fact_e);

    // 力: 每个分量都是 0.5 * 3.0 = 1.5
    for (int i = 0; i < nat; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            assert_double_eq("force(" + std::to_string(i) + "," + std::to_string(j) + ")",
                             force(i, j), 0.5 * fact_f);
        }
    }

    // virial: nat * 0.1 * 4.0 = 400.0 每个分量
    double expected_v = nat * 0.1 * fact_v;
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            assert_double_eq("virial(" + std::to_string(i) + "," + std::to_string(j) + ")",
                             virial(i, j), expected_v);
        }
    }
}

// ============================================================================
// Test 5: SoA 数据布局交叉验证 — 直接构造参考数据对比
// ============================================================================
void test_soa_layout()
{
    std::cout << "\n=== Test 5: SoA 数据布局交叉验证 ===" << std::endl;

    const int nat = 3;

    // 构造一个易于验证的力数据: 每个原子 i 的力为 (i*1.0+0.1, i*1.0+0.2, i*1.0+0.3)
    // SoA 格式: fx[0..2] = [0.1, 1.1, 2.1], fy[0..2] = [0.2, 1.2, 2.2], fz[0..2] = [0.3, 1.3, 2.3]
    std::vector<double> atomic_energy(nat, 0.0);
    std::vector<double> raw_force(3 * nat);
    std::vector<double> raw_virial(9 * nat, 0.0);

    for (int i = 0; i < nat; ++i)
    {
        raw_force[i] = i * 1.0 + 0.1;            // fx[i]
        raw_force[i + nat] = i * 1.0 + 0.2;       // fy[i]
        raw_force[i + 2 * nat] = i * 1.0 + 0.3;   // fz[i]
    }

    const double fact_e = 1.0;
    const double fact_f = 1.0;
    const double fact_v = 1.0;

    double potential = 0.0;
    ModuleBase::matrix force(nat, 3);
    ModuleBase::matrix virial(3, 3);

    postprocess_nep_cpu(nat, atomic_energy.data(), raw_force.data(),
                        raw_virial.data(), fact_e, fact_f, fact_v,
                        potential, force, virial);

    // 行主序 force(i,j) 索引 = i * 3 + j
    for (int i = 0; i < nat; ++i)
    {
        assert_double_eq("SoA force(" + std::to_string(i) + ",0)",
                         force(i, 0), i * 1.0 + 0.1);
        assert_double_eq("SoA force(" + std::to_string(i) + ",1)",
                         force(i, 1), i * 1.0 + 0.2);
        assert_double_eq("SoA force(" + std::to_string(i) + ",2)",
                         force(i, 2), i * 1.0 + 0.3);
    }
}

// ============================================================================
// Test 6: Virial SoA 布局验证 — 构造不同偏移值测试
// ============================================================================
void test_virial_soa()
{
    std::cout << "\n=== Test 6: Virial SoA 布局验证 ===" << std::endl;

    const int nat = 2;

    std::vector<double> atomic_energy(nat, 0.0);
    std::vector<double> raw_force(3 * nat, 0.0);

    // virial: 9 个分量, 每个分量 nat 个原子值
    // 第 j 个分量: [j*10 + 1, j*10 + 2]
    std::vector<double> raw_virial(9 * nat);
    for (int j = 0; j < 9; ++j)
    {
        for (int i = 0; i < nat; ++i)
        {
            raw_virial[j * nat + i] = j * 10.0 + (i + 1.0);
        }
    }

    const double fact_e = 1.0;
    const double fact_f = 1.0;
    const double fact_v = 1.0; // fact_v = 1 便于验证

    double potential = 0.0;
    ModuleBase::matrix force(nat, 3);
    ModuleBase::matrix virial(3, 3);

    postprocess_nep_cpu(nat, atomic_energy.data(), raw_force.data(),
                        raw_virial.data(), fact_e, fact_f, fact_v,
                        potential, force, virial);

    // virial(i,j) = sum_of_atoms(raw_virial[k*nat + :])  其中 k=3*i+j
    // 对于 nat=2, sum_of_atoms = (k*10+1) + (k*10+2) = k*20 + 3
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            int k = 3 * i + j;
            double expected = k * 20.0 + 3.0;
            assert_double_eq("virial(" + std::to_string(i) + "," + std::to_string(j)
                             + ") k=" + std::to_string(k),
                             virial(i, j), expected);
        }
    }
}

// ============================================================================
// main
// ============================================================================
int main()
{
    std::cout << "======================================================" << std::endl;
    std::cout << "  NEP Postprocess CPU Unit Test" << std::endl;
    std::cout << "  测试 esolver_nep_postprocess.cpp 中的" << std::endl;
    std::cout << "  postprocess_nep_cpu 函数" << std::endl;
    std::cout << "======================================================" << std::endl;

    test_single_atom();
    test_multi_atom();
    test_zero_input();
    test_large_system();
    test_soa_layout();
    test_virial_soa();

    std::cout << "\n======================================================" << std::endl;
    std::cout << "  Results: " << tests_passed << " passed, "
              << tests_failed << " failed" << std::endl;
    std::cout << "======================================================" << std::endl;

    if (tests_failed > 0)
    {
        std::cerr << "\n[FAIL] 存在 " << tests_failed << " 项测试失败!" << std::endl;
        return 1;
    }

    std::cout << "\n[PASS] 所有单元测试通过!" << std::endl;
    return 0;
}
