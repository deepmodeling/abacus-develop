/**
 * @file diago_mixed_precision_benchmark.cpp
 * @brief 混合精度特征值求解器的性能基准测试和正确性验证
 *
 * 测试内容:
 *   1. 不同矩阵尺寸下的性能对比（单精度/双精度/混合精度）
 *   2. 混合精度 vs 双精度的精度验证（误差 < 1e-6）
 *   3. 不同精度组合的正确性测试
 *   4. 边界情况测试（小矩阵、病态矩阵、不同稀疏度）
 */

#include "gtest/gtest.h"
#include "source_hsolver/diago_cg.h"
#include "source_hsolver/diago_david.h"
#include <complex>
#include <random>
#include <vector>
#include <chrono>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>

using Complex = std::complex<double>;
using namespace hsolver;

// ============================================================================
// 辅助函数
// ============================================================================

/// 生成随机 Hermitian 矩阵
static void make_hermitian(int n, std::vector<Complex>& H, unsigned seed = 12345)
{
    H.resize(static_cast<size_t>(n) * n);
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j <= i; ++j)
        {
            const double real = dist(rng);
            const double imag = (i == j ? 0.0 : dist(rng));
            H[static_cast<size_t>(i) * n + j] = Complex(real, imag);
            H[static_cast<size_t>(j) * n + i] = std::conj(H[static_cast<size_t>(i) * n + j]);
        }
    }
}

/// 生成具有可调条件数的 Hermitian 矩阵
static void make_hermitian_conditioned(int n, std::vector<Complex>& H, double cond_num, unsigned seed = 12345)
{
    H.resize(static_cast<size_t>(n) * n);
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);

    // 生成随机对角占优矩阵
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j <= i; ++j)
        {
            const double real = dist(rng);
            const double imag = (i == j ? 0.0 : dist(rng) * 0.1);
            H[static_cast<size_t>(i) * n + j] = Complex(real, imag);
            H[static_cast<size_t>(j) * n + i] = std::conj(H[static_cast<size_t>(i) * n + j]);
        }
    }

    // 调整对角元素以控制条件数
    double diag_scale = cond_num / n;
    for (int i = 0; i < n; ++i)
    {
        H[static_cast<size_t>(i) * n + i] += Complex(i * diag_scale, 0.0);
    }
}

/// 生成随机初始波函数
static void make_random_psi(int nband, int dim, std::vector<Complex>& psi, unsigned seed = 54321)
{
    psi.resize(static_cast<size_t>(nband) * dim);
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> dist(-0.5, 0.5);
    for (size_t i = 0; i < psi.size(); ++i)
    {
        psi[i] = Complex(dist(rng), dist(rng));
    }
}

/// 应用 Hamiltonian 矩阵
static void apply_hamiltonian(const std::vector<Complex>& H, int n,
                               const Complex* psi_in, Complex* hpsi_out,
                               int ld, int nvec)
{
    for (int v = 0; v < nvec; ++v)
    {
        const Complex* psi_vec = psi_in + static_cast<size_t>(v) * ld;
        Complex* out_vec = hpsi_out + static_cast<size_t>(v) * ld;
        for (int i = 0; i < n; ++i)
        {
            Complex sum = 0.0;
            for (int j = 0; j < n; ++j)
            {
                sum += H[static_cast<size_t>(i) * n + j] * psi_vec[j];
            }
            out_vec[i] = sum;
        }
    }
}

/// 恒等重叠矩阵
static void apply_overlap(const Complex* psi_in, Complex* spsi_out, int ld, int nvec)
{
    for (int i = 0; i < static_cast<size_t>(nvec) * ld; ++i)
    {
        spsi_out[i] = psi_in[i];
    }
}

/// 使用 LAPACK 计算参考特征值 (简化版: 仅用前 nband 个)
static std::vector<double> compute_reference_eigenvalues(const std::vector<Complex>& H, int n, int nband)
{
    // 复制 H 用于 LAPACK (zheev 会修改矩阵)
    std::vector<Complex> H_copy = H;
    std::vector<double> eigenvalues(n, 0.0);

    int lwork = 2 * n;
    std::vector<Complex> work(lwork);
    std::vector<double> rwork(3 * n - 2);
    int info = 0;
    char jobz = 'N'; // 只需要特征值
    char uplo = 'U';

    zheev_(&jobz, &uplo, &n, H_copy.data(), &n, eigenvalues.data(), work.data(), &lwork, rwork.data(), &info);

    if (info != 0)
    {
        std::cerr << "LAPACK zheev failed with info=" << info << std::endl;
    }

    // 返回前 nband 个特征值（zheev 返回升序排列）
    return std::vector<double>(eigenvalues.begin(), eigenvalues.begin() + nband);
}

/// 计时器辅助类
class ScopedTimer
{
  public:
    ScopedTimer(double& elapsed) : elapsed_(elapsed), start_(std::chrono::high_resolution_clock::now()) {}

    ~ScopedTimer()
    {
        auto end = std::chrono::high_resolution_clock::now();
        elapsed_ = std::chrono::duration<double>(end - start_).count();
    }

  private:
    double& elapsed_;
    std::chrono::high_resolution_clock::time_point start_;
};

// ============================================================================
// 测试 1: 混合精度正确性 — 不同矩阵尺寸
// ============================================================================

class MixedPrecisionCorrectnessTest : public ::testing::TestWithParam<int>
{
};

TEST_P(MixedPrecisionCorrectnessTest, CGMixedPrecisionMatchesDouble)
{
    const int dim = GetParam();
    const int nband = std::min(dim / 2, 8);
    const int ld_psi = dim;

    std::vector<Complex> H;
    make_hermitian(dim, H, 12345);

    std::vector<Complex> psi_initial;
    make_random_psi(nband, dim, psi_initial, 54321);

    std::vector<Complex> psi_double = psi_initial;
    std::vector<Complex> psi_mixed = psi_initial;
    std::vector<double> eigen_double(nband, 0.0);
    std::vector<double> eigen_mixed(nband, 0.0);

    auto hpsi_func = [&H, dim](Complex* psi_in, Complex* hpsi_out, const int ld, const int nvec) {
        apply_hamiltonian(H, dim, psi_in, hpsi_out, ld, nvec);
    };
    auto spsi_func = [](Complex* psi_in, Complex* spsi_out, const int ld, const int nvec) {
        apply_overlap(psi_in, spsi_out, ld, nvec);
    };

    std::vector<double> ethr_band(nband, 1e-6);

    // 双精度
    DiagoCG<Complex> cg_double("pw", "nscf", false, DiagoCG<Complex>::SubspaceFunc(), 1e-6, 200, 1, PrecisionMode::kDouble);
    cg_double.diag(hpsi_func, spsi_func, ld_psi, nband, dim,
                   psi_double.data(), eigen_double.data(), ethr_band, nullptr);

    // 混合精度
    DiagoCG<Complex> cg_mixed("pw", "nscf", false, DiagoCG<Complex>::SubspaceFunc(), 1e-6, 200, 1, PrecisionMode::kMixed);
    cg_mixed.diag(hpsi_func, spsi_func, ld_psi, nband, dim,
                  psi_mixed.data(), eigen_mixed.data(), ethr_band, nullptr);

    // 验证特征值一致性
    for (int i = 0; i < nband; ++i)
    {
        EXPECT_NEAR(eigen_double[i], eigen_mixed[i], 1e-6)
            << "Dim=" << dim << " Band=" << i
            << " double=" << eigen_double[i] << " mixed=" << eigen_mixed[i];
    }
}

INSTANTIATE_TEST_SUITE_P(VariousDimensions,
                         MixedPrecisionCorrectnessTest,
                         ::testing::Values(8, 16, 32, 64, 128));

// ============================================================================
// 测试 2: David 求解器混合精度正确性
// ============================================================================

class DavidMixedPrecisionTest : public ::testing::TestWithParam<int>
{
};

TEST_P(DavidMixedPrecisionTest, DavidMixedPrecisionMatchesDouble)
{
    const int dim = GetParam();
    const int nband = std::min(dim / 2, 8);
    const int ld_psi = dim;
    const int david_ndim = 4;

    std::vector<Complex> H;
    make_hermitian(dim, H, 23456);

    std::vector<Complex> psi_initial;
    make_random_psi(nband, dim, psi_initial, 65432);

    std::vector<Complex> psi_double = psi_initial;
    std::vector<Complex> psi_mixed = psi_initial;
    std::vector<double> eigen_double(nband, 0.0);
    std::vector<double> eigen_mixed(nband, 0.0);

    auto hpsi_func = [&H, dim](Complex* psi_in, Complex* hpsi_out, const int ld, const int nvec) {
        apply_hamiltonian(H, dim, psi_in, hpsi_out, ld, nvec);
    };
    auto spsi_func = [](Complex* psi_in, Complex* spsi_out, const int ld, const int nvec) {
        apply_overlap(psi_in, spsi_out, ld, nvec);
    };

    std::vector<double> ethr_band(nband, 1e-6);
    std::vector<double> precondition(dim, 1.0);

    diag_comm_info comm_info = {0, 1};

    // 双精度
    DiagoDavid<Complex> dav_double(precondition.data(), nband, dim, david_ndim, false, comm_info, PrecisionMode::kDouble);
    dav_double.diag(hpsi_func, spsi_func, ld_psi, psi_double.data(), eigen_double.data(),
                    ethr_band, 100, 5, 0);

    // 混合精度
    DiagoDavid<Complex> dav_mixed(precondition.data(), nband, dim, david_ndim, false, comm_info, PrecisionMode::kMixed);
    dav_mixed.diag(hpsi_func, spsi_func, ld_psi, psi_mixed.data(), eigen_mixed.data(),
                   ethr_band, 100, 5, 0);

    // 验证
    for (int i = 0; i < nband; ++i)
    {
        EXPECT_NEAR(eigen_double[i], eigen_mixed[i], 1e-5)
            << "David Dim=" << dim << " Band=" << i
            << " double=" << eigen_double[i] << " mixed=" << eigen_mixed[i];
    }
}

INSTANTIATE_TEST_SUITE_P(DavidVariousDimensions,
                         DavidMixedPrecisionTest,
                         ::testing::Values(8, 16, 32, 64));

// ============================================================================
// 测试 3: 性能基准测试
// ============================================================================

TEST(MixedPrecisionBenchmark, PerformanceComparison)
{
    const int dim = 128;
    const int nband = 8;
    const int ld_psi = dim;

    std::vector<Complex> H;
    make_hermitian(dim, H, 34567);

    auto hpsi_func = [&H, dim](Complex* psi_in, Complex* hpsi_out, const int ld, const int nvec) {
        apply_hamiltonian(H, dim, psi_in, hpsi_out, ld, nvec);
    };
    auto spsi_func = [](Complex* psi_in, Complex* spsi_out, const int ld, const int nvec) {
        apply_overlap(psi_in, spsi_out, ld, nvec);
    };
    std::vector<double> ethr_band(nband, 1e-6);

    std::vector<double> times(3, 0.0);
    std::vector<double> eigen_results[3];
    for (int i = 0; i < 3; ++i)
    {
        eigen_results[i].resize(nband);
    }

    // 双精度
    {
        std::vector<Complex> psi(nband * dim);
        make_random_psi(nband, dim, psi, 11111);
        DiagoCG<Complex> cg("pw", "nscf", false, DiagoCG<Complex>::SubspaceFunc(), 1e-6, 200, 1, PrecisionMode::kDouble);
        double elapsed = 0.0;
        {
            ScopedTimer timer(elapsed);
            cg.diag(hpsi_func, spsi_func, ld_psi, nband, dim, psi.data(), eigen_results[0].data(), ethr_band, nullptr);
        }
        times[0] = elapsed;
        std::cout << "[Benchmark] Double precision: " << elapsed << " s" << std::endl;
    }

    // 单精度
    {
        std::vector<Complex> psi(nband * dim);
        make_random_psi(nband, dim, psi, 11111);
        DiagoCG<Complex> cg("pw", "nscf", false, DiagoCG<Complex>::SubspaceFunc(), 1e-6, 200, 1, PrecisionMode::kFloat);
        double elapsed = 0.0;
        {
            ScopedTimer timer(elapsed);
            cg.diag(hpsi_func, spsi_func, ld_psi, nband, dim, psi.data(), eigen_results[1].data(), ethr_band, nullptr);
        }
        times[1] = elapsed;
        std::cout << "[Benchmark] Float precision:  " << elapsed << " s" << std::endl;
    }

    // 混合精度
    {
        std::vector<Complex> psi(nband * dim);
        make_random_psi(nband, dim, psi, 11111);
        DiagoCG<Complex> cg("pw", "nscf", false, DiagoCG<Complex>::SubspaceFunc(), 1e-6, 200, 1, PrecisionMode::kMixed);
        double elapsed = 0.0;
        {
            ScopedTimer timer(elapsed);
            cg.diag(hpsi_func, spsi_func, ld_psi, nband, dim, psi.data(), eigen_results[2].data(), ethr_band, nullptr);
        }
        times[2] = elapsed;
        std::cout << "[Benchmark] Mixed precision:  " << elapsed << " s" << std::endl;
    }

    // 计算加速比
    std::cout << "[Benchmark] Speedup (mixed/double): " << times[0] / times[2] << "x" << std::endl;
    std::cout << "[Benchmark] Speedup (float/double): " << times[0] / times[1] << "x" << std::endl;

    // 验证混合精度结果与双精度一致
    for (int i = 0; i < nband; ++i)
    {
        EXPECT_NEAR(eigen_results[0][i], eigen_results[2][i], 1e-6)
            << "Mixed precision eigenvalue differs at band " << i;
    }
}

// ============================================================================
// 测试 4: 精度切换边界情况
// ============================================================================

TEST(MixedPrecisionEdgeCases, SmallMatrix)
{
    // 测试 2x2 极小矩阵
    const int dim = 2;
    const int nband = 1;
    const int ld_psi = dim;

    std::vector<Complex> H = {Complex(1.0, 0.0), Complex(0.5, 0.1),
                               Complex(0.5, -0.1), Complex(2.0, 0.0)};

    std::vector<Complex> psi_double = {Complex(1.0, 0.0), Complex(0.0, 0.0)};
    std::vector<Complex> psi_mixed = {Complex(1.0, 0.0), Complex(0.0, 0.0)};
    std::vector<double> eigen_double(1, 0.0);
    std::vector<double> eigen_mixed(1, 0.0);

    auto hpsi_func = [&H, dim](Complex* psi_in, Complex* hpsi_out, const int ld, const int nvec) {
        apply_hamiltonian(H, dim, psi_in, hpsi_out, ld, nvec);
    };
    auto spsi_func = [](Complex* psi_in, Complex* spsi_out, const int ld, const int nvec) {
        apply_overlap(psi_in, spsi_out, ld, nvec);
    };
    std::vector<double> ethr_band(1, 1e-8);

    DiagoCG<Complex> cg_double("pw", "nscf", false, DiagoCG<Complex>::SubspaceFunc(), 1e-8, 200, 1, PrecisionMode::kDouble);
    cg_double.diag(hpsi_func, spsi_func, ld_psi, nband, dim,
                   psi_double.data(), eigen_double.data(), ethr_band, nullptr);

    DiagoCG<Complex> cg_mixed("pw", "nscf", false, DiagoCG<Complex>::SubspaceFunc(), 1e-8, 200, 1, PrecisionMode::kMixed);
    cg_mixed.diag(hpsi_func, spsi_func, ld_psi, nband, dim,
                  psi_mixed.data(), eigen_mixed.data(), ethr_band, nullptr);

    EXPECT_NEAR(eigen_double[0], eigen_mixed[0], 1e-6);
}

TEST(MixedPrecisionEdgeCases, IllConditionedMatrix)
{
    // 测试条件数较大的矩阵
    const int dim = 32;
    const int nband = 4;
    const int ld_psi = dim;

    std::vector<Complex> H;
    make_hermitian_conditioned(dim, H, 1e4, 99999);

    std::vector<Complex> psi_initial;
    make_random_psi(nband, dim, psi_initial, 77777);

    std::vector<Complex> psi_double = psi_initial;
    std::vector<Complex> psi_mixed = psi_initial;
    std::vector<double> eigen_double(nband, 0.0);
    std::vector<double> eigen_mixed(nband, 0.0);

    auto hpsi_func = [&H, dim](Complex* psi_in, Complex* hpsi_out, const int ld, const int nvec) {
        apply_hamiltonian(H, dim, psi_in, hpsi_out, ld, nvec);
    };
    auto spsi_func = [](Complex* psi_in, Complex* spsi_out, const int ld, const int nvec) {
        apply_overlap(psi_in, spsi_out, ld, nvec);
    };
    std::vector<double> ethr_band(nband, 1e-5);

    DiagoCG<Complex> cg_double("pw", "nscf", false, DiagoCG<Complex>::SubspaceFunc(), 1e-5, 500, 1, PrecisionMode::kDouble);
    cg_double.diag(hpsi_func, spsi_func, ld_psi, nband, dim,
                   psi_double.data(), eigen_double.data(), ethr_band, nullptr);

    DiagoCG<Complex> cg_mixed("pw", "nscf", false, DiagoCG<Complex>::SubspaceFunc(), 1e-5, 500, 1, PrecisionMode::kMixed);
    cg_mixed.diag(hpsi_func, spsi_func, ld_psi, nband, dim,
                  psi_mixed.data(), eigen_mixed.data(), ethr_band, nullptr);

    for (int i = 0; i < nband; ++i)
    {
        EXPECT_NEAR(eigen_double[i], eigen_mixed[i], 1e-5)
            << "Ill-conditioned matrix, band " << i;
    }
}

// ============================================================================
// 测试 5: 不同精度模式组合
// ============================================================================

TEST(MixedPrecisionCombinations, AllPrecisionModesCG)
{
    const int dim = 24;
    const int nband = 4;
    const int ld_psi = dim;

    std::vector<Complex> H;
    make_hermitian(dim, H, 11111);

    auto hpsi_func = [&H, dim](Complex* psi_in, Complex* hpsi_out, const int ld, const int nvec) {
        apply_hamiltonian(H, dim, psi_in, hpsi_out, ld, nvec);
    };
    auto spsi_func = [](Complex* psi_in, Complex* spsi_out, const int ld, const int nvec) {
        apply_overlap(psi_in, spsi_out, ld, nvec);
    };
    std::vector<double> ethr_band(nband, 1e-6);

    std::vector<double> eigen_double(nband, 0.0);
    std::vector<double> eigen_float(nband, 0.0);
    std::vector<double> eigen_mixed(nband, 0.0);

    // Double
    {
        std::vector<Complex> psi(nband * dim);
        make_random_psi(nband, dim, psi, 22222);
        DiagoCG<Complex> cg("pw", "nscf", false, DiagoCG<Complex>::SubspaceFunc(), 1e-6, 200, 1, PrecisionMode::kDouble);
        cg.diag(hpsi_func, spsi_func, ld_psi, nband, dim, psi.data(), eigen_double.data(), ethr_band, nullptr);
    }
    // Float
    {
        std::vector<Complex> psi(nband * dim);
        make_random_psi(nband, dim, psi, 22222);
        DiagoCG<Complex> cg("pw", "nscf", false, DiagoCG<Complex>::SubspaceFunc(), 1e-6, 200, 1, PrecisionMode::kFloat);
        cg.diag(hpsi_func, spsi_func, ld_psi, nband, dim, psi.data(), eigen_float.data(), ethr_band, nullptr);
    }
    // Mixed
    {
        std::vector<Complex> psi(nband * dim);
        make_random_psi(nband, dim, psi, 22222);
        DiagoCG<Complex> cg("pw", "nscf", false, DiagoCG<Complex>::SubspaceFunc(), 1e-6, 200, 1, PrecisionMode::kMixed);
        cg.diag(hpsi_func, spsi_func, ld_psi, nband, dim, psi.data(), eigen_mixed.data(), ethr_band, nullptr);
    }

    // Mixed 应在误差范围内匹配 Double
    for (int i = 0; i < nband; ++i)
    {
        EXPECT_NEAR(eigen_double[i], eigen_mixed[i], 1e-6)
            << "Mixed vs Double, band " << i;
    }

    // Float 可能有较大误差，但仍应在合理范围内
    for (int i = 0; i < nband; ++i)
    {
        double rel_err = std::abs(eigen_double[i] - eigen_float[i])
                         / std::max(1.0, std::abs(eigen_double[i]));
        EXPECT_LT(rel_err, 1e-3)
            << "Float vs Double relative error too large, band " << i
            << " rel_err=" << rel_err;
    }
}

// ============================================================================
// 测试 6: 收敛性验证
// ============================================================================

TEST(MixedPrecisionConvergence, ConvergenceTest)
{
    const int dim = 48;
    const int nband = 6;
    const int ld_psi = dim;

    std::vector<Complex> H;
    make_hermitian(dim, H, 33333);

    auto ref_eigen = compute_reference_eigenvalues(H, dim, nband);

    auto hpsi_func = [&H, dim](Complex* psi_in, Complex* hpsi_out, const int ld, const int nvec) {
        apply_hamiltonian(H, dim, psi_in, hpsi_out, ld, nvec);
    };
    auto spsi_func = [](Complex* psi_in, Complex* spsi_out, const int ld, const int nvec) {
        apply_overlap(psi_in, spsi_out, ld, nvec);
    };

    // 测试不同收敛阈值
    std::vector<double> thresholds = {1e-3, 1e-4, 1e-5, 1e-6};

    for (double thr : thresholds)
    {
        std::vector<double> ethr_band(nband, thr);
        std::vector<Complex> psi(nband * dim);
        make_random_psi(nband, dim, psi, 44444);
        std::vector<double> eigen_mixed(nband, 0.0);

        DiagoCG<Complex> cg("pw", "nscf", false, DiagoCG<Complex>::SubspaceFunc(), thr, 500, 1, PrecisionMode::kMixed);
        cg.diag(hpsi_func, spsi_func, ld_psi, nband, dim, psi.data(), eigen_mixed.data(), ethr_band, nullptr);

        for (int i = 0; i < nband; ++i)
        {
            double abs_err = std::abs(eigen_mixed[i] - ref_eigen[i]);
            EXPECT_LT(abs_err, std::max(thr * 10.0, 1e-5))
                << "Threshold=" << thr << " Band=" << i
                << " abs_err=" << abs_err
                << " mixed=" << eigen_mixed[i] << " ref=" << ref_eigen[i];
        }
    }
}

// ============================================================================
// 测试 7: 解析精度模式函数
// ============================================================================

TEST(PrecisionModeParsing, ParsePrecisionModeString)
{
    EXPECT_EQ(parse_precision_mode("double"), PrecisionMode::kDouble);
    EXPECT_EQ(parse_precision_mode("float"), PrecisionMode::kFloat);
    EXPECT_EQ(parse_precision_mode("single"), PrecisionMode::kFloat);
    EXPECT_EQ(parse_precision_mode("mixed"), PrecisionMode::kMixed);
    EXPECT_EQ(parse_precision_mode("auto"), PrecisionMode::kMixed);
    EXPECT_EQ(parse_precision_mode("unknown"), PrecisionMode::kDouble); // default
    EXPECT_EQ(parse_precision_mode(""), PrecisionMode::kDouble);
}

TEST(PrecisionModeToString, ConvertToString)
{
    EXPECT_EQ(precision_mode_to_string(PrecisionMode::kDouble), "double");
    EXPECT_EQ(precision_mode_to_string(PrecisionMode::kFloat), "float");
    EXPECT_EQ(precision_mode_to_string(PrecisionMode::kMixed), "mixed");
}
