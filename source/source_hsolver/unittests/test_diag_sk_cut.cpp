#include "source_hsolver/diag_sk_cut.h"

#include "source_base/matrix_block.h"
#include "source_base/parallel_blkcyc.h"
#include "source_base/parallel_common.h"
#include "source_base/parallel_global.h"
#include "source_base/parallel_reduce.h"
#include "source_hsolver/module_linalg/overlap_spectrum.h"
#include "source_psi/psi.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <string>
#include <vector>

#include <gtest/gtest.h>

namespace
{
TEST(OverlapSpectrum, ValidatesAndScalesRetainedModes)
{
    std::vector<double> weights;
    EXPECT_TRUE(hsolver::overlap_inverse_sqrt({1e-32, 1.0, 4.0}, 0.0, weights));
    ASSERT_EQ(weights.size(), 3);
    EXPECT_NEAR(weights[0], 1e16, 2);
    EXPECT_DOUBLE_EQ(weights[1], 1.0);
    EXPECT_DOUBLE_EQ(weights[2], 0.5);

    EXPECT_TRUE(hsolver::overlap_inverse_sqrt({-1e-15, 0.25, 1.0, 4.0}, 0.0625, weights));
    EXPECT_EQ(weights, (std::vector<double>{0.0, 0.0, 1.0, 0.5}));
    EXPECT_FALSE(hsolver::overlap_inverse_sqrt({0.0, 4.0}, 0.0, weights));
    EXPECT_FALSE(hsolver::overlap_inverse_sqrt({-1.0, 0.0}, 0.2, weights));
    EXPECT_FALSE(hsolver::overlap_inverse_sqrt({2.0, 1.0}, 0.1, weights));
    const double nan = std::numeric_limits<double>::quiet_NaN();
    EXPECT_FALSE(hsolver::overlap_inverse_sqrt({nan, 1.0}, 0.1, weights));
    EXPECT_TRUE(weights.empty());
}

template <typename T>
void check_projection(const std::string& method, bool distributed, T phase, int retained_rank, int block)
{
    const int order = 4;
    const int bands = retained_rank;
    const double small = 1e-8;
    const double small_h = -20.0 * small;
    const double large = 3.0;
    const double large_h = 3.5 * large;
    const double off_diag = 0.25;
    const double off_diag_h = off_diag * std::sqrt(2.0);
    const T phase_conj = T(1) / phase;
    std::vector<T> h_full(order * order);
    std::vector<T> s_full(order * order);

    // Rotate overlap modes of eigenvalues small and large into AO rows 0 and 3.
    h_full[0] = h_full[15] = (small_h + large_h) / 2.0;
    s_full[0] = s_full[15] = (small + large) / 2.0;
    h_full[12] = phase_conj * ((small_h - large_h) / 2.0);
    s_full[12] = phase_conj * ((small - large) / 2.0);
    h_full[3] = phase * ((small_h - large_h) / 2.0);
    s_full[3] = phase * ((small - large) / 2.0);
    s_full[5] = 1.0;
    s_full[10] = 2.0;
    h_full[5] = 1.0;
    h_full[10] = 4.0;
    h_full[9] = h_full[6] = off_diag_h;

    Parallel_Blkcyclic input_layout;
    Parallel_Blkcyclic output_layout;
    const int* desc = nullptr;
    if (distributed)
    {
        input_layout.init(order, order, block, MPI_COMM_WORLD, false);
        const int context = input_layout.ictxt();
        output_layout.bind_descriptor(order, bands, block, context);
        desc = input_layout.desc();
    }
    else
    {
        input_layout.set_serial(order, order);
        output_layout.set_serial(order, bands);
    }

    const int rows = input_layout.nrloc();
    const int cols = input_layout.ncloc();
    const std::size_t local_rows = rows;
    const std::size_t local_cols = cols;
    std::vector<T> h(input_layout.size());
    std::vector<T> s(input_layout.size());
    const double nan = std::numeric_limits<double>::quiet_NaN();
    for (int col = 0; col < cols; ++col)
    {
        const int gc = input_layout.icglb(col);
        for (int row = 0; row < rows; ++row)
        {
            const int gr = input_layout.irglb(row);
            const int global_idx = gr + gc * order;
            const int local_idx = row + col * rows;
            h[local_idx] = gr <= gc ? h_full[global_idx] : T(nan);
            s[local_idx] = gr <= gc ? s_full[global_idx] : T(nan);
        }
    }
    ModuleBase::MatrixBlock<T> hk{h.data(), local_rows, local_cols, desc};
    ModuleBase::MatrixBlock<T> sk{s.data(), local_rows, local_cols, desc};
    const std::size_t local_output = std::max<std::size_t>(1, output_layout.size());
    std::vector<T> vectors(local_output);
    std::vector<double> values(bands);
    const double cutoffs[] = {0.8, 0.55, 1e-6};
    const double cutoff = cutoffs[retained_rank - 1];
    const int local_bands = output_layout.ncloc();
    psi::Psi<T> psi(vectors.data(), 1, local_bands, rows, rows, true);
    const Parallel_Blkcyclic* pool = distributed ? &input_layout : nullptr;
    const hsolver::LcaoDiagOptions options{method, order, bands, false, pool, true};
    hsolver::diag_sk_cut(hk, sk, psi, values.data(), options, 0, cutoff, nullptr);

    if (retained_rank == 1)
    {
        EXPECT_NEAR(values[0], large_h / large, 2e-12);
    }
    else if (retained_rank == 2)
    {
        EXPECT_NEAR(values[0], 2.0, 2e-12);
        EXPECT_NEAR(values[1], large_h / large, 2e-12);
    }
    else
    {
        const double splitting = std::hypot(0.5, off_diag);
        EXPECT_NEAR(values[0], 1.5 - splitting, 2e-12);
        EXPECT_NEAR(values[1], 1.5 + splitting, 2e-12);
        EXPECT_NEAR(values[2], large_h / large, 2e-12);
    }
    std::vector<T> c(order * bands);
    for (int col = 0; col < output_layout.ncloc(); ++col)
    {
        const int gc = output_layout.icglb(col);
        for (int row = 0; row < rows; ++row)
        {
            const int gr = output_layout.irglb(row);
            c[gr + gc * order] = vectors[row + col * rows];
        }
    }
    if (distributed)
    {
        const int count = c.size();
        Parallel_Reduce::reduce_all(c.data(), count);
    }

    for (int j = 0; j < bands; ++j)
    {
        EXPECT_NEAR(std::abs(c[j * order] + phase_conj * c[3 + j * order]), 0, 1e-12);
        for (int i = 0; i < bands; ++i)
        {
            std::complex<double> gram(0.0);
            std::complex<double> projected_h(0.0);
            for (int a = 0; a < order; ++a)
            {
                for (int b = 0; b < order; ++b)
                {
                    const auto left = std::conj(std::complex<double>(c[a + i * order]));
                    const auto right = c[b + j * order];
                    gram += left * s_full[a + b * order] * right;
                    projected_h += left * h_full[a + b * order] * right;
                }
            }
            const double expected_gram = i == j ? 1.0 : 0.0;
            const double expected_h = i == j ? values[i] : 0.0;
            EXPECT_NEAR(std::abs(gram - expected_gram), 0, 2e-12);
            EXPECT_NEAR(std::abs(projected_h - expected_h), 0, 2e-12);
        }
    }
}

TEST(DiagSkCut, ReplicatedRealAndComplex)
{
    check_projection<double>("lapack", false, 1.0, 3, 1);
    using Complex = std::complex<double>;
    const Complex phase(0, 1);
    check_projection<Complex>("lapack", false, phase, 3, 1);
}

TEST(DiagSkCut, DistributedRealAndComplex)
{
    if (Parallel_Common::comm_size(MPI_COMM_WORLD) > 4)
    {
        GTEST_SKIP() << "Four-orbital test uses at most four ranks.";
    }
    check_projection<double>("scalapack_gvx", true, 1.0, 3, 1);
    using Complex = std::complex<double>;
    const Complex phase(0, 1);
    check_projection<Complex>("scalapack_gvx", true, phase, 3, 1);
}
#ifdef __ELPA
void check_full_rank_elpa_fallback(const std::string& method, bool compact)
{
    SCOPED_TRACE(method);
    SCOPED_TRACE(compact);
    constexpr int order = 4;
    constexpr int bands = 3;
    constexpr int matrix_size = order * order;
    constexpr int vector_size = order * bands;
    const double energies[order] = {-2.0, 0.5, 2.0, 5.0};
    std::vector<double> h_full(matrix_size);
    std::vector<double> s_full(matrix_size);
    // H = L diag(energies) L^T and S = L L^T for lower bidiagonal L.
    for (int col = 0; col < order; ++col)
    {
        h_full[col + col * order] = 4.0 * energies[col];
        s_full[col + col * order] = 4.0;
        if (col > 0)
        {
            h_full[col + col * order] += 0.25 * energies[col - 1];
            h_full[col - 1 + col * order] = h_full[col + (col - 1) * order] = energies[col - 1];
            s_full[col + col * order] += 0.25;
            s_full[col - 1 + col * order] = s_full[col + (col - 1) * order] = 1.0;
        }
    }

    Parallel_Blkcyclic layout;
    layout.init(order, order, 1, MPI_COMM_WORLD, false);
    const int* desc = layout.desc();
    const int context = layout.ictxt();
    Parallel_Blkcyclic output;
    output.bind_descriptor(order, bands, 1, context);
    const std::size_t rows = layout.nrloc();
    const std::size_t cols = layout.ncloc();
    const int local_bands = output.ncloc();
    const int stored_bands = compact ? local_bands : cols;
    const std::size_t stored_elements = rows * stored_bands;
    std::vector<double> h(layout.size());
    std::vector<double> s(layout.size());
    std::vector<double> vectors(stored_elements);
    std::vector<double> values(bands);
    ModuleBase::MatrixBlock<double> hk{h.data(), rows, cols, desc};
    ModuleBase::MatrixBlock<double> sk{s.data(), rows, cols, desc};
    psi::Psi<double> psi(vectors.data(), 1, stored_bands, rows, rows, true);
    const Parallel_Blkcyclic* pool = compact ? &layout : nullptr;
    const hsolver::LcaoDiagOptions options{method, order, bands, false, pool, true};
    for (int col = 0; col < layout.ncloc(); ++col)
    {
        const int global_col = layout.icglb(col);
        for (int row = 0; row < layout.nrloc(); ++row)
        {
            const int global_row = layout.irglb(row);
            s[row + col * rows] = s_full[global_row + global_col * order];
        }
    }
    const auto original_overlap = s;

    for (int step = 0; step < 2; ++step)
    {
        SCOPED_TRACE(step);
        for (int col = 0; col < layout.ncloc(); ++col)
        {
            const int global_col = layout.icglb(col);
            for (int row = 0; row < layout.nrloc(); ++row)
            {
                const int global_row = layout.irglb(row);
                const int idx = global_row + global_col * order;
                h[row + col * rows] = h_full[idx] + step * s_full[idx];
            }
        }
        hsolver::diag_sk_cut(hk, sk, psi, values.data(), options, 0, 1e-5, nullptr);
        EXPECT_EQ(s, original_overlap);

        std::vector<double> c(vector_size);
        for (int col = 0; col < local_bands; ++col)
        {
            const int global_col = output.icglb(col);
            for (int row = 0; row < output.nrloc(); ++row)
            {
                const int global_row = output.irglb(row);
                c[global_row + global_col * order] = vectors[row + col * rows];
            }
        }
        const int count = c.size();
        Parallel_Reduce::reduce_all(c.data(), count);
        for (int band = 0; band < bands; ++band)
        {
            const double expected_energy = energies[band] + step;
            EXPECT_NEAR(values[band], expected_energy, 2e-12);
            std::vector<double> sc(order);
            for (int row = 0; row < order; ++row)
            {
                double hc = 0.0;
                for (int col = 0; col < order; ++col)
                {
                    const int idx = row + col * order;
                    const double coeff = c[col + band * order];
                    sc[row] += s_full[idx] * coeff;
                    hc += (h_full[idx] + step * s_full[idx]) * coeff;
                }
                const double expected_hc = values[band] * sc[row];
                EXPECT_NEAR(hc, expected_hc, 2e-12);
            }
            for (int other = 0; other < bands; ++other)
            {
                double gram = 0.0;
                for (int row = 0; row < order; ++row)
                {
                    gram += c[row + other * order] * sc[row];
                }
                const double expected_gram = other == band ? 1.0 : 0.0;
                EXPECT_NEAR(gram, expected_gram, 2e-12);
            }
        }
    }
}

TEST(DiagSkCut, FullRankElpaFallbackAcrossSolves)
{
    if (Parallel_Common::comm_size(MPI_COMM_WORLD) > 4)
    {
        GTEST_SKIP() << "Four-orbital test uses at most four ranks.";
    }
    check_full_rank_elpa_fallback("genelpa", false);
    check_full_rank_elpa_fallback("genelpa", true);
    check_full_rank_elpa_fallback("elpa", false);
    check_full_rank_elpa_fallback("elpa", true);
}

TEST(DiagSkCut, ElpaBackends)
{
    if (Parallel_Common::comm_size(MPI_COMM_WORLD) > 4)
    {
        GTEST_SKIP() << "Four-orbital test uses at most four ranks.";
    }
    using Complex = std::complex<double>;
    const Complex phase(0, 1);
    check_projection<double>("genelpa", true, 1.0, 3, 1);
    check_projection<Complex>("genelpa", true, phase, 3, 1);
    check_projection<double>("elpa", true, 1.0, 3, 1);
    check_projection<Complex>("elpa", true, phase, 3, 1);
}

TEST(DiagSkCut, ReducedRankSmallerThanProcessGrid)
{
    if (Parallel_Common::comm_size(MPI_COMM_WORLD) > 4)
    {
        GTEST_SKIP() << "Four-orbital test uses at most four ranks.";
    }
    using Complex = std::complex<double>;
    const Complex phase(0, 1);
    check_projection<Complex>("elpa", true, phase, 1, 1);
    check_projection<Complex>("elpa", true, phase, 2, 1);
    check_projection<Complex>("elpa", true, phase, 2, 2);
}
#endif
} // namespace

int main(int argc, char** argv)
{
    int nproc = 0;
    int nthread = 0;
    int rank = 0;
    Parallel_Global::read_pal_param(argc, argv, nproc, nthread, rank);
    POOL_WORLD = MPI_COMM_NULL;
    KP_WORLD = MPI_COMM_NULL;
    INT_BGROUP = MPI_COMM_NULL;
    BP_WORLD = MPI_COMM_NULL;
    GRID_WORLD = MPI_COMM_NULL;
    DIAG_WORLD = MPI_COMM_NULL;
    testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
    Parallel_Global::finalize_mpi();
    return result;
}
