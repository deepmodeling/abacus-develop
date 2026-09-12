#include <gtest/gtest.h>

// Keep the test-only RED handoff runnable before the opt-in helper is added.
#if __has_include("source_lcao/module_ri/sternheimer_weak_matrix_audit.h")
#include "source_lcao/module_ri/sternheimer_weak_matrix_audit.h"
#include "source_lcao/module_ri/sternheimer_abfs_perturbation.h"
#include "source_base/module_external/blas_connector.h"

#include <array>
#include <complex>
#include <limits>
#include <stdexcept>
#include <vector>

using ModuleRI::SternheimerWeakAuditShard;

TEST(SternheimerWeakMatrixAudit, DefaultsOwnEveryColumn)
{
    const auto shard = SternheimerWeakAuditShard::parse(nullptr, nullptr, 7);
    EXPECT_EQ(shard.index, 0);
    EXPECT_EQ(shard.count, 1);
    EXPECT_EQ(shard.columns(), (std::vector<int>{0, 1, 2, 3, 4, 5, 6}));
}

TEST(SternheimerWeakMatrixAudit, SixteenShardsPartitionRealAuxiliaryCounts)
{
    for (const int channels : {1113, 928})
    {
        std::vector<int> visits(channels, 0);
        for (int index = 0; index < 16; ++index)
        {
            const auto raw_index = std::to_string(index);
            const auto shard = SternheimerWeakAuditShard::parse(raw_index.c_str(), "16", channels);
            const auto owned = shard.columns();
            ASSERT_FALSE(owned.empty());
            EXPECT_EQ(owned.size(), static_cast<std::size_t>((channels - 1 - index) / 16 + 1));
            for (const int column : owned)
            {
                EXPECT_EQ(column % 16, index);
                EXPECT_TRUE(shard.owns(column));
                ++visits.at(column);
            }
            for (int column = 0; column < channels; ++column)
                EXPECT_EQ(shard.owns(column), column % 16 == index);
        }
        for (const int count : visits) EXPECT_EQ(count, 1);
    }
}

TEST(SternheimerWeakMatrixAudit, OneShardPerChannelAndUnevenLastShard)
{
    EXPECT_EQ(SternheimerWeakAuditShard::parse("4", "5", 5).columns(), (std::vector<int>{4}));
    EXPECT_EQ(SternheimerWeakAuditShard::parse("2", "3", 8).columns(), (std::vector<int>{2, 5}));
    const auto shard = SternheimerWeakAuditShard::parse("0", "1", 3);
    EXPECT_FALSE(shard.owns(-1));
    EXPECT_FALSE(shard.owns(3));
}

TEST(SternheimerWeakMatrixAudit, RejectsInvalidShardConfiguration)
{
    for (const char* index : {"-1", "2", "3", "", "x", "1x", " 1", "+1", "1.0",
                              "9999999999999999999999999999"})
        EXPECT_THROW(SternheimerWeakAuditShard::parse(index, "2", 5), std::invalid_argument);
    for (const char* count : {"0", "-1", "6", "", "x", "2x", " 2", "+2", "2.0",
                              "9999999999999999999999999999"})
        EXPECT_THROW(SternheimerWeakAuditShard::parse("0", count, 5), std::invalid_argument);
    EXPECT_THROW(SternheimerWeakAuditShard::parse(nullptr, nullptr, 0), std::invalid_argument);
    EXPECT_THROW(SternheimerWeakAuditShard::parse(nullptr, nullptr, -1), std::invalid_argument);
}

TEST(SternheimerWeakMatrixAudit, SingleChannelSamplingPreservesGlobalMetadataAndValues)
{
    ModuleRI::SternheimerRadialPerturbation s;
    s.radial_grid = {0.0, 0.5, 1.5};
    s.radial_values = {1.0, 0.7, 0.0};
    s.label = "toy_s";
    auto p = s;
    p.angular_momentum = 1;
    p.radial_index = 2;
    p.label = "toy_p";
    const std::vector<std::vector<ModuleRI::SternheimerRadialPerturbation>> radials{{s, p}, {s}};
    const std::vector<int> types{0, 1, 0};
    const std::vector<ModuleBase::Vector3<double>> positions{{0.0, 0.0, 0.0},
                                                          {0.3, 0.2, 0.1}, {0.8, 0.6, 0.4}};
    ModuleRI::SternheimerFDHamiltonian::Grid grid;
    grid.nx = 3; grid.ny = 2; grid.nz = 2;
    grid.hx = 0.7; grid.hy = 0.8; grid.hz = 0.9;
    grid.periodic = true;
    const ModuleRI::SternheimerReducedKPoint q{0.25, 0.0, 0.0};
    const auto full = ModuleRI::sample_sternheimer_abf_bloch_grid_channels(
        radials, types, positions, grid, q, -1);
    const auto description = ModuleRI::describe_sternheimer_abf_grid_channels(radials, types, positions, -1);
    ASSERT_EQ(full.size(), 9);
    ASSERT_EQ(description.size(), full.size());
    EXPECT_EQ(ModuleRI::sample_sternheimer_abf_bloch_grid_channels(
                  radials, types, positions, grid, q, 2).size(), 2);
    for (int j = 0; j < static_cast<int>(full.size()); ++j)
    {
        auto single = ModuleRI::sample_sternheimer_abf_bloch_grid_channels(
            radials, types, positions, grid, q, 1, j);
        ASSERT_EQ(single.size(), 1);
        const auto& channel = single.front();
        EXPECT_EQ(channel.channel_index, j);
        EXPECT_EQ(channel.atom_index, full[j].atom_index);
        EXPECT_EQ(channel.atom_local_index, full[j].atom_local_index);
        EXPECT_EQ(channel.type_index, full[j].type_index);
        EXPECT_EQ(channel.angular_momentum, full[j].angular_momentum);
        EXPECT_EQ(channel.magnetic_index, full[j].magnetic_index);
        EXPECT_EQ(channel.radial_index, full[j].radial_index);
        EXPECT_EQ(channel.label, full[j].label);
        EXPECT_EQ(channel.potential_r, full[j].potential_r);
        EXPECT_EQ(channel.max_abs, full[j].max_abs);
        EXPECT_EQ(description[j].channel_index, channel.channel_index);
        EXPECT_EQ(description[j].atom_local_index, channel.atom_local_index);
        EXPECT_TRUE(description[j].potential_r.empty());
        auto reference = std::vector<ModuleRI::SternheimerABFBlochGridChannel>{full[j]};
        ModuleRI::solve_sternheimer_abf_periodic_full_coulomb_in_place(single, grid, q, 0.0);
        ModuleRI::solve_sternheimer_abf_periodic_full_coulomb_in_place(reference, grid, q, 0.0);
        EXPECT_EQ(single.front().potential_r, reference.front().potential_r);
    }
    EXPECT_TRUE(ModuleRI::sample_sternheimer_abf_bloch_grid_channels(
                    radials, types, positions, grid, q, 1, 9).empty());
    EXPECT_THROW(ModuleRI::sample_sternheimer_abf_bloch_grid_channels(
                     radials, types, positions, grid, q, 1, -1), std::invalid_argument);
}

TEST(SternheimerWeakMatrixAudit, ComplexBlasAllRowsMatchExplicitContractionAndTwoSignSum)
{
    using Complex = std::complex<double>;
    using Vector = std::vector<Complex>;
    constexpr int channels = 3;
    constexpr int nf = 2;
    constexpr int nw = 4;
    // Independent logical left vertices; the packed BLAS matrices are 2x3 and 4x3.
    const std::array<Vector, channels> gf{{{{1.0, 0.5}, {-0.2, 0.8}},
                                          {{-0.7, 0.3}, {1.4, -0.6}},
                                          {{0.9, -1.2}, {0.4, 0.7}}}};
    const std::array<Vector, channels> gw{{{{0.3, 0.9}, {-0.8, 0.2}, {1.1, -0.4}, {0.6, 0.7}},
                                          {{-0.4, 0.5}, {0.2, -1.3}, {0.8, 0.6}, {-0.9, 0.4}},
                                          {{1.2, -0.8}, {-0.5, -0.1}, {0.7, 0.3}, {0.1, -0.6}}}};
    Vector packed_f(nf * channels), packed_w(nw * channels);
    for (int row = 0; row < channels; ++row)
    {
        for (int u = 0; u < nf; ++u) packed_f[u + nf * row] = gf[row][u];
        for (int u = 0; u < nw; ++u) packed_w[u + nw * row] = gw[row][u];
    }
    // Separate complex solutions, not negatives/conjugates that could hide a sum error.
    const std::array<Vector, 2> solutions{{{{0.4, 0.9}, {-0.7, 0.2}, {1.1, -0.5},
                                           {0.3, 0.8}, {-0.6, -0.4}, {0.9, 0.1}},
                                          {{-0.2, 0.6}, {0.8, -0.3}, {-0.4, 0.7},
                                           {1.2, -0.2}, {0.5, 0.9}, {-0.1, -0.8}}}};
    Vector signed_sum(channels, 0.0), explicit_sum(channels, 0.0);
    for (int branch = 0; branch < 2; ++branch)
    {
        SCOPED_TRACE(branch == 0 ? "sign +1" : "sign -1");
        const auto& coefficients = solutions[branch];
        // Nonzero initial output also checks beta=0 for F and beta=1 for W.
        Vector block_column(channels, Complex(8.0, -3.0));
        BlasConnector::gemv('C', nf, channels, Complex(0.5), packed_f.data(), nf,
                            coefficients.data(), 1, Complex(0.0), block_column.data(), 1);
        BlasConnector::gemv('C', nw, channels, Complex(0.5), packed_w.data(), nw,
                            coefficients.data() + nf, 1, Complex(1.0), block_column.data(), 1);
        for (int row = 0; row < channels; ++row)
        {
            SCOPED_TRACE(row);
            Complex expected = 0.0;
            for (int u = 0; u < nf; ++u) expected += 0.5 * std::conj(gf[row][u]) * coefficients[u];
            for (int u = 0; u < nw; ++u) expected += 0.5 * std::conj(gw[row][u]) * coefficients[nf + u];
            EXPECT_NEAR(block_column[row].real(), expected.real(), 1e-13);
            EXPECT_NEAR(block_column[row].imag(), expected.imag(), 1e-13);
            signed_sum[row] += block_column[row];
            explicit_sum[row] += expected;
        }
    }
    for (int row = 0; row < channels; ++row)
    {
        SCOPED_TRACE(row);
        EXPECT_NEAR(signed_sum[row].real(), explicit_sum[row].real(), 1e-13);
        EXPECT_NEAR(signed_sum[row].imag(), explicit_sum[row].imag(), 1e-13);
        EXPECT_GT(std::abs(explicit_sum[row].imag()), 1e-3);
    }
}

#else
TEST(SternheimerWeakMatrixAudit, OptInImplementationMustBePresent)
{
    FAIL() << "Missing sternheimer_weak_matrix_audit.h: expected test-first RED before implementation";
}
#endif
