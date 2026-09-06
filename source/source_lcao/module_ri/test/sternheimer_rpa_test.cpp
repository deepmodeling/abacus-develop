#include "source_lcao/module_ri/sternheimer_rpa.h"

#include <algorithm>
#include <atomic>
#include <complex>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <future>
#include <gtest/gtest.h>
#include <string>
#include <thread>
#include <type_traits>
#include <vector>

namespace
{

using Complex = std::complex<double>;
using Vector = ModuleRI::SternheimerRPA::Vector;

Complex dot(const Vector& lhs, const Vector& rhs)
{
    return ModuleRI::SternheimerRPA::local_grid_dot(lhs, rhs, 1.0);
}

void write_coulomb_v1_block(const std::string& filename,
                            const int pair_index,
                            const Complex value)
{
    std::ofstream out(filename.c_str(), std::ios::binary | std::ios::trunc);
    const std::int32_t header[] = {-20129433, 1, 2, 1, 2, 1, 1, 1, pair_index};
    out.write(reinterpret_cast<const char*>(header), sizeof(header));
    const std::int64_t offset = static_cast<std::int64_t>(sizeof(header) + sizeof(std::int64_t));
    out.write(reinterpret_cast<const char*>(&offset), sizeof(offset));
    out.write(reinterpret_cast<const char*>(&value), sizeof(value));
}

} // namespace

TEST(SternheimerRPA, BuildRhsFromHartreePerturbation)
{
    const std::vector<double> v_r = {2.0, -1.0};
    const Vector psi_r = {Complex(1.0, 1.0), Complex(2.0, 0.0)};

    Vector rhs;
    ModuleRI::SternheimerRPA::build_rhs_from_hartree_perturbation(v_r, psi_r, rhs);

    ASSERT_EQ(rhs.size(), 2);
    EXPECT_DOUBLE_EQ(rhs[0].real(), -2.0);
    EXPECT_DOUBLE_EQ(rhs[0].imag(), -2.0);
    EXPECT_DOUBLE_EQ(rhs[1].real(), 2.0);
    EXPECT_DOUBLE_EQ(rhs[1].imag(), 0.0);
}

TEST(SternheimerRPA, ComplexQPotentialBuildsRhsAndUsesMinusQProbe)
{
    const Vector v_q = {Complex(2.0, 1.0), Complex(-1.0, 0.5)};
    const Vector psi_k = {Complex(1.0, -1.0), Complex(0.5, 2.0)};
    const Vector delta_psi_kq = {Complex(0.25, 0.5), Complex(-1.0, 0.75)};
    constexpr double grid_weight = 0.2;

    Vector rhs;
    ModuleRI::SternheimerRPA::build_rhs_from_hartree_perturbation(v_q, psi_k, rhs);
    ASSERT_EQ(rhs.size(), psi_k.size());
    for (std::size_t ir = 0; ir != rhs.size(); ++ir)
    {
        EXPECT_NEAR(rhs[ir].real(), (-v_q[ir] * psi_k[ir]).real(), 1.0e-14);
        EXPECT_NEAR(rhs[ir].imag(), (-v_q[ir] * psi_k[ir]).imag(), 1.0e-14);
    }

    const Complex value = ModuleRI::SternheimerRPA::accumulate_polarizability_grid_element(
        v_q, psi_k, delta_psi_kq, grid_weight);
    const Complex expected
        = grid_weight
          * (std::conj(psi_k[0]) * std::conj(v_q[0]) * delta_psi_kq[0]
             + std::conj(psi_k[1]) * std::conj(v_q[1]) * delta_psi_kq[1]);
    EXPECT_NEAR(value.real(), expected.real(), 1.0e-14);
    EXPECT_NEAR(value.imag(), expected.imag(), 1.0e-14);
}

TEST(SternheimerRPA, AccumulatePolarizabilityGridElement)
{
    const std::vector<double> v_r = {2.0, -1.0};
    const Vector psi_r = {Complex(1.0, 1.0), Complex(2.0, 0.0)};
    const Vector delta_psi_r = {Complex(0.5, 0.0), Complex(0.0, 1.0)};
    constexpr double grid_weight = 0.25;

    const Complex value
        = ModuleRI::SternheimerRPA::accumulate_polarizability_grid_element(v_r, psi_r, delta_psi_r, grid_weight);
    const Complex expected
        = grid_weight * (std::conj(psi_r[0]) * v_r[0] * delta_psi_r[0] + std::conj(psi_r[1]) * v_r[1] * delta_psi_r[1]);

    EXPECT_NEAR(value.real(), expected.real(), 1.0e-14);
    EXPECT_NEAR(value.imag(), expected.imag(), 1.0e-14);
}

TEST(SternheimerRPA, AccumulateChi0BranchColumnIncludesOccupation)
{
    const std::vector<std::vector<double>> potentials = {{2.0, -1.0}, {0.5, 3.0}};
    const Vector psi_r = {Complex(1.0, 1.0), Complex(2.0, 0.0)};
    const Vector delta_psi_r = {Complex(0.5, 0.0), Complex(0.0, 1.0)};
    constexpr double grid_weight = 0.25;
    constexpr double occupation = 2.0;

    std::vector<Complex> branch(4, Complex(0.0, 0.0));
    ModuleRI::SternheimerRPA::accumulate_chi0_branch_column(potentials,
                                                            psi_r,
                                                            delta_psi_r,
                                                            grid_weight,
                                                            occupation,
                                                            1,
                                                            branch);

    const Complex expected0 = occupation
                              * ModuleRI::SternheimerRPA::accumulate_polarizability_grid_element(potentials[0],
                                                                                                 psi_r,
                                                                                                 delta_psi_r,
                                                                                                 grid_weight);
    const Complex expected1 = occupation
                              * ModuleRI::SternheimerRPA::accumulate_polarizability_grid_element(potentials[1],
                                                                                                 psi_r,
                                                                                                 delta_psi_r,
                                                                                                 grid_weight);
    EXPECT_EQ(branch[0], Complex(0.0, 0.0));
    EXPECT_EQ(branch[2], Complex(0.0, 0.0));
    EXPECT_NEAR(branch[1].real(), expected0.real(), 1.0e-14);
    EXPECT_NEAR(branch[1].imag(), expected0.imag(), 1.0e-14);
    EXPECT_NEAR(branch[3].real(), expected1.real(), 1.0e-14);
    EXPECT_NEAR(branch[3].imag(), expected1.imag(), 1.0e-14);
}

TEST(SternheimerRPA, CollinearSpinBranchesAddBeforeConjugateSymmetrization)
{
    const std::vector<std::vector<double>> potentials = {{2.0, -1.0}, {0.5, 3.0}};
    const Vector psi_up = {Complex(1.0, 0.5), Complex(0.5, -0.25)};
    const Vector delta_up = {Complex(-0.2, 0.3), Complex(0.4, 0.1)};
    const Vector psi_down = {Complex(0.3, -0.4), Complex(-0.8, 0.2)};
    const Vector delta_down = {Complex(0.6, 0.2), Complex(-0.1, 0.5)};
    constexpr double grid_weight = 0.125;

    std::vector<Complex> spin_up(4, Complex(0.0, 0.0));
    std::vector<Complex> spin_down(4, Complex(0.0, 0.0));
    std::vector<Complex> combined(4, Complex(0.0, 0.0));
    ModuleRI::SternheimerRPA::accumulate_chi0_branch_column(potentials, psi_up, delta_up, grid_weight, 1.0, 0, spin_up);
    ModuleRI::SternheimerRPA::accumulate_chi0_branch_column(potentials,
                                                            psi_down,
                                                            delta_down,
                                                            grid_weight,
                                                            1.0,
                                                            0,
                                                            spin_down);
    ModuleRI::SternheimerRPA::accumulate_chi0_branch_column(potentials,
                                                            psi_up,
                                                            delta_up,
                                                            grid_weight,
                                                            1.0,
                                                            0,
                                                            combined);
    ModuleRI::SternheimerRPA::accumulate_chi0_branch_column(potentials,
                                                            psi_down,
                                                            delta_down,
                                                            grid_weight,
                                                            1.0,
                                                            0,
                                                            combined);

    for (std::size_t index = 0; index != combined.size(); ++index)
    {
        EXPECT_NEAR(combined[index].real(), (spin_up[index] + spin_down[index]).real(), 1.0e-14);
        EXPECT_NEAR(combined[index].imag(), (spin_up[index] + spin_down[index]).imag(), 1.0e-14);
    }

    const auto symmetrized_combined = ModuleRI::SternheimerRPA::symmetrize_chi0_imaginary_frequency(combined, 2);
    const auto symmetrized_up = ModuleRI::SternheimerRPA::symmetrize_chi0_imaginary_frequency(spin_up, 2);
    const auto symmetrized_down = ModuleRI::SternheimerRPA::symmetrize_chi0_imaginary_frequency(spin_down, 2);
    for (std::size_t index = 0; index != symmetrized_combined.size(); ++index)
    {
        EXPECT_NEAR(symmetrized_combined[index].real(),
                    (symmetrized_up[index] + symmetrized_down[index]).real(),
                    1.0e-14);
        EXPECT_NEAR(symmetrized_combined[index].imag(),
                    (symmetrized_up[index] + symmetrized_down[index]).imag(),
                    1.0e-14);
    }
}

TEST(SternheimerRPA, AccumulateComplexQChi0BranchColumnReadsChannelStorageDirectly)
{
    std::vector<ModuleRI::SternheimerABFBlochGridChannel> channels(2);
    channels[0].potential_r = {Complex(2.0, 1.0), Complex(-1.0, 0.5)};
    channels[1].potential_r = {Complex(0.5, -0.25), Complex(3.0, 2.0)};
    const Vector psi_k = {Complex(1.0, 1.0), Complex(2.0, -0.5)};
    const Vector delta_psi_kq = {Complex(0.5, -0.25), Complex(-0.75, 1.0)};
    constexpr double grid_weight = 0.25;
    constexpr double occupation = 1.5;

    std::vector<Complex> direct(4, Complex(0.0, 0.0));
    ModuleRI::SternheimerRPA::accumulate_chi0_branch_column(channels,
                                                            psi_k,
                                                            delta_psi_kq,
                                                            grid_weight,
                                                            occupation,
                                                            0,
                                                            direct);

    const std::vector<Vector> copied = {channels[0].potential_r, channels[1].potential_r};
    std::vector<Complex> expected(4, Complex(0.0, 0.0));
    ModuleRI::SternheimerRPA::accumulate_chi0_branch_column(copied,
                                                            psi_k,
                                                            delta_psi_kq,
                                                            grid_weight,
                                                            occupation,
                                                            0,
                                                            expected);
    EXPECT_EQ(direct, expected);
}

TEST(SternheimerRPA, AccumulateComplexQChi0BranchColumnsMatchScalarColumns)
{
    std::vector<ModuleRI::SternheimerABFBlochGridChannel> channels(3);
    channels[0].potential_r = {Complex(2.0, 1.0), Complex(-1.0, 0.5)};
    channels[1].potential_r = {Complex(0.5, -0.25), Complex(3.0, 2.0)};
    channels[2].potential_r = {Complex(-0.75, 0.25), Complex(1.5, -2.0)};
    const Vector psi_k = {Complex(1.0, 1.0), Complex(2.0, -0.5)};
    const std::vector<Vector> delta_psi_batch{{Complex(0.5, -0.25), Complex(-0.75, 1.0)},
                                              {Complex(-0.25, 0.75), Complex(1.25, 0.5)}};
    const std::vector<const Vector*> delta_psi_batch_views{&delta_psi_batch[0], &delta_psi_batch[1]};
    constexpr double grid_weight = 0.25;
    constexpr double occupation = 1.5;

    std::vector<Complex> batched(9, Complex(0.0, 0.0));
    ModuleRI::SternheimerRPA::accumulate_chi0_branch_columns(channels,
                                                             psi_k,
                                                             delta_psi_batch_views,
                                                             grid_weight,
                                                             occupation,
                                                             1,
                                                             batched);

    std::vector<Complex> scalar(9, Complex(0.0, 0.0));
    for (std::size_t offset = 0; offset != delta_psi_batch.size(); ++offset)
    {
        ModuleRI::SternheimerRPA::accumulate_chi0_branch_column(channels,
                                                                psi_k,
                                                                delta_psi_batch[offset],
                                                                grid_weight,
                                                                occupation,
                                                                static_cast<int>(offset) + 1,
                                                                scalar);
    }
    for (std::size_t index = 0; index != batched.size(); ++index)
    {
        EXPECT_NEAR(batched[index].real(), scalar[index].real(), 1.0e-14);
        EXPECT_NEAR(batched[index].imag(), scalar[index].imag(), 1.0e-14);
    }
}

TEST(SternheimerRPA, SymmetrizeImaginaryFrequencyChi0AddsAdjointBranch)
{
    const std::vector<Complex> branch = {Complex(1.0, 2.0), Complex(3.0, -4.0), Complex(5.0, 6.0), Complex(-7.0, 8.0)};

    const std::vector<Complex> chi0 = ModuleRI::SternheimerRPA::symmetrize_chi0_imaginary_frequency(branch, 2);

    ASSERT_EQ(chi0.size(), 4);
    EXPECT_EQ(chi0[0], branch[0] + std::conj(branch[0]));
    EXPECT_EQ(chi0[1], branch[1] + std::conj(branch[2]));
    EXPECT_EQ(chi0[2], branch[2] + std::conj(branch[1]));
    EXPECT_EQ(chi0[3], branch[3] + std::conj(branch[3]));
    EXPECT_NE(chi0[1], Complex(2.0 * branch[1].real(), 0.0));
}

TEST(SternheimerRPA, WriteChi0V1FileUsesAtomPairBlocks)
{
    ModuleRI::SternheimerRPA::Chi0V1Metadata metadata;
    metadata.iq = 1;
    metadata.ifrequency = 2;
    metadata.omega = 0.5;
    metadata.weight = 0.125;
    metadata.atom_naux = {1, 1};

    const std::vector<ModuleRI::SternheimerRPA::AuxiliaryChannel> channels = {{0, 0, 0}, {1, 1, 0}};
    const std::vector<Complex> chi0 = {Complex(1.0, 0.0), Complex(2.0, -0.5), Complex(2.0, 0.5), Complex(3.0, 0.0)};

    const std::string filename = ::testing::TempDir() + "/sternheimer_chi0_v1_test.dat";
    ModuleRI::SternheimerRPA::write_chi0_v1_file(filename, metadata, channels, chi0);

    std::ifstream in(filename.c_str(), std::ios::binary);
    ASSERT_TRUE(in.good());

    auto read_i32 = [&in]() {
        std::int32_t value = 0;
        in.read(reinterpret_cast<char*>(&value), sizeof(value));
        return value;
    };
    auto read_i64 = [&in]() {
        std::int64_t value = 0;
        in.read(reinterpret_cast<char*>(&value), sizeof(value));
        return value;
    };
    auto read_double = [&in]() {
        double value = 0.0;
        in.read(reinterpret_cast<char*>(&value), sizeof(value));
        return value;
    };

    EXPECT_EQ(read_i32(), ModuleRI::SternheimerRPA::chi0_v1_marker());
    EXPECT_EQ(read_i32(), 1);
    EXPECT_EQ(read_i32(), 2);
    EXPECT_EQ(read_i32(), 2);
    EXPECT_EQ(read_i32(), 1);
    EXPECT_EQ(read_i32(), 2);
    EXPECT_DOUBLE_EQ(read_double(), 0.5);
    EXPECT_DOUBLE_EQ(read_double(), 0.125);
    EXPECT_EQ(read_i32(), 3);
    EXPECT_EQ(read_i32(), 1);
    EXPECT_EQ(read_i32(), 1);
    EXPECT_EQ(read_i32(), 0);
    const std::int64_t first_offset = read_i64();
    EXPECT_EQ(read_i32(), 1);
    const std::int64_t second_offset = read_i64();
    EXPECT_EQ(read_i32(), 2);
    const std::int64_t third_offset = read_i64();
    EXPECT_LT(first_offset, second_offset);
    EXPECT_LT(second_offset, third_offset);

    std::remove(filename.c_str());
}

TEST(SternheimerRPA, ReadsCoulombV1BlocksAcrossRankFiles)
{
    const std::string prefix = ::testing::TempDir() + "/sternheimer_coulomb_v1_rank";
    const std::vector<Complex> blocks{Complex(1.0, 0.0), Complex(2.0, 3.0), Complex(4.0, 0.0)};
    std::vector<std::string> filenames;
    for (int pair_index = 0; pair_index != 3; ++pair_index)
    {
        const std::string filename = prefix + std::to_string(pair_index) + ".dat";
        write_coulomb_v1_block(filename, pair_index, blocks[static_cast<std::size_t>(pair_index)]);
        filenames.push_back(filename);
    }

    const auto matrix = ModuleRI::SternheimerRPA::read_coulomb_v1_files(filenames);

    EXPECT_EQ(matrix.iq, 1);
    EXPECT_EQ(matrix.atom_naux, std::vector<int>({1, 1}));
    EXPECT_EQ(matrix.values,
              std::vector<Complex>({Complex(1.0, 0.0),
                                    Complex(2.0, 3.0),
                                    Complex(2.0, -3.0),
                                    Complex(4.0, 0.0)}));
    for (const std::string& filename: filenames)
    {
        std::remove(filename.c_str());
    }
}

TEST(SternheimerRPA, ProjectOutSubspace)
{
    const std::vector<Vector> occupied
        = {{Complex(1.0, 0.0), Complex(0.0, 0.0)}, {Complex(0.0, 0.0), Complex(0.0, 1.0)}};
    Vector vec = {Complex(2.0, 1.0), Complex(3.0, -4.0)};

    ModuleRI::SternheimerRPA::project_out_subspace(occupied, dot, vec);

    EXPECT_NEAR(std::abs(dot(occupied[0], vec)), 0.0, 1.0e-14);
    EXPECT_NEAR(std::abs(dot(occupied[1], vec)), 0.0, 1.0e-14);
}

TEST(SternheimerRPA, SubspaceProjectorAcceptsNonowningPointerViews)
{
    using Projector = ModuleRI::SternheimerSubspaceProjector;
    using Views = std::vector<const Vector*>;
    EXPECT_TRUE((std::is_constructible<Projector, Views, Projector::Dot>::value));
    EXPECT_TRUE((std::is_constructible<Projector, Views, double>::value));
    EXPECT_TRUE((std::is_constructible<Projector, Views, double, bool>::value));
}

TEST(SternheimerRPA, SubspaceProjectorRejectsMismatchedBasisSizesAtConstruction)
{
    const std::vector<Vector> invalid{{Complex(1.0, 0.0)}, {Complex(1.0, 0.0), Complex(0.0, 0.0)}};
    EXPECT_THROW(ModuleRI::SternheimerSubspaceProjector(invalid, dot), std::invalid_argument);
    EXPECT_THROW(ModuleRI::SternheimerSubspaceProjector(invalid, 1.0), std::invalid_argument);
}

TEST(SternheimerRPA, SubspaceProjectorBorrowsOriginalVectorsAndOwnsOnlyViewList)
{
    using Projector = ModuleRI::SternheimerSubspaceProjector;
    const Vector first{Complex(1.0, 0.0), Complex(0.0, 1.0), Complex(0.0, 0.0)};
    const Vector second{Complex(0.0, 1.0), Complex(1.0, 0.0), Complex(0.0, 0.0)};
    int dot_calls = 0;
    const auto identity_dot = [&](const Vector& lhs, const Vector& rhs) {
        EXPECT_TRUE(&lhs == &first || &lhs == &second);
        ++dot_calls;
        return dot(lhs, rhs);
    };
    const Projector projector = [&]() {
        std::vector<const Vector*> views{&first, &second};
        return Projector(views, identity_dot);
    }();
    EXPECT_EQ(dot_calls, 2);
    Vector scalar{Complex(2.0, 1.0), Complex(-3.0, 0.5), Complex(4.0, -2.0)};
    std::vector<Vector> batch{scalar};
    projector.project(scalar);
    projector.project_batch(batch);
    EXPECT_EQ(dot_calls, 6);
    EXPECT_EQ(batch.front(), scalar);
    EXPECT_NEAR(std::abs(scalar[0]), 0.0, 1.0e-14);
    EXPECT_NEAR(std::abs(scalar[1]), 0.0, 1.0e-14);
    EXPECT_EQ(scalar[2], Complex(4.0, -2.0));
}

TEST(SternheimerRPA, SubspaceProjectorLegacyConstructorAlsoBorrowsOriginalVectors)
{
    const std::vector<Vector> basis{{Complex(1.0, 0.0), Complex(0.0, 0.0)}};
    const ModuleRI::SternheimerSubspaceProjector projector(
        basis,
        [&](const Vector& lhs, const Vector& rhs) {
            EXPECT_EQ(&lhs, &basis.front());
            return dot(lhs, rhs);
        });
    Vector value{Complex(2.0, -1.0), Complex(0.5, 0.25)};
    projector.project(value);
    EXPECT_NEAR(std::abs(value.front()), 0.0, 1.0e-14);
}

TEST(SternheimerRPA, SubspaceProjectorPointerViewsRejectNullAndInconsistentSizes)
{
    using Projector = ModuleRI::SternheimerSubspaceProjector;
    const Vector first(3, Complex(1.0, 0.0));
    const Vector shorter(2, Complex(1.0, 0.0));
    const std::vector<std::vector<const Vector*>> invalid_views{
        {nullptr}, {&first, nullptr}, {&first, &shorter}};
    for (const auto& views: invalid_views)
    {
        EXPECT_THROW(Projector(views, dot), std::invalid_argument);
        EXPECT_THROW(Projector(views, 1.0), std::invalid_argument);
        EXPECT_THROW(Projector(views, 1.0, true), std::invalid_argument);
    }
    const std::vector<const Vector*> views{&first};
    EXPECT_THROW(Projector(views, Projector::Dot{}), std::invalid_argument);
    EXPECT_THROW(Projector(views, 0.0), std::invalid_argument);
    for (const bool packed: {false, true})
    {
        const Projector projector(views, 1.0, packed);
        Vector invalid(2);
        EXPECT_THROW(projector.project(invalid), std::invalid_argument);
        std::vector<Vector> invalid_batch{first, invalid};
        EXPECT_THROW(projector.project_batch(invalid_batch), std::invalid_argument);
    }
}

TEST(SternheimerRPA, SubspaceProjectorPointerViewsAllowEmptyAndZeroNormBases)
{
    using Projector = ModuleRI::SternheimerSubspaceProjector;
    const Vector zero(3, Complex(0.0, 0.0));
    const Vector input{Complex(1.0, 2.0), Complex(3.0, 4.0), Complex(-1.0, 0.5)};
    for (const auto& views: {std::vector<const Vector*>{}, std::vector<const Vector*>{&zero}})
    {
        for (const bool packed: {false, true})
        {
            const Projector projector(views, 0.5, packed);
            Vector scalar = input;
            std::vector<Vector> batch{input, input};
            projector.project(scalar);
            projector.project_batch(batch);
            EXPECT_EQ(scalar, input);
            EXPECT_EQ(batch, (std::vector<Vector>{input, input}));
            std::vector<Vector> empty;
            projector.project_batch(empty);
            EXPECT_TRUE(empty.empty());
        }
    }
}

TEST(SternheimerRPA, SubspaceProjectorOwnsOnePackedBasisAndNoScalarGridCopy)
{
    using Projector = ModuleRI::SternheimerSubspaceProjector;
    constexpr std::size_t grid_size = 1024;
    Vector first(grid_size, Complex(0.0, 0.0));
    Vector second = first;
    first[0] = Complex(1.0, 0.0);
    second[1] = Complex(0.0, 1.0);
    const std::vector<const Vector*> views{&first, &second};
    const Projector scalar(views, 1.0);
    const Projector packed(views, 1.0, true);
    EXPECT_EQ(scalar.storage_bytes(), views.size() * (sizeof(const Vector*) + sizeof(Complex)));
    EXPECT_EQ(packed.storage_bytes() - scalar.storage_bytes(),
              views.size() * grid_size * sizeof(Complex));
    const auto bytes_before = packed.storage_bytes();
    for (int repeat = 0; repeat != 8; ++repeat)
    {
        std::vector<Vector> batch{first, second};
        packed.project_batch(batch);
        EXPECT_EQ(packed.storage_bytes(), bytes_before);
    }
}

TEST(SternheimerRPA, SubspaceProjectorReportsBasisSizeIncludingZeroNormViews)
{
    using Projector = ModuleRI::SternheimerSubspaceProjector;
    const Vector zero(3, Complex(0.0, 0.0));
    const Vector first(3, Complex(1.0, 0.0));
    const Projector empty(std::vector<const Vector*>{}, 1.0);
    const Projector views(std::vector<const Vector*>{&zero, &first}, 1.0, true);
    const std::vector<Vector> basis{first};
    const Projector legacy(basis, dot);
    EXPECT_EQ(empty.basis_size(), 0U);
    EXPECT_EQ(views.basis_size(), 2U);
    EXPECT_EQ(legacy.basis_size(), 1U);
}

TEST(SternheimerRPA, SubspaceProjectorPointerViewsSupportConcurrentScalarAndBatchReuse)
{
    using Projector = ModuleRI::SternheimerSubspaceProjector;
    const Vector first{Complex(1.0, 0.0), Complex(0.0, 1.0), Complex(0.0, 0.0), Complex(0.0, 0.0)};
    const Vector second{Complex(0.0, 1.0), Complex(1.0, 0.0), Complex(0.0, 0.0), Complex(0.0, 0.0)};
    const Vector zero(first.size(), Complex(0.0, 0.0));
    const std::vector<const Vector*> views{&first, &second, &zero};
    const std::vector<Vector> basis{first, second, zero};
    const std::vector<Vector> inputs{
        {Complex(2.0, 1.0), Complex(3.0, -4.0), Complex(0.5, 0.25), Complex(-0.2, 0.7)},
        {Complex(-1.0, 2.0), Complex(0.5, 0.25), Complex(3.0, -0.75), Complex(0.8, -0.6)},
        {Complex(0.1, -0.2), Complex(-0.3, 0.4), Complex(0.5, -0.6), Complex(0.9, 1.1)}};
    auto expected = inputs;
    for (auto& column: expected)
    {
        ModuleRI::SternheimerRPA::project_out_subspace(basis, dot, column);
    }
    for (int mode = 0; mode != 3; ++mode)
    {
        const Projector projector = mode == 0 ? Projector(views, dot) : Projector(views, 0.25, mode == 2);
        const auto bytes_before = projector.storage_bytes();
        std::atomic<int> ready{0};
        std::atomic<bool> start{false};
        std::vector<std::future<double>> workers;
        for (int worker = 0; worker != 8; ++worker)
        {
            workers.push_back(std::async(std::launch::async, [&]() {
                ++ready;
                while (!start.load())
                {
                    std::this_thread::yield();
                }
                double max_error = 0.0;
                for (int repeat = 0; repeat != 12; ++repeat)
                {
                    auto scalar = inputs;
                    auto batch = inputs;
                    for (auto& column: scalar)
                    {
                        projector.project(column);
                    }
                    projector.project_batch(batch);
                    for (std::size_t column = 0; column != inputs.size(); ++column)
                    {
                        for (std::size_t ir = 0; ir != first.size(); ++ir)
                        {
                            max_error = std::max(max_error, std::abs(scalar[column][ir] - expected[column][ir]));
                            max_error = std::max(max_error, std::abs(batch[column][ir] - expected[column][ir]));
                        }
                    }
                }
                return max_error;
            }));
        }
        while (ready.load() != 8)
        {
            std::this_thread::yield();
        }
        start = true;
        for (auto& worker: workers)
        {
            EXPECT_LE(worker.get(), 1.0e-13);
        }
        EXPECT_EQ(projector.storage_bytes(), bytes_before);
    }
}

TEST(SternheimerRPA, CachedSubspaceProjectorReusesBasisNorms)
{
    int dot_calls = 0;
    const auto counting_dot = [&dot_calls](const Vector& lhs, const Vector& rhs) {
        ++dot_calls;
        return dot(lhs, rhs);
    };
    const std::vector<Vector> occupied
        = {{Complex(2.0, 0.0), Complex(0.0, 0.0)}, {Complex(0.0, 0.0), Complex(0.0, 3.0)}};
    const ModuleRI::SternheimerSubspaceProjector projector(occupied, counting_dot);

    EXPECT_EQ(dot_calls, 2);
    Vector first = {Complex(2.0, 1.0), Complex(3.0, -4.0)};
    projector.project(first);
    EXPECT_EQ(dot_calls, 4);
    Vector second = {Complex(-1.0, 2.0), Complex(0.5, 0.25)};
    projector.project(second);
    EXPECT_EQ(dot_calls, 6);
    EXPECT_NEAR(std::abs(dot(occupied[0], first)), 0.0, 1.0e-14);
    EXPECT_NEAR(std::abs(dot(occupied[1], first)), 0.0, 1.0e-14);
    EXPECT_NEAR(std::abs(dot(occupied[0], second)), 0.0, 1.0e-14);
    EXPECT_NEAR(std::abs(dot(occupied[1], second)), 0.0, 1.0e-14);
}

TEST(SternheimerRPA, CachedSubspaceProjectorBatchMatchesScalar)
{
    const std::vector<Vector> occupied
        = {{Complex(2.0, 0.0), Complex(0.0, 0.0), Complex(1.0, -0.5)},
           {Complex(0.0, 0.0), Complex(0.0, 3.0), Complex(-0.25, 0.75)}};
    const ModuleRI::SternheimerSubspaceProjector projector(occupied, dot);
    std::vector<Vector> expected
        = {{Complex(2.0, 1.0), Complex(3.0, -4.0), Complex(0.5, 0.25)},
           {Complex(-1.0, 2.0), Complex(0.5, 0.25), Complex(3.0, -0.75)},
           {Complex(0.1, -0.2), Complex(-0.3, 0.4), Complex(0.5, -0.6)},
           {Complex(-0.7, 0.8), Complex(0.9, -1.0), Complex(-1.1, 1.2)}};
    std::vector<Vector> actual = expected;
    for (Vector& column: expected)
    {
        projector.project(column);
    }

    projector.project_batch(actual);

    ASSERT_EQ(actual.size(), expected.size());
    for (std::size_t column = 0; column != expected.size(); ++column)
    {
        ASSERT_EQ(actual[column].size(), expected[column].size());
        for (std::size_t ir = 0; ir != expected[column].size(); ++ir)
        {
            EXPECT_NEAR(actual[column][ir].real(), expected[column][ir].real(), 1.0e-13);
            EXPECT_NEAR(actual[column][ir].imag(), expected[column][ir].imag(), 1.0e-13);
        }
    }

    std::vector<Vector> empty;
    projector.project_batch(empty);
    EXPECT_TRUE(empty.empty());

    std::vector<Vector> one{{Complex(2.0, 1.0), Complex(3.0, -4.0), Complex(0.5, 0.25)}};
    projector.project_batch(one);
    ASSERT_EQ(one.size(), 1U);
    EXPECT_EQ(one.front(), expected.front());

    std::vector<Vector> invalid{{Complex(1.0, 0.0)}};
    EXPECT_THROW(projector.project_batch(invalid), std::invalid_argument);
}

TEST(SternheimerRPA, WeightedSubspaceProjectorBatchMatchesScalar)
{
    const double grid_weight = 0.125;
    const std::vector<Vector> occupied
        = {{Complex(2.0, 0.0), Complex(0.0, 0.0), Complex(1.0, -0.5)},
           {Complex(0.0, 0.0), Complex(0.0, 3.0), Complex(-0.25, 0.75)}};
    const ModuleRI::SternheimerSubspaceProjector scalar_projector(
        occupied,
        [grid_weight](const Vector& lhs, const Vector& rhs) { return grid_weight * dot(lhs, rhs); });
    const ModuleRI::SternheimerSubspaceProjector batch_projector(occupied, grid_weight);
    std::vector<Vector> expected
        = {{Complex(2.0, 1.0), Complex(3.0, -4.0), Complex(0.5, 0.25)},
           {Complex(-1.0, 2.0), Complex(0.5, 0.25), Complex(3.0, -0.75)},
           {Complex(0.1, -0.2), Complex(-0.3, 0.4), Complex(0.5, -0.6)},
           {Complex(-0.7, 0.8), Complex(0.9, -1.0), Complex(-1.1, 1.2)}};
    std::vector<Vector> actual = expected;
    for (Vector& column: expected)
    {
        scalar_projector.project(column);
    }

    batch_projector.project_batch(actual);

    ASSERT_EQ(actual.size(), expected.size());
    for (std::size_t column = 0; column != expected.size(); ++column)
    {
        ASSERT_EQ(actual[column].size(), expected[column].size());
        for (std::size_t ir = 0; ir != expected[column].size(); ++ir)
        {
            EXPECT_NEAR(actual[column][ir].real(), expected[column][ir].real(), 1.0e-13);
            EXPECT_NEAR(actual[column][ir].imag(), expected[column][ir].imag(), 1.0e-13);
        }
    }

    EXPECT_THROW(ModuleRI::SternheimerSubspaceProjector(occupied, 0.0), std::invalid_argument);
}

TEST(SternheimerRPA, OrthonormalSubspaceProjectorBatchMatchesScalar)
{
    const double grid_weight = 0.25;
    const std::vector<Vector> occupied
        = {{Complex(2.0, 0.0), Complex(0.0, 0.0), Complex(0.0, 0.0), Complex(0.0, 0.0)},
           {Complex(0.0, 0.0), Complex(0.0, 2.0), Complex(0.0, 0.0), Complex(0.0, 0.0)}};
    const ModuleRI::SternheimerSubspaceProjector scalar_projector(
        occupied,
        [grid_weight](const Vector& lhs, const Vector& rhs) { return grid_weight * dot(lhs, rhs); });
    const ModuleRI::SternheimerSubspaceProjector batch_projector(occupied, grid_weight, true);
    std::vector<Vector> expected
        = {{Complex(2.0, 1.0), Complex(3.0, -4.0), Complex(0.5, 0.25), Complex(-0.2, 0.7)},
           {Complex(-1.0, 2.0), Complex(0.5, 0.25), Complex(3.0, -0.75), Complex(0.8, -0.6)},
           {Complex(0.1, -0.2), Complex(-0.3, 0.4), Complex(0.5, -0.6), Complex(0.9, 1.1)}};
    std::vector<Vector> actual = expected;
    for (Vector& column: expected)
    {
        scalar_projector.project(column);
    }

    batch_projector.project_batch(actual);

    ASSERT_EQ(actual.size(), expected.size());
    for (std::size_t column = 0; column != expected.size(); ++column)
    {
        ASSERT_EQ(actual[column].size(), expected[column].size());
        for (std::size_t ir = 0; ir != expected[column].size(); ++ir)
        {
            EXPECT_NEAR(actual[column][ir].real(), expected[column][ir].real(), 1.0e-13);
            EXPECT_NEAR(actual[column][ir].imag(), expected[column][ir].imag(), 1.0e-13);
        }
    }
}

TEST(SternheimerRPA, SpectralPreconditionerIsDefault)
{
    const ModuleRI::SternheimerRPA::SolverOptions options{};

    EXPECT_TRUE(options.use_fd_spectral_preconditioner);
    EXPECT_DOUBLE_EQ(options.fd_spectral_preconditioner_regularization, 0.0);
}

TEST(SternheimerRPA, ApplyKineticPreconditioner)
{
    const std::vector<double> kinetic = {1.0, 3.0};
    const Vector input = {Complex(2.0, 0.0), Complex(0.0, 4.0)};
    Vector output;

    ModuleRI::SternheimerRPA::apply_kinetic_preconditioner(kinetic, 0.5, 2.0, 0.1, input, output);

    ASSERT_EQ(output.size(), input.size());
    EXPECT_NEAR(output[0].real(), (input[0] / Complex(0.6, 2.0)).real(), 1.0e-14);
    EXPECT_NEAR(output[0].imag(), (input[0] / Complex(0.6, 2.0)).imag(), 1.0e-14);
    EXPECT_NEAR(output[1].real(), (input[1] / Complex(2.6, 2.0)).real(), 1.0e-14);
    EXPECT_NEAR(output[1].imag(), (input[1] / Complex(2.6, 2.0)).imag(), 1.0e-14);
}

TEST(SternheimerRPA, TransitionEnergyWindowUsesOccupiedToEmptyPairsAndReturnsHartree)
{
    const std::vector<double> eigenvalues_ry = {-1.0, -0.5, 0.25, 1.5};
    const std::vector<double> occupations = {2.0, 2.0, 0.0, 0.0};

    const auto window
        = ModuleRI::SternheimerRPA::transition_energy_window_from_eigenvalues_ry(eigenvalues_ry, occupations);

    EXPECT_DOUBLE_EQ(window.emin_ha, 0.375);
    EXPECT_DOUBLE_EQ(window.emax_ha, 1.25);
}

TEST(SternheimerRPA, OptionalTransitionEnergyWindowDoesNotBlockExternalFrequencyGrid)
{
    using Window = ModuleRI::SternheimerRPA::TransitionEnergyWindow;
    Window window;

    EXPECT_FALSE(ModuleRI::SternheimerRPA::try_transition_energy_window_from_eigenvalues_ry(
        {-1.0, -0.5}, {2.0, 2.0}, window));
    EXPECT_THROW(ModuleRI::SternheimerRPA::try_transition_energy_window_from_eigenvalues_ry(
                     {-1.0}, {2.0, 0.0}, window),
                 std::invalid_argument);

    EXPECT_TRUE(ModuleRI::SternheimerRPA::try_transition_energy_window_from_eigenvalues_ry(
        {-1.0, 0.5}, {2.0, 0.0}, window));
    EXPECT_DOUBLE_EQ(window.emin_ha, 0.75);
    EXPECT_DOUBLE_EQ(window.emax_ha, 0.75);
}

TEST(SternheimerRPA, MergeTransitionEnergyWindowsCoversAllSpinChannels)
{
    using Window = ModuleRI::SternheimerRPA::TransitionEnergyWindow;
    const auto merged = ModuleRI::SternheimerRPA::merge_transition_energy_windows({Window{0.4, 2.0}, Window{0.2, 3.5}});

    EXPECT_DOUBLE_EQ(merged.emin_ha, 0.2);
    EXPECT_DOUBLE_EQ(merged.emax_ha, 3.5);
    EXPECT_THROW(ModuleRI::SternheimerRPA::merge_transition_energy_windows({}), std::invalid_argument);
    EXPECT_THROW(ModuleRI::SternheimerRPA::merge_transition_energy_windows({Window{0.0, 1.0}}), std::invalid_argument);
}

TEST(SternheimerRPA, GreenXMinimaxFrequencyGridMatchesReference)
{
    const auto grid = ModuleRI::SternheimerRPA::generate_greenx_minimax_frequency_grid(10, 2.0, 30.0);

    const std::vector<double> expected_nodes = {0.37397426,
                                                1.17910655,
                                                2.16745299,
                                                3.50267933,
                                                5.42267468,
                                                8.30517692,
                                                12.81707194,
                                                20.35677335,
                                                34.86933737,
                                                75.79619716};
    const std::vector<double> expected_weights = {0.75729004,
                                                  0.87356340,
                                                  1.12983421,
                                                  1.57927297,
                                                  2.32153469,
                                                  3.54877702,
                                                  5.68714236,
                                                  9.94234857,
                                                  21.24030777,
                                                  80.07961331};

    ASSERT_EQ(grid.omega_ha.size(), expected_nodes.size());
    ASSERT_EQ(grid.weights_ha.size(), expected_weights.size());
    for (std::size_t i = 0; i != expected_nodes.size(); ++i)
    {
        EXPECT_NEAR(grid.omega_ha[i], expected_nodes[i], 1.0e-8);
        EXPECT_NEAR(grid.weights_ha[i], expected_weights[i], 1.0e-8);
    }
}

TEST(SternheimerRPA, ReadFrequencyGridFileAcceptsTwoAndThreeColumnRows)
{
    const std::string filename = ::testing::TempDir() + "/sternheimer_frequency_grid.dat";
    {
        std::ofstream out(filename.c_str());
        out << "# omega_Ha weight_Ha\n";
        out << "0.25 0.75\n";
        out << "2 0.50 1.25 # optional one-based index column\n";
    }

    const auto grid = ModuleRI::SternheimerRPA::read_frequency_grid_file(filename, 2);

    ASSERT_EQ(grid.omega_ha.size(), 2);
    ASSERT_EQ(grid.weights_ha.size(), 2);
    EXPECT_DOUBLE_EQ(grid.omega_ha[0], 0.25);
    EXPECT_DOUBLE_EQ(grid.weights_ha[0], 0.75);
    EXPECT_DOUBLE_EQ(grid.omega_ha[1], 0.50);
    EXPECT_DOUBLE_EQ(grid.weights_ha[1], 1.25);

    std::remove(filename.c_str());
}

TEST(SternheimerRPA, ReadFrequencyGridFileRejectsInvalidRows)
{
    const std::string count_mismatch = ::testing::TempDir() + "/sternheimer_frequency_grid_count.dat";
    {
        std::ofstream out(count_mismatch.c_str());
        out << "0.25 0.75\n";
    }
    EXPECT_THROW(ModuleRI::SternheimerRPA::read_frequency_grid_file(count_mismatch, 2), std::runtime_error);
    std::remove(count_mismatch.c_str());

    const std::string invalid_value = ::testing::TempDir() + "/sternheimer_frequency_grid_invalid.dat";
    {
        std::ofstream out(invalid_value.c_str());
        out << "-0.25 0.75\n";
    }
    EXPECT_THROW(ModuleRI::SternheimerRPA::read_frequency_grid_file(invalid_value, 1), std::invalid_argument);
    std::remove(invalid_value.c_str());
}

TEST(SternheimerRPA, FrequencyOwnerRankSupportsRankShift)
{
    EXPECT_EQ(ModuleRI::SternheimerRPA::frequency_owner_rank(0, 6, 0), 0);
    EXPECT_EQ(ModuleRI::SternheimerRPA::frequency_owner_rank(5, 6, 0), 5);
    EXPECT_EQ(ModuleRI::SternheimerRPA::frequency_owner_rank(0, 6, 1), 1);
    EXPECT_EQ(ModuleRI::SternheimerRPA::frequency_owner_rank(5, 6, 1), 0);
    EXPECT_EQ(ModuleRI::SternheimerRPA::frequency_owner_rank(0, 1, 3), 0);
    EXPECT_EQ(ModuleRI::SternheimerRPA::frequency_owner_rank(1, 4, -1), 0);
    EXPECT_THROW(ModuleRI::SternheimerRPA::frequency_owner_rank(-1, 6, 0), std::invalid_argument);
    EXPECT_THROW(ModuleRI::SternheimerRPA::frequency_owner_rank(0, 0, 0), std::invalid_argument);
}

TEST(SternheimerRPA, AssignsChannelMPIWithinFrequencyGroups)
{
    const auto leader = ModuleRI::SternheimerRPA::frequency_mpi_assignment(0, 16, 32, 2, 1, true);
    EXPECT_TRUE(leader.owns_frequency);
    EXPECT_EQ(leader.frequency_leader_rank, 2);
    EXPECT_EQ(leader.frequency_group_size, 2);
    EXPECT_EQ(leader.frequency_group_local_rank, 0);

    const auto partner = ModuleRI::SternheimerRPA::frequency_mpi_assignment(0, 16, 32, 3, 1, true);
    EXPECT_TRUE(partner.owns_frequency);
    EXPECT_EQ(partner.frequency_leader_rank, 2);
    EXPECT_EQ(partner.frequency_group_size, 2);
    EXPECT_EQ(partner.frequency_group_local_rank, 1);

    const auto other = ModuleRI::SternheimerRPA::frequency_mpi_assignment(0, 16, 32, 4, 1, true);
    EXPECT_FALSE(other.owns_frequency);
    EXPECT_EQ(other.frequency_group_local_rank, -1);

    const auto wrapped = ModuleRI::SternheimerRPA::frequency_mpi_assignment(15, 16, 32, 1, 1, true);
    EXPECT_TRUE(wrapped.owns_frequency);
    EXPECT_EQ(wrapped.frequency_leader_rank, 0);
    EXPECT_EQ(wrapped.frequency_group_local_rank, 1);

    EXPECT_EQ(ModuleRI::SternheimerRPA::channel_group_owner(0, 0, 428, 2), 0);
    EXPECT_EQ(ModuleRI::SternheimerRPA::channel_group_owner(0, 1, 428, 2), 1);
    EXPECT_EQ(ModuleRI::SternheimerRPA::channel_group_owner(1, 0, 428, 2), 0);
    EXPECT_EQ(ModuleRI::SternheimerRPA::channel_group_owner(1, 1, 428, 2), 1);

    EXPECT_THROW(ModuleRI::SternheimerRPA::frequency_mpi_assignment(0, 16, 30, 0, 0, true),
                 std::invalid_argument);
    EXPECT_THROW(ModuleRI::SternheimerRPA::channel_group_owner(0, 0, 428, 0), std::invalid_argument);
}

TEST(SternheimerRPA, AssignsGlobalEquationsUniquelyAndEvenly)
{
    EXPECT_EQ(ModuleRI::SternheimerRPA::global_equation_owner(0, 0, 0, 16, 428, 32, 0), 0);
    EXPECT_EQ(ModuleRI::SternheimerRPA::global_equation_owner(0, 0, 31, 16, 428, 32, 0), 31);
    EXPECT_EQ(ModuleRI::SternheimerRPA::global_equation_owner(0, 0, 32, 16, 428, 32, 0), 0);
    EXPECT_EQ(ModuleRI::SternheimerRPA::global_equation_owner(0, 0, 0, 16, 428, 32, 1), 1);

    constexpr int occupied_count = 2;
    constexpr int frequency_count = 3;
    constexpr int channel_count = 5;
    constexpr int mpi_ranks = 4;
    std::vector<int> owner_counts(mpi_ranks, 0);
    for (int occupied = 0; occupied != occupied_count; ++occupied)
    {
        for (int frequency = 0; frequency != frequency_count; ++frequency)
        {
            for (int channel = 0; channel != channel_count; ++channel)
            {
                const int owner = ModuleRI::SternheimerRPA::global_equation_owner(occupied,
                                                                                  frequency,
                                                                                  channel,
                                                                                  frequency_count,
                                                                                  channel_count,
                                                                                  mpi_ranks,
                                                                                  0);
                ASSERT_GE(owner, 0);
                ASSERT_LT(owner, mpi_ranks);
                ++owner_counts[owner];
            }
        }
    }
    int minimum_count = owner_counts.front();
    int maximum_count = owner_counts.front();
    int total_count = 0;
    for (const int count: owner_counts)
    {
        minimum_count = std::min(minimum_count, count);
        maximum_count = std::max(maximum_count, count);
        total_count += count;
    }
    EXPECT_EQ(total_count, occupied_count * frequency_count * channel_count);
    EXPECT_LE(maximum_count - minimum_count, 1);

    EXPECT_THROW(ModuleRI::SternheimerRPA::global_equation_owner(-1, 0, 0, 3, 5, 4, 0),
                 std::invalid_argument);
    EXPECT_THROW(ModuleRI::SternheimerRPA::global_equation_owner(0, -1, 0, 3, 5, 4, 0),
                 std::invalid_argument);
    EXPECT_THROW(ModuleRI::SternheimerRPA::global_equation_owner(0, 3, 0, 3, 5, 4, 0),
                 std::invalid_argument);
    EXPECT_THROW(ModuleRI::SternheimerRPA::global_equation_owner(0, 0, 5, 3, 5, 4, 0),
                 std::invalid_argument);
    EXPECT_THROW(ModuleRI::SternheimerRPA::global_equation_owner(0, 0, 0, 0, 5, 4, 0),
                 std::invalid_argument);
    EXPECT_THROW(ModuleRI::SternheimerRPA::global_equation_owner(0, 0, 0, 3, 0, 4, 0),
                 std::invalid_argument);
    EXPECT_THROW(ModuleRI::SternheimerRPA::global_equation_owner(0, 0, 0, 3, 5, 0, 0),
                 std::invalid_argument);
}

TEST(SternheimerRPA, ValidatesGlobalEquationMPILayoutRequirements)
{
    EXPECT_NO_THROW(ModuleRI::SternheimerRPA::validate_mpi_layout("frequency_grouped",
                                                                   true,
                                                                   true,
                                                                   true,
                                                                   false,
                                                                   16,
                                                                   32));
    EXPECT_NO_THROW(ModuleRI::SternheimerRPA::validate_mpi_layout("frequency_grouped",
                                                                   true,
                                                                   true,
                                                                   false,
                                                                   true,
                                                                   16,
                                                                   32));
    EXPECT_THROW(ModuleRI::SternheimerRPA::validate_mpi_layout("frequency_grouped",
                                                                true,
                                                                true,
                                                                true,
                                                                false,
                                                                16,
                                                                30),
                 std::invalid_argument);

    EXPECT_NO_THROW(ModuleRI::SternheimerRPA::validate_mpi_layout("global_equation",
                                                                   true,
                                                                   true,
                                                                   true,
                                                                   false,
                                                                   6,
                                                                   4));
    EXPECT_NO_THROW(ModuleRI::SternheimerRPA::validate_mpi_layout("global_equation",
                                                                   true,
                                                                   true,
                                                                   false,
                                                                   true,
                                                                   6,
                                                                   4));
    EXPECT_NO_THROW(ModuleRI::SternheimerRPA::validate_mpi_layout("global_equation",
                                                                   true,
                                                                   true,
                                                                   true,
                                                                   true,
                                                                   6,
                                                                   4));
    EXPECT_THROW(ModuleRI::SternheimerRPA::validate_mpi_layout("global_equation",
                                                                false,
                                                                true,
                                                                true,
                                                                false,
                                                                6,
                                                                4),
                 std::invalid_argument);
    EXPECT_THROW(ModuleRI::SternheimerRPA::validate_mpi_layout("global_equation",
                                                                true,
                                                                false,
                                                                true,
                                                                false,
                                                                6,
                                                                4),
                 std::invalid_argument);
    EXPECT_THROW(ModuleRI::SternheimerRPA::validate_mpi_layout("global_equation",
                                                                true,
                                                                true,
                                                                false,
                                                                false,
                                                                6,
                                                                4),
                 std::invalid_argument);
    EXPECT_THROW(ModuleRI::SternheimerRPA::validate_mpi_layout("invalid",
                                                                true,
                                                                true,
                                                                true,
                                                                false,
                                                                6,
                                                                4),
                 std::invalid_argument);
}

TEST(SternheimerRPA, SolveBiCGStabDiagonalComplexSystem)
{
    const Vector diagonal = {Complex(2.0, 1.0), Complex(3.0, -0.5), Complex(4.0, 2.0)};
    const Vector exact = {Complex(1.0, -0.5), Complex(-2.0, 1.0), Complex(0.25, 0.75)};
    Vector rhs(exact.size());
    for (std::size_t i = 0; i != exact.size(); ++i)
    {
        rhs[i] = diagonal[i] * exact[i];
    }

    ModuleRI::SternheimerRPA::LinearProblem problem;
    problem.apply = [&diagonal](const Vector& input, Vector& output) {
        output.resize(input.size());
        for (std::size_t i = 0; i != input.size(); ++i)
        {
            output[i] = diagonal[i] * input[i];
        }
    };
    problem.dot = dot;

    ModuleRI::SternheimerRPA::SolverOptions options;
    options.max_iter = 50;
    options.residual_tol = 1.0e-12;

    Vector solution(exact.size(), Complex(0.0, 0.0));
    const auto result = ModuleRI::SternheimerRPA::solve_bicgstab(problem, rhs, solution, options);

    EXPECT_TRUE(result.converged);
    EXPECT_LE(result.relative_residual, options.residual_tol);
    for (std::size_t i = 0; i != exact.size(); ++i)
    {
        EXPECT_NEAR(solution[i].real(), exact[i].real(), 1.0e-10);
        EXPECT_NEAR(solution[i].imag(), exact[i].imag(), 1.0e-10);
    }
}

TEST(SternheimerRPA, SolveGMRESDiagonalComplexSystem)
{
    const Vector diagonal = {Complex(2.0, 1.0), Complex(3.0, -0.5), Complex(4.0, 2.0)};
    const Vector exact = {Complex(1.0, -0.5), Complex(-2.0, 1.0), Complex(0.25, 0.75)};
    Vector rhs(exact.size());
    for (std::size_t i = 0; i != exact.size(); ++i)
    {
        rhs[i] = diagonal[i] * exact[i];
    }

    ModuleRI::SternheimerRPA::LinearProblem problem;
    problem.apply = [&diagonal](const Vector& input, Vector& output) {
        output.resize(input.size());
        for (std::size_t i = 0; i != input.size(); ++i)
        {
            output[i] = diagonal[i] * input[i];
        }
    };
    problem.dot = dot;

    ModuleRI::SternheimerRPA::SolverOptions options;
    options.max_iter = 10;
    options.residual_tol = 1.0e-12;

    Vector solution(exact.size(), Complex(0.0, 0.0));
    const auto result = ModuleRI::SternheimerRPA::solve_gmres(problem, rhs, solution, options, 3);

    EXPECT_TRUE(result.converged);
    EXPECT_LE(result.relative_residual, options.residual_tol);
    for (std::size_t i = 0; i != exact.size(); ++i)
    {
        EXPECT_NEAR(solution[i].real(), exact[i].real(), 1.0e-10);
        EXPECT_NEAR(solution[i].imag(), exact[i].imag(), 1.0e-10);
    }
}

TEST(SternheimerRPA, SolveGMRESBatchMatchesIndependentScalarHistories)
{
    using Matrix = std::vector<Vector>;
    const Vector diagonal
        = {Complex(1.5, 0.2), Complex(2.5, -0.3), Complex(4.0, 0.7), Complex(7.0, -0.4)};
    const Matrix exact
        = {{Complex(0.5, -0.25), Complex(-1.0, 0.75), Complex(0.2, 0.4), Complex(-0.3, 0.1)},
           {Complex(1.25, -0.5), Complex(0.0, 0.0), Complex(0.0, 0.0), Complex(0.0, 0.0)},
           {Complex(-0.75, 0.2), Complex(0.6, -0.9), Complex(1.1, 0.3), Complex(-0.4, 0.8)}};
    Matrix rhs = exact;
    for (Vector& column: rhs)
    {
        for (std::size_t i = 0; i != column.size(); ++i)
        {
            column[i] *= diagonal[i];
        }
    }

    std::vector<std::size_t> batch_apply_widths;
    ModuleRI::SternheimerRPA::BatchLinearProblem batch_problem;
    batch_problem.apply = [&diagonal, &batch_apply_widths](const Matrix& input, Matrix& output) {
        batch_apply_widths.push_back(input.size());
        output = input;
        for (Vector& column: output)
        {
            for (std::size_t i = 0; i != column.size(); ++i)
            {
                column[i] *= diagonal[i];
            }
        }
    };
    batch_problem.dot = dot;

    ModuleRI::SternheimerRPA::LinearProblem scalar_problem;
    scalar_problem.apply = [&diagonal](const Vector& input, Vector& output) {
        output = input;
        for (std::size_t i = 0; i != output.size(); ++i)
        {
            output[i] *= diagonal[i];
        }
    };
    scalar_problem.dot = dot;

    ModuleRI::SternheimerRPA::SolverOptions options;
    options.max_iter = 30;
    options.residual_tol = 1.0e-11;
    Matrix scalar_solutions
        = {exact[0], Vector(diagonal.size(), Complex(0.0, 0.0)), Vector(diagonal.size(), Complex(0.0, 0.0))};
    Matrix batch_solutions = scalar_solutions;
    std::vector<ModuleRI::SternheimerRPA::SolverResult> scalar_results;
    for (std::size_t column = 0; column != rhs.size(); ++column)
    {
        scalar_results.push_back(ModuleRI::SternheimerRPA::solve_gmres(
            scalar_problem, rhs[column], scalar_solutions[column], options, 2));
    }

    const auto batch_results
        = ModuleRI::SternheimerRPA::solve_gmres_batch(batch_problem, rhs, batch_solutions, options, 2);

    ASSERT_EQ(batch_results.size(), scalar_results.size());
    ASSERT_EQ(batch_solutions.size(), scalar_solutions.size());
    ASSERT_FALSE(batch_apply_widths.empty());
    EXPECT_EQ(batch_apply_widths.front(), 3U);
    EXPECT_NE(std::find(batch_apply_widths.begin(), batch_apply_widths.end(), 2U), batch_apply_widths.end());
    EXPECT_NE(std::find(batch_apply_widths.begin(), batch_apply_widths.end(), 1U), batch_apply_widths.end());
    for (std::size_t column = 0; column != rhs.size(); ++column)
    {
        EXPECT_EQ(batch_results[column].converged, scalar_results[column].converged);
        EXPECT_EQ(batch_results[column].iterations, scalar_results[column].iterations);
        EXPECT_NEAR(batch_results[column].absolute_residual,
                    scalar_results[column].absolute_residual,
                    1.0e-13);
        EXPECT_NEAR(batch_results[column].relative_residual,
                    scalar_results[column].relative_residual,
                    1.0e-13);
        ASSERT_EQ(batch_solutions[column].size(), scalar_solutions[column].size());
        for (std::size_t i = 0; i != scalar_solutions[column].size(); ++i)
        {
            EXPECT_NEAR(batch_solutions[column][i].real(), scalar_solutions[column][i].real(), 1.0e-12);
            EXPECT_NEAR(batch_solutions[column][i].imag(), scalar_solutions[column][i].imag(), 1.0e-12);
        }
    }
}
