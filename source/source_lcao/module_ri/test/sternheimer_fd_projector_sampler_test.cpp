#include "source_lcao/module_ri/sternheimer_fd_projector_sampler.h"

#include "source_base/matrix.h"
#include "source_base/vector3.h"

#include <cmath>
#include <complex>
#include <gtest/gtest.h>
#include <limits>
#include <stdexcept>
#include <vector>

namespace
{

constexpr double y00()
{
    return 0.28209479177387814347;
}

ModuleRI::SternheimerFDHamiltonian::Grid line_grid(const bool periodic)
{
    ModuleRI::SternheimerFDHamiltonian::Grid grid;
    grid.nx = 3;
    grid.ny = 1;
    grid.nz = 1;
    grid.hx = 1.0;
    grid.hy = 10.0;
    grid.hz = 10.0;
    grid.periodic = periodic;
    return grid;
}

} // namespace

TEST(SternheimerFDProjectorSampler, SamplesSProjectorWithoutMinimumImage)
{
    ModuleRI::SternheimerFDRadialProjectorSet radial_set;
    radial_set.radial_grid = {0.0, 1.0, 2.0};
    radial_set.beta_radials = {{2.0, 4.0, 6.0}};
    radial_set.angular_momenta = {0};
    radial_set.d_radial = {{{1.5, 0.0}}};

    const auto block = ModuleRI::sample_sternheimer_fd_projector_block(
        radial_set, line_grid(false), ModuleBase::Vector3<double>(0.0, 0.0, 0.0));

    ASSERT_EQ(block.projectors.size(), 1);
    ASSERT_EQ(block.projectors[0].size(), 3);
    EXPECT_NEAR(block.projectors[0][0].real(), 2.0 * y00(), 1.0e-14);
    EXPECT_NEAR(block.projectors[0][1].real(), 4.0 * y00(), 1.0e-14);
    EXPECT_NEAR(block.projectors[0][2].real(), 6.0 * y00(), 1.0e-14);
    EXPECT_NEAR(block.d_matrix[0][0].real(), 1.5, 1.0e-14);
}

TEST(SternheimerFDProjectorSampler, WritesProjectorsInABACUSRealSpaceIndexOrder)
{
    ModuleRI::SternheimerFDRadialProjectorSet radial_set;
    radial_set.radial_grid = {0.0, 1.0, 2.0};
    radial_set.beta_radials = {{0.0, 1.0, 2.0}};
    radial_set.angular_momenta = {0};
    radial_set.d_radial = {{{1.0, 0.0}}};

    ModuleRI::SternheimerFDHamiltonian::Grid grid;
    grid.nx = 2;
    grid.ny = 2;
    grid.nz = 2;
    grid.hx = 1.0;
    grid.hy = 1.0;
    grid.hz = 1.0;
    grid.periodic = false;

    const auto block = ModuleRI::sample_sternheimer_fd_projector_block(
        radial_set, grid, ModuleBase::Vector3<double>(0.0, 0.0, 0.0));

    ASSERT_EQ(block.projectors.size(), 1);
    ASSERT_EQ(block.projectors[0].size(), 8);
    EXPECT_NEAR(block.projectors[0][4].real(), 1.0 * y00(), 1.0e-14); // (ix, iy, iz) = (1, 0, 0)
    EXPECT_NEAR(block.projectors[0][1].real(), 1.0 * y00(), 1.0e-14); // (ix, iy, iz) = (0, 0, 1)
    EXPECT_NEAR(block.projectors[0][6].real(), std::sqrt(2.0) * y00(), 1.0e-14); // (1, 1, 0)
}

TEST(SternheimerFDProjectorSampler, BuildsRadialSetFromABACUSMatrices)
{
    ModuleBase::matrix betar(2, 3);
    betar(0, 0) = 0.0;
    betar(0, 1) = 2.0;
    betar(0, 2) = 6.0;
    betar(1, 0) = 0.0;
    betar(1, 1) = 4.0;
    betar(1, 2) = 6.0;

    ModuleBase::matrix dion(2, 2);
    dion(0, 0) = 7.0;
    dion(0, 1) = 8.0;
    dion(1, 0) = 9.0;
    dion(1, 1) = 10.0;

    const auto radial_set = ModuleRI::make_sternheimer_fd_radial_projector_set_from_abacus_matrices(
        {0.0, 1.0, 2.0}, betar, {0, 1}, dion);

    EXPECT_EQ(radial_set.radial_grid.size(), 3);
    ASSERT_EQ(radial_set.beta_radials.size(), 2);
    EXPECT_DOUBLE_EQ(radial_set.beta_radials[0][0], 2.0);
    EXPECT_DOUBLE_EQ(radial_set.beta_radials[0][2], 3.0);
    EXPECT_DOUBLE_EQ(radial_set.beta_radials[1][0], 0.0);
    EXPECT_DOUBLE_EQ(radial_set.beta_radials[1][2], 3.0);
    ASSERT_EQ(radial_set.d_radial.size(), 2);
    EXPECT_DOUBLE_EQ(radial_set.d_radial[0][1].real(), 8.0);
    EXPECT_EQ(radial_set.angular_momenta[1], 1);
}

TEST(SternheimerFDProjectorSampler, PeriodizesAllProjectorImagesAtGamma)
{
    ModuleRI::SternheimerFDRadialProjectorSet radial_set;
    radial_set.radial_grid = {0.0, 1.0, 2.0};
    radial_set.beta_radials = {{0.0, 1.0, 2.0}};
    radial_set.angular_momenta = {0};
    radial_set.d_radial = {{{1.0, 0.0}}};

    const auto block = ModuleRI::sample_sternheimer_fd_projector_block(
        radial_set, line_grid(true), ModuleBase::Vector3<double>(0.0, 0.0, 0.0));

    ASSERT_EQ(block.projectors.size(), 1);
    EXPECT_NEAR(block.projectors[0][0].real(), 0.0, 1.0e-14);
    EXPECT_NEAR(block.projectors[0][1].real(), 3.0 * y00(), 1.0e-14);
    EXPECT_NEAR(block.projectors[0][2].real(), 3.0 * y00(), 1.0e-14);
}

TEST(SternheimerFDProjectorSampler, AppliesBlochPhaseToPeriodicImages)
{
    ModuleRI::SternheimerFDRadialProjectorSet radial_set;
    radial_set.radial_grid = {0.0, 0.5, 1.0};
    radial_set.beta_radials = {{0.0, 1.0, 2.0}};
    radial_set.angular_momenta = {0};
    radial_set.d_radial = {{{1.0, 0.0}}};

    auto grid = line_grid(true);
    grid.kpoint = {0.25, 0.0, 0.0};
    const auto block = ModuleRI::sample_sternheimer_fd_projector_block(
        radial_set, grid, ModuleBase::Vector3<double>(0.0, 0.0, 0.0));

    ASSERT_EQ(block.projectors.size(), 1);
    EXPECT_NEAR(block.projectors[0][2].real(), 0.0, 1.0e-14);
    EXPECT_NEAR(block.projectors[0][2].imag(), 2.0 * y00(), 1.0e-14);
}

TEST(SternheimerFDProjectorSampler, RejectsInvalidBlochKPoint)
{
    ModuleRI::SternheimerFDRadialProjectorSet radial_set;
    radial_set.radial_grid = {0.0, 1.0};
    radial_set.beta_radials = {{0.0, 1.0}};
    radial_set.angular_momenta = {0};
    radial_set.d_radial = {{{1.0, 0.0}}};

    auto grid = line_grid(false);
    grid.kpoint = {0.25, 0.0, 0.0};
    EXPECT_THROW(ModuleRI::sample_sternheimer_fd_projector_block(
                     radial_set, grid, ModuleBase::Vector3<double>(0.0, 0.0, 0.0)),
                 std::invalid_argument);

    grid.periodic = true;
    grid.kpoint = {std::numeric_limits<double>::quiet_NaN(), 0.0, 0.0};
    EXPECT_THROW(ModuleRI::sample_sternheimer_fd_projector_block(
                     radial_set, grid, ModuleBase::Vector3<double>(0.0, 0.0, 0.0)),
                 std::invalid_argument);
}

TEST(SternheimerFDProjectorSampler, ComplexBlochProjectorsGiveHermitianNonlocalOperator)
{
    using Complex = std::complex<double>;
    using Vector = ModuleRI::SternheimerFDNonlocalProjector::Vector;

    ModuleRI::SternheimerFDRadialProjectorSet radial_set;
    radial_set.radial_grid = {0.0, 0.5, 1.0};
    radial_set.beta_radials = {{0.0, 1.0, 2.0}, {0.0, 2.0, 4.0}};
    radial_set.angular_momenta = {0, 0};
    radial_set.d_radial = {{Complex(2.0, 0.0), Complex(0.5, 0.25)},
                           {Complex(0.5, -0.25), Complex(3.0, 0.0)}};

    auto grid = line_grid(true);
    grid.kpoint = {0.25, 0.0, 0.0};
    const auto block = ModuleRI::sample_sternheimer_fd_projector_block(
        radial_set, grid, ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
    const int size = grid.nx * grid.ny * grid.nz;
    ModuleRI::SternheimerFDNonlocalProjector projector(size, 1.0, {block});

    ASSERT_EQ(block.d_matrix.size(), 2);
    EXPECT_EQ(block.d_matrix[0][1], std::conj(block.d_matrix[1][0]));
    EXPECT_GT(std::abs(block.projectors[0][2].imag()), 1.0e-12);

    std::vector<Vector> matrix(size, Vector(size, Complex(0.0, 0.0)));
    for (int column = 0; column != size; ++column)
    {
        Vector basis(size, Complex(0.0, 0.0));
        basis[column] = Complex(1.0, 0.0);
        Vector output;
        projector.apply(basis, output);
        for (int row = 0; row != size; ++row)
        {
            matrix[row][column] = output[row];
        }
    }
    for (int row = 0; row != size; ++row)
    {
        for (int column = 0; column != size; ++column)
        {
            const Complex difference = matrix[row][column] - std::conj(matrix[column][row]);
            EXPECT_NEAR(difference.real(), 0.0, 1.0e-12);
            EXPECT_NEAR(difference.imag(), 0.0, 1.0e-12);
        }
    }
}

TEST(SternheimerFDProjectorSampler, ExpandsDMatrixByAngularChannel)
{
    ModuleRI::SternheimerFDRadialProjectorSet radial_set;
    radial_set.radial_grid = {0.0, 1.0};
    radial_set.beta_radials = {{0.0, 1.0}};
    radial_set.angular_momenta = {1};
    radial_set.d_radial = {{{5.0, 0.0}}};

    const auto block = ModuleRI::sample_sternheimer_fd_projector_block(
        radial_set, line_grid(false), ModuleBase::Vector3<double>(0.0, 0.0, 0.0));

    ASSERT_EQ(block.projectors.size(), 3);
    ASSERT_EQ(block.d_matrix.size(), 3);
    for (int i = 0; i != 3; ++i)
    {
        ASSERT_EQ(block.d_matrix[i].size(), 3);
        for (int j = 0; j != 3; ++j)
        {
            EXPECT_NEAR(block.d_matrix[i][j].real(), i == j ? 5.0 : 0.0, 1.0e-14);
        }
    }
}
