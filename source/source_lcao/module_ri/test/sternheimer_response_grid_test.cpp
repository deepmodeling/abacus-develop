#include "source_lcao/module_ri/sternheimer_response_grid.h"

#include "source_base/constants.h"
#include "source_base/matrix3.h"
#include "source_basis/module_pw/pw_basis.h"

#include <cmath>
#include <gtest/gtest.h>
#include <memory>
#include <stdexcept>
#include <vector>

#ifdef __MPI
#include <mpi.h>
#endif

namespace
{

std::unique_ptr<ModulePW::PW_Basis> make_basis(const int nx, const int ny, const int nz)
{
    auto basis = std::make_unique<ModulePW::PW_Basis>("cpu", "double");
#ifdef __MPI
    basis->initmpi(1, 0, MPI_COMM_SELF);
#endif
    basis->initgrids(1.0, ModuleBase::Matrix3(8.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 8.0), nx, ny, nz);
    basis->initparameters(false, 1.0e12);
    basis->fft_bundle.initfftmode(0);
    basis->setuptransform();
    basis->collect_local_pw();
    return basis;
}

std::unique_ptr<ModulePW::PW_Basis> make_cutoff_basis(const double ecutwfc)
{
    auto basis = std::make_unique<ModulePW::PW_Basis>("cpu", "double");
#ifdef __MPI
    basis->initmpi(1, 0, MPI_COMM_SELF);
#endif
    basis->initgrids(1.0,
                     ModuleBase::Matrix3(8.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 8.0),
                     4.0 * ecutwfc);
    basis->initparameters(false, 4.0 * ecutwfc);
    basis->fft_bundle.initfftmode(0);
    basis->setuptransform();
    basis->collect_local_pw();
    return basis;
}

std::size_t index(const int ix, const int iy, const int iz, const int ny, const int nz)
{
    return static_cast<std::size_t>((ix * ny + iy) * nz + iz);
}

std::vector<double> cosine_mode(const int nx, const int ny, const int nz, const int mx, const int my, const int mz)
{
    std::vector<double> values(static_cast<std::size_t>(nx) * ny * nz, 0.0);
    for (int ix = 0; ix != nx; ++ix)
    {
        for (int iy = 0; iy != ny; ++iy)
        {
            for (int iz = 0; iz != nz; ++iz)
            {
                const double phase = ModuleBase::TWO_PI
                                     * (static_cast<double>(mx * ix) / nx + static_cast<double>(my * iy) / ny
                                        + static_cast<double>(mz * iz) / nz);
                values[index(ix, iy, iz, ny, nz)] = std::cos(phase);
            }
        }
    }
    return values;
}

double l2_norm(const std::vector<double>& values)
{
    double norm2 = 0.0;
    for (const double value: values)
    {
        norm2 += value * value;
    }
    return std::sqrt(norm2);
}

double relative_l2_error(const std::vector<double>& actual, const std::vector<double>& expected)
{
    EXPECT_EQ(actual.size(), expected.size());
    std::vector<double> difference(actual.size(), 0.0);
    for (std::size_t index = 0; index != actual.size(); ++index)
    {
        difference[index] = actual[index] - expected[index];
    }
    return l2_norm(difference) / l2_norm(expected);
}

} // namespace

TEST(SternheimerResponseGrid, CutoffDecisionKeepsCompatibilityAndRejectsUpsampling)
{
    EXPECT_FALSE(ModuleRI::sternheimer_uses_independent_response_grid(0.0, 80.0));
    EXPECT_TRUE(ModuleRI::sternheimer_uses_independent_response_grid(30.0, 80.0));
    EXPECT_THROW(ModuleRI::sternheimer_uses_independent_response_grid(-1.0, 80.0), std::invalid_argument);
    EXPECT_THROW(ModuleRI::sternheimer_uses_independent_response_grid(90.0, 80.0), std::invalid_argument);
}

TEST(SternheimerResponseGrid, SameGridRestrictionReturnsInputExactly)
{
    const auto basis = make_basis(6, 6, 6);
    std::vector<double> input(static_cast<std::size_t>(basis->nxyz), 0.0);
    for (std::size_t ir = 0; ir != input.size(); ++ir)
    {
        input[ir] = 0.25 + static_cast<double>(ir) / input.size();
    }

    const auto output = ModuleRI::restrict_sternheimer_real_field(*basis, *basis, input);

    EXPECT_EQ(output, input);
}

TEST(SternheimerResponseGrid, EqualCutoffKeepsTheOriginalPBEGridPath)
{
    const auto pbe_basis = make_cutoff_basis(80.0);

    const auto response_grid = ModuleRI::make_sternheimer_response_grid(*pbe_basis, 80.0, 80.0, 0);

    EXPECT_FALSE(response_grid.independent);
    EXPECT_EQ(response_grid.basis, pbe_basis.get());
    EXPECT_EQ(response_grid.serial_fine_basis, nullptr);
    EXPECT_EQ(response_grid.serial_response_basis, nullptr);
}

TEST(SternheimerResponseGrid, ReducedCutoffBuildsSerialFineAndResponseBases)
{
    const auto pbe_basis = make_cutoff_basis(80.0);

    const auto response_grid = ModuleRI::make_sternheimer_response_grid(*pbe_basis, 80.0, 30.0, 0);

    ASSERT_TRUE(response_grid.independent);
    ASSERT_NE(response_grid.serial_fine_basis, nullptr);
    ASSERT_NE(response_grid.serial_response_basis, nullptr);
    EXPECT_EQ(response_grid.serial_fine_basis->nx, pbe_basis->nx);
    EXPECT_EQ(response_grid.serial_fine_basis->ny, pbe_basis->ny);
    EXPECT_EQ(response_grid.serial_fine_basis->nz, pbe_basis->nz);
    EXPECT_LT(response_grid.serial_response_basis->nxyz, pbe_basis->nxyz);
    EXPECT_EQ(response_grid.basis, response_grid.serial_response_basis.get());
}

TEST(SternheimerResponseGrid, RestrictionPreservesConstant)
{
    const auto fine = make_basis(8, 8, 8);
    const auto coarse = make_basis(4, 4, 4);
    const std::vector<double> input(static_cast<std::size_t>(fine->nxyz), 3.25);

    const auto output = ModuleRI::restrict_sternheimer_real_field(*fine, *coarse, input);

    ASSERT_EQ(output.size(), static_cast<std::size_t>(coarse->nxyz));
    for (const double value: output)
    {
        EXPECT_NEAR(value, 3.25, 1.0e-12);
    }
}

TEST(SternheimerResponseGrid, RestrictionPreservesRetainedLowGMode)
{
    const auto fine = make_basis(8, 8, 8);
    const auto coarse = make_basis(4, 4, 4);
    const auto input = cosine_mode(8, 8, 8, 1, 1, 0);
    const auto expected = cosine_mode(4, 4, 4, 1, 1, 0);

    const auto output = ModuleRI::restrict_sternheimer_real_field(*fine, *coarse, input);

    EXPECT_LT(relative_l2_error(output, expected), 1.0e-11);
}

TEST(SternheimerResponseGrid, RestrictionRemovesFineOnlyHighGMode)
{
    const auto fine = make_basis(8, 8, 8);
    const auto coarse = make_basis(4, 4, 4);
    const auto input = cosine_mode(8, 8, 8, 3, 0, 0);

    const auto output = ModuleRI::restrict_sternheimer_real_field(*fine, *coarse, input);

    EXPECT_LT(l2_norm(output), 1.0e-11 * l2_norm(input));
}

TEST(SternheimerResponseGrid, RestrictionDoesNotAliasRemovedModeIntoRetainedMode)
{
    const auto fine = make_basis(8, 8, 8);
    const auto coarse = make_basis(4, 4, 4);
    auto input = cosine_mode(8, 8, 8, 1, 0, 0);
    const auto high_mode = cosine_mode(8, 8, 8, 3, 0, 0);
    for (std::size_t ir = 0; ir != input.size(); ++ir)
    {
        input[ir] += 0.3 * high_mode[ir];
    }
    const auto expected = cosine_mode(4, 4, 4, 1, 0, 0);

    const auto output = ModuleRI::restrict_sternheimer_real_field(*fine, *coarse, input);

    EXPECT_LT(relative_l2_error(output, expected), 1.0e-11);
}

int main(int argc, char** argv)
{
#ifdef __MPI
    MPI_Init(&argc, &argv);
#endif
    testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
#ifdef __MPI
    MPI_Finalize();
#endif
    return result;
}
