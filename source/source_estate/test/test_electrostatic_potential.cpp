#include "source_estate/module_pot/electrostatic_potential.h"
#include "source_estate/module_pot/efield.h"
#include "source_estate/module_pot/h_hartree_pw.h"
#include "source_base/parallel_global.h"

#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <stdexcept>

// Only the public cell/grid storage is needed by the averaging function.
UnitCell::UnitCell() {}
UnitCell::~UnitCell() {}
Magnetism::Magnetism() {}
Magnetism::~Magnetism() {}
SepPot::SepPot() {}
SepPot::~SepPot() {}
Sep_Cell::Sep_Cell() noexcept {}
Sep_Cell::~Sep_Cell() noexcept {}

// The same production translation unit also constructs potentials. These
// dependencies must remain unused: the tests supply a hand-computed grid.
ModuleBase::matrix elecstate::H_Hartree_pw::v_hartree(
    const UnitCell&, ModulePW::PW_Basis*, const int&, const double* const* const)
{
    throw std::logic_error("unexpected Hartree construction in averaging test");
}
ModuleBase::matrix elecstate::Efield::add_efield(
    const UnitCell&, const ModulePW::PW_Basis*, const int&, const double* const* const, const surchem&)
{
    throw std::logic_error("unexpected dipole construction in averaging test");
}

namespace
{
int process_count = 1;
int process_rank = 0;

class ElectrostaticPotentialTest : public testing::Test
{
  protected:
    void set_grid(int planes)
    {
        cell.lat0_angstrom = 4.0;
        const double z_length = planes / 4.0;
        cell.latvec = ModuleBase::Matrix3(1.0, 0.0, 0.0,
                                         0.0, 1.0, 0.0,
                                         0.0, 0.0, z_length);
        basis.nx = 4;
        basis.ny = 4;
        basis.nz = planes;
        const int remainder = planes % process_count;
        const int extra = process_rank < remainder ? 1 : 0;
        basis.nplane = planes / process_count + extra;
        basis.startz_current = process_rank * (planes / process_count) + std::min(process_rank, remainder);
        basis.nrxx = basis.nx * basis.ny * basis.nplane;
#ifdef __MPI
        basis.pool_world = MPI_COMM_WORLD;
#endif
        potential.clear();
        for (int x = 0; x < basis.nx; ++x)
        {
            for (int y = 0; y < basis.ny; ++y)
            {
                for (int z = basis.startz_current; z < basis.startz_current + basis.nplane; ++z)
                {
                    const double value = 100.0 * x + 10.0 * y + z;
                    potential.push_back(value);
                }
            }
        }
    }

    void expect_average(int axis, double lower, double upper, double expected)
    {
        const double actual = elecstate::planar_potential_average(cell, basis, potential, axis, lower, upper);
        EXPECT_NEAR(actual, expected, 1.0e-12);
    }

    UnitCell cell;
    ModulePW::PW_Basis basis;
    std::vector<double> potential;
};

TEST_F(ElectrostaticPotentialTest, HalfOpenWindowsOnAllAxes)
{
    set_grid(5);
    // V = 100*x + 10*y + z; unrestricted means are x=y=1.5, z=2.
    expect_average(0, 1.0, 2.0, 117.0);
    expect_average(1, 1.0, 2.0, 162.0);
    expect_average(2, 1.0, 3.0, 166.5);
}

TEST_F(ElectrostaticPotentialTest, UnequalPartitionsUsePointCounts)
{
    set_grid(5);
    // Four ranks hold 2/1/1/1 planes; averaging rank means would be wrong.
    expect_average(2, 0.0, 5.0, 167.0);
}

TEST_F(ElectrostaticPotentialTest, EmptyGlobalWindowFailsOnEveryRank)
{
    set_grid(5);
    EXPECT_THROW(elecstate::planar_potential_average(cell, basis, potential, 2, 0.1, 0.2),
                 std::runtime_error);
}

TEST_F(ElectrostaticPotentialTest, EmptyLocalPlanesContributeZero)
{
    set_grid(2);
    // With four ranks, two ranks own no planes.
    if (process_rank >= 2)
    {
        EXPECT_EQ(basis.nplane, 0);
        EXPECT_TRUE(potential.empty());
    }
    expect_average(2, 0.0, 2.0, 165.5);
}

TEST_F(ElectrostaticPotentialTest, RanksOutsideWindowStillParticipate)
{
    set_grid(5);
    expect_average(2, 4.0, 5.0, 169.0);
}

TEST_F(ElectrostaticPotentialTest, TriclinicWindowsUseNormalHeights)
{
    set_grid(5);
    cell.latvec = ModuleBase::Matrix3(1.0, 0.0, 0.0,
                                     0.75, 1.0, 0.0,
                                     0.75, 0.75, 1.0);
    // For this cell the normal heights are 64/sqrt(409), 3.2, and 4 A.
    const double heights[3] = {64.0 / std::sqrt(409.0), 3.2, 4.0};
    const double centers[3] = {0.25, 0.25, 0.2};
    const double expected[3] = {117.0, 162.0, 166.0};
    for (int axis = 0; axis < 3; ++axis)
    {
        const double lower = (centers[axis] - 0.01) * heights[axis];
        const double upper = (centers[axis] + 0.01) * heights[axis];
        expect_average(axis, lower, upper, expected[axis]);
    }
}
} // namespace

int main(int argc, char** argv)
{
    int thread_count = 1;
    Parallel_Global::read_pal_param(argc, argv, process_count, thread_count, process_rank);
#ifdef __MPI
    POOL_WORLD = MPI_COMM_NULL;
    KP_WORLD = MPI_COMM_NULL;
    INT_BGROUP = MPI_COMM_NULL;
    BP_WORLD = MPI_COMM_NULL;
    GRID_WORLD = MPI_COMM_NULL;
    DIAG_WORLD = MPI_COMM_NULL;
#endif
    testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
#ifdef __MPI
    Parallel_Global::finalize_mpi();
#endif
    return result;
}
