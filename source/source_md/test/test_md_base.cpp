#include "source_md/md_base.h"
#include "source_base/constants.h"
#include "source_base/parallel_cell.h"

#include <gtest/gtest.h>

#include <fstream>
#include <vector>

namespace
{
MDCell make_position_test_cell()
{
    ModuleBase::Matrix3 lattice;
    lattice.e11 = 1.0;
    lattice.e22 = 1.0;
    lattice.e33 = 1.0;
    const ModuleBase::Matrix3 inverse = lattice.Inverse();
    const ModuleBase::CommunicationDomain domain = ModuleBase::world_comm_domain();
    std::vector<LocalAtom> atoms(2);
    for (int i = 0; i < 2; ++i)
    {
        atoms[i].type_index = i;
        atoms[i].mass = 1.0;
        atoms[i].mbl.set(1, 1, 1);
        atoms[i].frac.set(0.25, 0.5, 0.75);
        atoms[i].cart = atoms[i].frac * lattice;
    }
    const std::vector<std::string> labels(1, "X");
    const std::vector<double> masses(1, 1.0);
    const std::vector<std::int64_t> counts(1, 2);
    MDCell cell;
    cell.initialize_from_owned_atoms(lattice, inverse, 1.0, 1.0, 2, atoms, labels, masses, counts, 0.0, domain);
    return cell;
}

TEST(MdPositionTest, TinyNegativeDisplacementRemainsInsideHalfOpenCell)
{
    MDCell cell = make_position_test_cell();
    MD_para parameters;
    parameters.md_dt = ModuleBase::AU_to_FS;
    parameters.md_tfirst = 0.0;
    parameters.md_tlast = 0.0;
    MD_base integrator(parameters, false, true, 0, cell);
    LocalAtom& moving = cell.owned_atoms()[0];
    moving.frac.set(0.0, 0.0, 0.0);
    moving.cart.set(0.0, 0.0, 0.0);
    moving.vel.set(-1.0e-21, -2.0e-21, -3.0e-21);
    std::ofstream log;
    integrator.first_half(log);
    for (int axis = 0; axis < 3; ++axis)
    {
        EXPECT_GE(moving.frac[axis], 0.0);
        EXPECT_LT(moving.frac[axis], 1.0);
        EXPECT_DOUBLE_EQ(moving.frac[axis], 0.0);
        EXPECT_DOUBLE_EQ(moving.cart[axis], 0.0);
        EXPECT_EQ(moving.boundary_shift[axis], 0);
        EXPECT_DOUBLE_EQ(moving.displacement[axis], moving.vel[axis]);
    }
}

TEST(MdPositionTest, FiniteBoundaryCrossingsAndFrozenComponentsKeepTheirPositions)
{
    MDCell cell = make_position_test_cell();
    MD_para parameters;
    parameters.md_dt = ModuleBase::AU_to_FS;
    parameters.md_tfirst = 0.0;
    parameters.md_tlast = 0.0;
    MD_base integrator(parameters, false, true, 0, cell);
    cell.owned_atoms()[0].vel.set(-0.5, 0.75, -1.25);
    cell.owned_atoms()[1].mbl.set(0, 0, 0);
    cell.owned_atoms()[1].vel.set(-0.5, 1.75, -1.25);
    std::ofstream log;
    integrator.first_half(log);
    const LocalAtom& moving = cell.owned_atoms()[0];
    EXPECT_DOUBLE_EQ(moving.frac.x, 0.75);
    EXPECT_DOUBLE_EQ(moving.frac.y, 0.25);
    EXPECT_DOUBLE_EQ(moving.frac.z, 0.5);
    EXPECT_EQ(moving.boundary_shift.x, 1);
    EXPECT_EQ(moving.boundary_shift.y, -1);
    EXPECT_EQ(moving.boundary_shift.z, 1);
    const LocalAtom& frozen = cell.owned_atoms()[1];
    EXPECT_DOUBLE_EQ(frozen.frac.x, 0.25);
    EXPECT_DOUBLE_EQ(frozen.frac.y, 0.5);
    EXPECT_DOUBLE_EQ(frozen.frac.z, 0.75);
}
TEST(MdPositionTest, ExactUpperBoundaryIsWrappedAndTheNextStepClearsItsShift)
{
    MDCell cell = make_position_test_cell();
    MD_para parameters;
    parameters.md_dt = ModuleBase::AU_to_FS;
    parameters.md_tfirst = 0.0;
    parameters.md_tlast = 0.0;
    MD_base integrator(parameters, false, true, 0, cell);
    LocalAtom& moving = cell.owned_atoms()[0];
    moving.vel.set(0.75, 0.0, 0.0);
    std::ofstream log;
    integrator.first_half(log);
    EXPECT_DOUBLE_EQ(moving.frac.x, 0.0);
    EXPECT_DOUBLE_EQ(moving.displacement.x, 0.75);
    EXPECT_EQ(moving.boundary_shift.x, -1);
    moving.vel.set(0.0, 0.0, 0.0);
    integrator.first_half(log);
    EXPECT_DOUBLE_EQ(moving.frac.x, 0.0);
    EXPECT_DOUBLE_EQ(moving.displacement.x, 0.0);
    EXPECT_EQ(moving.boundary_shift.x, 0);
}

TEST(MdPositionTest, RejectsMovementBeyondTheLegacySingleCellRange)
{
    MDCell cell = make_position_test_cell();
    MD_para parameters;
    parameters.md_dt = ModuleBase::AU_to_FS;
    parameters.md_tfirst = 0.0;
    parameters.md_tlast = 0.0;
    MD_base integrator(parameters, false, true, 0, cell);
    LocalAtom& moving = cell.owned_atoms()[0];
    std::ofstream log;
    for (const double velocity : {-1.5, 2.0})
    {
        moving.frac.x = 0.25;
        moving.vel.x = velocity;
        EXPECT_EXIT(integrator.first_half(log), testing::ExitedWithCode(1), "");
    }
}
} // namespace
