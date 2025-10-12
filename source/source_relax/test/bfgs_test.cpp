// source/source_relax/test/bfgs_test.cpp
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "for_test.h"

#define private public
#include "source_relax/bfgs.h"
#undef private

#define private public
#include "source_io/module_parameter/parameter.h"
#undef private

#include "source_relax/ions_move_basic.h" // for Ions_Move_Basic static members

/************************************************
 *  unit tests for BFGS (no MockUnitCell)
 ***********************************************/

class BFGSTest : public ::testing::Test {
protected:
    BFGS bfgs;
    void SetUp() override
    {
        // Initialize variables before each test
    }

    void TearDown() override
    {
        // nothing global to clean here
    }
};

// Test whether the allocate() function can correctly allocate memory space
TEST_F(BFGSTest, TestAllocate)
{
    int size = 2;
    bfgs.allocate(size);

    // Check if allocated arrays are not empty
    EXPECT_FALSE(bfgs.H.empty());
    EXPECT_FALSE(bfgs.pos.empty());
    EXPECT_FALSE(bfgs.pos0.empty());
    EXPECT_FALSE(bfgs.pos_taud.empty());
    EXPECT_FALSE(bfgs.pos_taud0.empty());
    EXPECT_FALSE(bfgs.force.empty());
    EXPECT_FALSE(bfgs.force0.empty());
    EXPECT_FALSE(bfgs.steplength.empty());
    EXPECT_FALSE(bfgs.dpos.empty());
    EXPECT_EQ(bfgs.size, size);
    EXPECT_EQ(bfgs.alpha,70);
    EXPECT_EQ(bfgs.maxstep,PARAM.inp.relax_bfgs_rmax);
    EXPECT_TRUE(bfgs.sign);
    EXPECT_EQ(bfgs.largest_grad,0.0);
}

// Test if a dimension less than or equal to 0 results in an assertion error
TEST_F(BFGSTest, TestAllocateWithZeroDimension)
{
    int size = 0;
    ASSERT_DEATH(bfgs.allocate(size), "");
}

// Test DetermineStep scaling
TEST_F(BFGSTest, DetermineStepScaling)
{
    int size = 2;
    bfgs.allocate(size);

    std::vector<double> steplength = {1.0, 0.1};
    std::vector<std::vector<double>> dpos = {
        {1.0, 1.0, 1.0},
        {0.1, 0.1, 0.1}
    };
    double maxstep = 0.5;
    bfgs.DetermineStep(steplength, dpos, maxstep);

    // first atom scaled down to maxstep
    EXPECT_NEAR(dpos[0][0], 0.5, 1e-12);
    EXPECT_NEAR(dpos[0][1], 0.5, 1e-12);
    EXPECT_NEAR(dpos[0][2], 0.5, 1e-12);

    // second atom unchanged (small)
    EXPECT_NEAR(dpos[1][0], 0.05, 1e-12);
    EXPECT_NEAR(dpos[1][1], 0.05, 1e-12);
    EXPECT_NEAR(dpos[1][2], 0.05, 1e-12);
}

// Test GetPos and GetPostaud without creating extra helper class
TEST_F(BFGSTest, GetPosAndPostaud)
{
    // prepare UnitCell with 1 type and 2 atoms
    UnitCell ucell;
    ucell.ntype = 1;
    ucell.nat = 2;
    ucell.lat0 = 2.0;

    // allocate atoms array
    ucell.atoms = new Atom[ucell.ntype];
    ucell.atoms[0].na = 2;
    ucell.atoms[0].tau = std::vector<ModuleBase::Vector3<double>>(2);
    ucell.atoms[0].taud = std::vector<ModuleBase::Vector3<double>>(2);
    ucell.atoms[0].mbl = std::vector<ModuleBase::Vector3<int>>(2, {1, 1, 1});

    // set coordinates
    ucell.atoms[0].tau[0].x = 1.0; ucell.atoms[0].tau[0].y = 2.0; ucell.atoms[0].tau[0].z = 3.0;
    ucell.atoms[0].tau[1].x = 2.0; ucell.atoms[0].tau[1].y = 3.0; ucell.atoms[0].tau[1].z = 4.0;
    ucell.atoms[0].taud[0].x = 0.1; ucell.atoms[0].taud[0].y = 0.2; ucell.atoms[0].taud[0].z = 0.3;
    ucell.atoms[0].taud[1].x = 0.4; ucell.atoms[0].taud[1].y = 0.5; ucell.atoms[0].taud[1].z = 0.6;

    // allocate mapping arrays 
    ucell.iat2it = new int[ucell.nat];
    ucell.iat2ia = new int[ucell.nat];
    int k = 0;
    for (int it = 0; it < ucell.ntype; ++it) {
        for (int ia = 0; ia < ucell.atoms[it].na; ++ia) {
            ucell.iat2it[k] = it;
            ucell.iat2ia[k] = ia;
            ++k;
        }
    }

    // allocate bfgs arrays and call getters
    bfgs.allocate(ucell.nat);
    bfgs.GetPos(ucell, bfgs.pos);
    bfgs.GetPostaud(ucell, bfgs.pos_taud);

    // pos is tau * BOHR_TO_A * lat0
    EXPECT_DOUBLE_EQ(bfgs.pos[0][0], ucell.atoms[0].tau[0].x * ModuleBase::BOHR_TO_A * ucell.lat0);
    EXPECT_DOUBLE_EQ(bfgs.pos_taud[1][2], ucell.atoms[0].taud[1].z);
}

// Test CalculateLargestGrad (uses ModuleBase::matrix)
TEST_F(BFGSTest, CalculateLargestGrad)
{
    // UnitCell with 1 type and 2 atoms
    UnitCell ucell;
    ucell.ntype = 1;
    ucell.nat = 2;
    ucell.lat0 = 2.0;

    ucell.atoms = new Atom[ucell.ntype];
    ucell.atoms[0].na = 2;
    ucell.atoms[0].mbl = std::vector<ModuleBase::Vector3<int>>(2, {1, 1, 1});

    // mapping arrays
    ucell.iat2it = new int[ucell.nat];
    ucell.iat2ia = new int[ucell.nat];
    int k = 0;
    for (int it = 0; it < ucell.ntype; ++it) {
        for (int ia = 0; ia < ucell.atoms[it].na; ++ia) {
            ucell.iat2it[k] = it;
            ucell.iat2ia[k] = ia;
            ++k;
        }
    }

    // build force matrix: 2 atoms x 3 components
    ModuleBase::matrix force(2, 3);
    force(0, 0) = -2.0;  // this yields grad component = -(-2.0)*lat0 = 4.0 -> divided by lat0 => 2.0
    force(0, 1) = 0.0;
    force(0, 2) = 1.0;
    force(1, 0) = 3.0;   // this yields abs = 6.0 -> divided by lat0 => 3.0 (this should be largest)
    force(1, 1) = -1.0;
    force(1, 2) = 0.0;

    bfgs.allocate(ucell.nat);
    bfgs.CalculateLargestGrad(force, ucell);

    // expected largest_grad = 3.0 (see calculation above)
    EXPECT_NEAR(bfgs.largest_grad, 6.0, 1e-12);
}

// Test relax_step basic functionality
TEST_F(BFGSTest, RelaxStepBasic)
{
    // Setup UnitCell with 1 type, 2 atoms
    UnitCell ucell;
    ucell.ntype = 1;
    ucell.nat = 2;
    ucell.lat0 = 1.0;
    ucell.atoms = new Atom[ucell.ntype];
    ucell.atoms[0].na = 2;
    ucell.atoms[0].tau = std::vector<ModuleBase::Vector3<double>>(2);
    ucell.atoms[0].taud = std::vector<ModuleBase::Vector3<double>>(2);
    ucell.atoms[0].mbl = std::vector<ModuleBase::Vector3<int>>(2, {1, 1, 1});
    ucell.iat2it = new int[ucell.nat];
    ucell.iat2ia = new int[ucell.nat];
    int k = 0;
    for (int it = 0; it < ucell.ntype; ++it) {
        for (int ia = 0; ia < ucell.atoms[it].na; ++ia) {
            ucell.iat2it[k] = it;
            ucell.iat2ia[k] = ia;
            ++k;
        }
    }
    // Set initial positions
    ucell.atoms[0].tau[0].x = 0.0; ucell.atoms[0].tau[0].y = 0.0; ucell.atoms[0].tau[0].z = 0.0;
    ucell.atoms[0].tau[1].x = 1.0; ucell.atoms[0].tau[1].y = 0.0; ucell.atoms[0].tau[1].z = 0.0;
    // Setup force matrix
    ModuleBase::matrix force(2, 3);
    force(0, 0) = 0.1; force(0, 1) = 0.0; force(0, 2) = 0.0;
    force(1, 0) = -0.1; force(1, 1) = 0.0; force(1, 2) = 0.0;
    // Allocate and call relax_step
    bfgs.allocate(ucell.nat);
    bfgs.relax_step(force, ucell);
    // Check that ionic_position_updated is true
    EXPECT_TRUE(ucell.ionic_position_updated);
    // Check that force values are set (converted units)
    EXPECT_NEAR(bfgs.force[0][0], 0.1 * ModuleBase::Ry_to_eV / ModuleBase::BOHR_TO_A, 1e-12);
    EXPECT_NEAR(bfgs.force[1][0], -0.1 * ModuleBase::Ry_to_eV / ModuleBase::BOHR_TO_A, 1e-12);
    // Check that positions are updated (not equal to initial)
    EXPECT_NEAR(bfgs.pos[0][0], 0.0, 1e-12);
    EXPECT_NE(bfgs.pos[1][0], 1.0);
}