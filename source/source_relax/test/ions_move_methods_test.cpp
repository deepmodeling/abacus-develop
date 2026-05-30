#include "for_test.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#define private public
#include "source_relax/ions_move_methods.h"
#undef private
/************************************************
 *  unit tests of class Ions_Move_Methods
 ***********************************************/

/**
 * - Tested Functions:
 *   - Ions_Move_Methods::allocate()
 *   - Ions_Move_Methods::cal_movement()
 *   - Ions_Move_Methods::get_converged()
 *   - Ions_Move_Methods::get_ediff()
 *   - Ions_Move_Methods::get_largest_grad()
 *   - Ions_Move_Methods::get_trust_radius()
 *   - Ions_Move_Methods::get_update_iter()
 */

 // Mock the remake_cell function from update_cell.h
namespace unitcell
{
    // Track if remake_cell was called and with what lattice
    static bool remake_cell_called = false;
    static std::string remake_cell_latName;
    static ModuleBase::Matrix3 remake_cell_latvec;

    void remake_cell(Lattice& lat)
    {
        remake_cell_called = true;
        remake_cell_latName = lat.latName;
        remake_cell_latvec = lat.latvec;

        // Mock implementation: enforce simple cubic structure for "sc"
        if (lat.latName == "sc")
        {
            double celldm = std::sqrt(lat.latvec.e11 * lat.latvec.e11 +
                                     lat.latvec.e12 * lat.latvec.e12 +
                                     lat.latvec.e13 * lat.latvec.e13);
            lat.latvec.Zero();
            lat.latvec.e11 = celldm;
            lat.latvec.e22 = celldm;
            lat.latvec.e33 = celldm;
        }
        // Mock implementation: enforce FCC structure for "fcc"
        else if (lat.latName == "fcc")
        {
            double celldm = std::sqrt(lat.latvec.e11 * lat.latvec.e11 +
                                     lat.latvec.e12 * lat.latvec.e12 +
                                     lat.latvec.e13 * lat.latvec.e13) / std::sqrt(2.0);
            lat.latvec.e11 = -celldm;
            lat.latvec.e12 = 0.0;
            lat.latvec.e13 = celldm;
            lat.latvec.e21 = 0.0;
            lat.latvec.e22 = celldm;
            lat.latvec.e23 = celldm;
            lat.latvec.e31 = -celldm;
            lat.latvec.e32 = celldm;
            lat.latvec.e33 = 0.0;
        }
    }

    // Helper function to reset mock state
    void reset_remake_cell_mock()
    {
        remake_cell_called = false;
        remake_cell_latName = "";
        remake_cell_latvec.Zero();
    }

    // Helper function to check if remake_cell was called
    bool was_remake_cell_called()
    {
        return remake_cell_called;
    }
}

// Define a fixture for the tests
class IonsMoveMethodsTest : public ::testing::Test
{
  protected:
    Ions_Move_Methods imm;
    const int natom = 2;

    virtual void SetUp()
    {
        // Initialize variables before each test
    }

    virtual void TearDown()
    {
        // Clean up after each test
    }
};

// Test the allocate() function
TEST_F(IonsMoveMethodsTest, Allocate)
{
    Ions_Move_Basic::relax_method[0] = "bfgs";
    imm.allocate(natom);
    EXPECT_EQ(Ions_Move_Basic::dim, 6);

    Ions_Move_Basic::relax_method[0] = "sd";
    imm.allocate(natom);
    EXPECT_EQ(Ions_Move_Basic::dim, 6);

    Ions_Move_Basic::relax_method[0] = "cg";
    imm.allocate(natom);
    EXPECT_EQ(Ions_Move_Basic::dim, 6);

    Ions_Move_Basic::relax_method[0] = "cg_bfgs";
    imm.allocate(natom);
    EXPECT_EQ(Ions_Move_Basic::dim, 6);
}

// Test the allocate() function warning quit
TEST_F(IonsMoveMethodsTest, AllocateWarningQuit)
{
    Ions_Move_Basic::relax_method[0] = "none";
    GlobalV::ofs_warning.open("log");
    imm.allocate(natom);
    GlobalV::ofs_warning.close();

    std::ifstream ifs("log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(output, testing::HasSubstr("the parameter Ions_Move_Basic::relax_method is not correct."));
    ifs.close();
    std::remove("log");
}

// Test the cal_movement() function
TEST_F(IonsMoveMethodsTest, CalMovement)
{
    const int istep = 0;
    const int force_step = 1;
    const ModuleBase::matrix f(3, 3);
    const double etot = 0.0;
    UnitCell ucell;

    Ions_Move_Basic::relax_method[0] = "bfgs";
    imm.allocate(natom);
    imm.cal_movement(istep, force_step, f, etot, ucell);
    EXPECT_EQ(Ions_Move_Basic::istep, force_step);

    Ions_Move_Basic::relax_method[0] = "sd";
    imm.allocate(natom);
    imm.cal_movement(istep, force_step, f, etot, ucell);
    EXPECT_EQ(Ions_Move_Basic::istep, force_step);

    Ions_Move_Basic::relax_method[0] = "cg";
    imm.allocate(natom);
    imm.cal_movement(istep, force_step, f, etot, ucell);
    EXPECT_EQ(Ions_Move_Basic::istep, force_step);

    Ions_Move_Basic::relax_method[0] = "cg_bfgs";
    imm.allocate(natom);
    imm.cal_movement(istep, force_step, f, etot, ucell);
    EXPECT_EQ(Ions_Move_Basic::istep, force_step);
}

// Test the cal_movement() function warning quit
TEST_F(IonsMoveMethodsTest, CalMovementWarningQuit)
{
    const int istep = 0;
    const int force_step = 1;
    const ModuleBase::matrix f(3, 3);
    const double etot = 0.0;
    UnitCell ucell;
    Ions_Move_Basic::relax_method[0] = "none";
    imm.allocate(natom);

    GlobalV::ofs_warning.open("log");
    imm.cal_movement(istep, force_step, f, etot, ucell);
    GlobalV::ofs_warning.close();

    std::ifstream ifs("log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(output, testing::HasSubstr("the parameter Ions_Move_Basic::relax_method is not correct."));
    ifs.close();
    std::remove("log");
}

// Test the get_converged() function
TEST_F(IonsMoveMethodsTest, GetConverged)
{
    Ions_Move_Basic::converged = true;

    EXPECT_EQ(imm.get_converged(), true);
}

// Test the get_ediff() function
TEST_F(IonsMoveMethodsTest, GetEdiff)
{
    Ions_Move_Basic::ediff = 1.0;

    EXPECT_DOUBLE_EQ(imm.get_ediff(), 1.0);
}

// Test the get_largest_grad() function
TEST_F(IonsMoveMethodsTest, GetLargestGrad)
{
    Ions_Move_Basic::largest_grad = 2.0;

    EXPECT_DOUBLE_EQ(imm.get_largest_grad(), 2.0);
}

// Test the get_trust_radius() function
TEST_F(IonsMoveMethodsTest, GetTrustRadius)
{
    Ions_Move_Basic::trust_radius = 3.0;

    EXPECT_DOUBLE_EQ(imm.get_trust_radius(), 3.0);
}

// Test the get_update_iter() function
TEST_F(IonsMoveMethodsTest, GetUpdateIter)
{
    Ions_Move_Basic::update_iter = 4;

    EXPECT_EQ(imm.get_update_iter(), 4);
}