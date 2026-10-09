#include "source_relax/relax_criteria.h"
#include "for_test.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>
#include "source_relax/ions_move_methods.h"
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

    void update_pos_tau(const Lattice&, const double*, const int, const int, Atom*)
    {
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
  public:
    Relax_Criteria criteria;

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
    imm.allocate(natom, "bfgs", "1");
    EXPECT_EQ(Ions_Move_Basic::dim, 6);

    imm.allocate(natom, "sd", "1");
    EXPECT_EQ(Ions_Move_Basic::dim, 6);

    imm.allocate(natom, "cg", "1");
    EXPECT_EQ(Ions_Move_Basic::dim, 6);

    imm.allocate(natom, "cg_bfgs", "1");
    EXPECT_EQ(Ions_Move_Basic::dim, 6);
}

// Test the allocate() function warning quit
TEST_F(IonsMoveMethodsTest, AllocateWarningQuit)
{
    GlobalV::ofs_warning.open("log");
    imm.allocate(natom, "none", "1");
    GlobalV::ofs_warning.close();

    std::ifstream ifs("log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(output, testing::HasSubstr("the parameter relax_method is not correct."));
    ifs.close();
    std::remove("log");
}

// Test the cal_movement() function
TEST_F(IonsMoveMethodsTest, CalMovement)
{
    const int istep = 0;
    const int force_step = 1;
    ModuleBase::matrix f(natom, 3);
    f(0, 0) = 0.1;
    f(1, 1) = -0.1;
    const double etot = 0.0;
    UnitCell ucell;
    std::ofstream ofs;

    std::vector<std::string> relax_method;

    relax_method = {"bfgs", "1"};
    imm.allocate(natom, "bfgs", "1");
    imm.cal_movement(istep, force_step, f, etot, ucell, ofs, relax_method, criteria);

    relax_method = {"sd", "1"};
    imm.allocate(natom, "sd", "1");
    imm.cal_movement(istep, force_step, f, etot, ucell, ofs, relax_method, criteria);

    relax_method = {"cg", "1"};
    imm.allocate(natom, "cg", "1");
    imm.cal_movement(istep, force_step, f, etot, ucell, ofs, relax_method, criteria);

    relax_method = {"cg_bfgs", "1"};
    imm.allocate(natom, "cg_bfgs", "1");
    imm.cal_movement(istep, force_step, f, etot, ucell, ofs, relax_method, criteria);
}

// Test the cal_movement() function warning quit
TEST_F(IonsMoveMethodsTest, CalMovementWarningQuit)
{
    const int istep = 0;
    const int force_step = 1;
    const ModuleBase::matrix f(3, 3);
    const double etot = 0.0;
    UnitCell ucell;
    std::ofstream ofs;
    std::vector<std::string> relax_method = {"none", "1"};
    imm.allocate(natom, "none", "1");

    GlobalV::ofs_warning.open("log");
    imm.cal_movement(istep, force_step, f, etot, ucell, ofs, relax_method, criteria);
    GlobalV::ofs_warning.close();

    std::ifstream ifs("log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(output, testing::HasSubstr("the parameter relax_method is not correct."));
    ifs.close();
    std::remove("log");
}



// Test the get_ediff() function
TEST_F(IonsMoveMethodsTest, GetEdiff)
{
    // etot_info_ is initialized to {0.0, 0.0, 0.0} by default
    EXPECT_DOUBLE_EQ(imm.get_ediff(), 0.0);
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
    imm.set_update_iter(4);

    EXPECT_EQ(imm.get_update_iter(), 4);
}

TEST_F(IonsMoveMethodsTest, ResetAfterCellChange)
{
    const std::string log_file = "reset_after_cell_change.log";
    std::ofstream ofs(log_file);

    imm.allocate(natom, "bfgs", "2");
    imm.set_converged(true);
    imm.set_update_iter(4);
    imm.get_etot_info() = {-1.0, -2.0};
    Ions_Move_BFGS& bfgs = imm.get_bfgs();
    bfgs.get_first_step() = false;
    bfgs.get_save_flag() = true;
    bfgs.get_tr_min_hit() = true;
    std::fill(bfgs.get_pos().begin(), bfgs.get_pos().end(), 1.0);
    std::fill(bfgs.get_pos_p().begin(), bfgs.get_pos_p().end(), 2.0);
    std::fill(bfgs.get_grad().begin(), bfgs.get_grad().end(), 3.0);
    std::fill(bfgs.get_grad_p().begin(), bfgs.get_grad_p().end(), 4.0);
    std::fill(bfgs.get_move().begin(), bfgs.get_move().end(), 5.0);
    std::fill(bfgs.get_move_p().begin(), bfgs.get_move_p().end(), 6.0);
    Ions_Move_Basic::trust_radius = 0.3;
    Ions_Move_Basic::trust_radius_old = 0.2;

    imm.reset_after_cell_change({"bfgs", "2"}, ofs);

    EXPECT_FALSE(imm.get_converged());
    EXPECT_EQ(imm.get_update_iter(), 0);
    EXPECT_THAT(imm.get_etot_info(), testing::Each(0.0));
    EXPECT_TRUE(bfgs.get_first_step());
    EXPECT_FALSE(bfgs.get_save_flag());
    EXPECT_FALSE(bfgs.get_tr_min_hit());
    EXPECT_THAT(bfgs.get_pos(), testing::Each(0.0));
    EXPECT_THAT(bfgs.get_pos_p(), testing::Each(0.0));
    EXPECT_THAT(bfgs.get_grad(), testing::Each(0.0));
    EXPECT_THAT(bfgs.get_grad_p(), testing::Each(0.0));
    EXPECT_THAT(bfgs.get_move(), testing::Each(0.0));
    EXPECT_THAT(bfgs.get_move_p(), testing::Each(0.0));
    for (int i = 0; i < Ions_Move_Basic::dim; ++i)
    {
        for (int j = 0; j < Ions_Move_Basic::dim; ++j)
        {
            EXPECT_DOUBLE_EQ(bfgs.get_inv_hess()(i, j), i == j ? 1.0 : 0.0);
        }
    }
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::trust_radius, 0.0);
    EXPECT_DOUBLE_EQ(Ions_Move_Basic::trust_radius_old, 0.0);

    imm.allocate(natom, "bfgs", "1");
    imm.set_converged(true);
    imm.set_update_iter(3);
    imm.get_etot_info() = {-3.0, -4.0};
    ASSERT_TRUE(imm.get_bfgs_trad().get_is_initialized());

    imm.reset_after_cell_change({"bfgs", "1"}, ofs);
    ofs.close();

    EXPECT_FALSE(imm.get_converged());
    EXPECT_EQ(imm.get_update_iter(), 0);
    EXPECT_THAT(imm.get_etot_info(), testing::Each(0.0));
    EXPECT_FALSE(imm.get_bfgs_trad().get_is_initialized());

    std::ifstream ifs(log_file);
    const std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(output, testing::HasSubstr("Reset ionic BFGS history after cell change."));
    EXPECT_THAT(output, testing::HasSubstr("Reset traditional ionic BFGS history after cell change."));
    ifs.close();
    std::remove(log_file.c_str());
}
