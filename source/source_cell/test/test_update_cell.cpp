#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include "source_base/global_variable.h"
#include "source_base/mathzone.h"
#include "source_cell/unitcell.h"
#include "source_cell/update_cell.h"
#include "prepare_unitcell.h"

// The test-only cell_info object library does not contain magnetism.cpp,
// so the Magnetism constructor/destructor must be provided locally.
Magnetism::Magnetism()
{
    this->tot_mag = 0.0;
    this->abs_mag = 0.0;
}
Magnetism::~Magnetism()
{
}

/************************************************
 *  unit test of update_cell.cpp
 ***********************************************/

/**
 * - Tested Functions:
 *   - RemakeCell
 *     - remake_cell(): rebuild cell according to its latName
 *   - RemakeCellWarnings
 *     - remake_cell(): deliver warnings when find wrong latname or cos12
 *   - PeriodicBoundaryAdjustment
 *     - periodic_boundary_adjustment(): move atoms inside the unitcell after relaxation
 *   - UpdateVel
 *     - update_vel(const ModuleBase::Vector3<double>* vel_in)
 */

class UpdateCellTest : public ::testing::Test
{
  protected:
    std::unique_ptr<UnitCell> ucell{new UnitCell};
    std::string output;
};

using UpdateCellDeathTest = UpdateCellTest;


TEST_F(UpdateCellTest, RemakeCell)
{
    std::vector<std::string> latname_in = {"sc",
                                           "fcc",
                                           "bcc",
                                           "hexagonal",
                                           "trigonal",
                                           "st",
                                           "bct",
                                           "so",
                                           "baco",
                                           "fco",
                                           "bco",
                                           "sm",
                                           "bacm",
                                           "triclinic"};
    for (int i = 0; i < latname_in.size(); ++i)
    {
        ucell->latvec.e11 = 10.0;
        ucell->latvec.e12 = 0.00;
        ucell->latvec.e13 = 0.00;
        ucell->latvec.e21 = 0.00;
        ucell->latvec.e22 = 10.0;
        ucell->latvec.e23 = 0.00;
        ucell->latvec.e31 = 0.00;
        ucell->latvec.e32 = 0.00;
        ucell->latvec.e33 = 10.0;
        ucell->latName = latname_in[i];
        unitcell::remake_cell(ucell->lat);
        if (latname_in[i] == "sc")
        {
            double celldm
                = std::sqrt(pow(ucell->latvec.e11, 2) + pow(ucell->latvec.e12, 2) + pow(ucell->latvec.e13, 2));
            EXPECT_DOUBLE_EQ(ucell->latvec.e11, celldm);
        }
        else if (latname_in[i] == "fcc")
        {
            double celldm = std::sqrt(pow(ucell->latvec.e11, 2) + pow(ucell->latvec.e12, 2) + pow(ucell->latvec.e13, 2))
                            / std::sqrt(2.0);
            EXPECT_DOUBLE_EQ(ucell->latvec.e11, -celldm);
            EXPECT_DOUBLE_EQ(ucell->latvec.e12, 0.0);
            EXPECT_DOUBLE_EQ(ucell->latvec.e13, celldm);
            EXPECT_DOUBLE_EQ(ucell->latvec.e21, 0.0);
            EXPECT_DOUBLE_EQ(ucell->latvec.e22, celldm);
            EXPECT_DOUBLE_EQ(ucell->latvec.e23, celldm);
            EXPECT_DOUBLE_EQ(ucell->latvec.e31, -celldm);
            EXPECT_DOUBLE_EQ(ucell->latvec.e32, celldm);
            EXPECT_DOUBLE_EQ(ucell->latvec.e33, 0.0);
        }
        else if (latname_in[i] == "bcc")
        {
            double celldm = std::sqrt(pow(ucell->latvec.e11, 2) + pow(ucell->latvec.e12, 2) + pow(ucell->latvec.e13, 2))
                            / std::sqrt(3.0);
            EXPECT_DOUBLE_EQ(ucell->latvec.e11, celldm);
            EXPECT_DOUBLE_EQ(ucell->latvec.e12, celldm);
            EXPECT_DOUBLE_EQ(ucell->latvec.e13, celldm);
            EXPECT_DOUBLE_EQ(ucell->latvec.e21, -celldm);
            EXPECT_DOUBLE_EQ(ucell->latvec.e22, celldm);
            EXPECT_DOUBLE_EQ(ucell->latvec.e23, celldm);
            EXPECT_DOUBLE_EQ(ucell->latvec.e31, -celldm);
            EXPECT_DOUBLE_EQ(ucell->latvec.e32, -celldm);
            EXPECT_DOUBLE_EQ(ucell->latvec.e33, celldm);
        }
        else if (latname_in[i] == "hexagonal")
        {
            double celldm1
                = std::sqrt(pow(ucell->latvec.e11, 2) + pow(ucell->latvec.e12, 2) + pow(ucell->latvec.e13, 2));
            double celldm3
                = std::sqrt(pow(ucell->latvec.e31, 2) + pow(ucell->latvec.e32, 2) + pow(ucell->latvec.e33, 2));
            double mathfoo = sqrt(3.0) / 2.0;
            EXPECT_DOUBLE_EQ(ucell->latvec.e11, celldm1);
            EXPECT_DOUBLE_EQ(ucell->latvec.e12, 0.0);
            EXPECT_DOUBLE_EQ(ucell->latvec.e13, 0.0);
            EXPECT_DOUBLE_EQ(ucell->latvec.e21, -0.5 * celldm1);
            EXPECT_DOUBLE_EQ(ucell->latvec.e22, celldm1 * mathfoo);
            EXPECT_DOUBLE_EQ(ucell->latvec.e23, 0.0);
            EXPECT_DOUBLE_EQ(ucell->latvec.e31, 0.0);
            EXPECT_DOUBLE_EQ(ucell->latvec.e32, 0.0);
            EXPECT_DOUBLE_EQ(ucell->latvec.e33, celldm3);
        }
        else if (latname_in[i] == "trigonal")
        {
            double a1 = std::sqrt(pow(ucell->latvec.e11, 2) + pow(ucell->latvec.e12, 2) + pow(ucell->latvec.e13, 2));
            double a2 = std::sqrt(pow(ucell->latvec.e21, 2) + pow(ucell->latvec.e22, 2) + pow(ucell->latvec.e23, 2));
            double a1da2 = (ucell->latvec.e11 * ucell->latvec.e21 + ucell->latvec.e12 * ucell->latvec.e22
                            + ucell->latvec.e13 * ucell->latvec.e23);
            double cosgamma = a1da2 / (a1 * a2);
            double tx = std::sqrt((1.0 - cosgamma) / 2.0);
            double ty = std::sqrt((1.0 - cosgamma) / 6.0);
            double tz = std::sqrt((1.0 + 2.0 * cosgamma) / 3.0);
            EXPECT_DOUBLE_EQ(ucell->latvec.e11, a1 * tx);
            EXPECT_DOUBLE_EQ(ucell->latvec.e12, -a1 * ty);
            EXPECT_DOUBLE_EQ(ucell->latvec.e13, a1 * tz);
            EXPECT_DOUBLE_EQ(ucell->latvec.e21, 0.0);
            EXPECT_DOUBLE_EQ(ucell->latvec.e22, 2.0 * a1 * ty);
            EXPECT_DOUBLE_EQ(ucell->latvec.e23, a1 * tz);
            EXPECT_DOUBLE_EQ(ucell->latvec.e31, -a1 * tx);
            EXPECT_DOUBLE_EQ(ucell->latvec.e32, -a1 * ty);
            EXPECT_DOUBLE_EQ(ucell->latvec.e33, a1 * tz);
        }
        else if (latname_in[i] == "st")
        {
            double a1 = std::sqrt(pow(ucell->latvec.e11, 2) + pow(ucell->latvec.e12, 2) + pow(ucell->latvec.e13, 2));
            double a3 = std::sqrt(pow(ucell->latvec.e31, 2) + pow(ucell->latvec.e32, 2) + pow(ucell->latvec.e33, 2));
            EXPECT_DOUBLE_EQ(ucell->latvec.e11, a1);
            EXPECT_DOUBLE_EQ(ucell->latvec.e12, 0.0);
            EXPECT_DOUBLE_EQ(ucell->latvec.e13, 0.0);
            EXPECT_DOUBLE_EQ(ucell->latvec.e21, 0.0);
            EXPECT_DOUBLE_EQ(ucell->latvec.e22, a1);
            EXPECT_DOUBLE_EQ(ucell->latvec.e23, 0.0);
            EXPECT_DOUBLE_EQ(ucell->latvec.e31, 0.0);
            EXPECT_DOUBLE_EQ(ucell->latvec.e32, 0.0);
            EXPECT_DOUBLE_EQ(ucell->latvec.e33, a3);
        }
        else if (latname_in[i] == "bct")
        {
            double d1 = std::abs(ucell->latvec.e11);
            double d2 = std::abs(ucell->latvec.e13);
            EXPECT_DOUBLE_EQ(ucell->latvec.e11, d1);
            EXPECT_DOUBLE_EQ(ucell->latvec.e12, -d1);
            EXPECT_DOUBLE_EQ(ucell->latvec.e13, d2);
            EXPECT_DOUBLE_EQ(ucell->latvec.e21, d1);
            EXPECT_DOUBLE_EQ(ucell->latvec.e22, d1);
            EXPECT_DOUBLE_EQ(ucell->latvec.e23, d2);
            EXPECT_DOUBLE_EQ(ucell->latvec.e31, -d1);
            EXPECT_DOUBLE_EQ(ucell->latvec.e32, -d1);
            EXPECT_DOUBLE_EQ(ucell->latvec.e33, d2);
        }
        else if (latname_in[i] == "so")
        {
            double a1 = std::sqrt(pow(ucell->latvec.e11, 2) + pow(ucell->latvec.e12, 2) + pow(ucell->latvec.e13, 2));
            double a2 = std::sqrt(pow(ucell->latvec.e21, 2) + pow(ucell->latvec.e22, 2) + pow(ucell->latvec.e23, 2));
            double a3 = std::sqrt(pow(ucell->latvec.e31, 2) + pow(ucell->latvec.e32, 2) + pow(ucell->latvec.e33, 2));
            EXPECT_DOUBLE_EQ(ucell->latvec.e11, a1);
            EXPECT_DOUBLE_EQ(ucell->latvec.e12, 0.0);
            EXPECT_DOUBLE_EQ(ucell->latvec.e13, 0.0);
            EXPECT_DOUBLE_EQ(ucell->latvec.e21, 0.0);
            EXPECT_DOUBLE_EQ(ucell->latvec.e22, a2);
            EXPECT_DOUBLE_EQ(ucell->latvec.e23, 0.0);
            EXPECT_DOUBLE_EQ(ucell->latvec.e31, 0.0);
            EXPECT_DOUBLE_EQ(ucell->latvec.e32, 0.0);
            EXPECT_DOUBLE_EQ(ucell->latvec.e33, a3);
        }
        else if (latname_in[i] == "baco")
        {
            double d1 = std::abs(ucell->latvec.e11);
            double d2 = std::abs(ucell->latvec.e22);
            double d3 = std::abs(ucell->latvec.e33);
            EXPECT_DOUBLE_EQ(ucell->latvec.e11, d1);
            EXPECT_DOUBLE_EQ(ucell->latvec.e12, d2);
            EXPECT_DOUBLE_EQ(ucell->latvec.e13, 0.0);
            EXPECT_DOUBLE_EQ(ucell->latvec.e21, -d1);
            EXPECT_DOUBLE_EQ(ucell->latvec.e22, d2);
            EXPECT_DOUBLE_EQ(ucell->latvec.e23, 0.0);
            EXPECT_DOUBLE_EQ(ucell->latvec.e31, 0.0);
            EXPECT_DOUBLE_EQ(ucell->latvec.e32, 0.0);
            EXPECT_DOUBLE_EQ(ucell->latvec.e33, d3);
        }
        else if (latname_in[i] == "fco")
        {
            double d1 = std::abs(ucell->latvec.e11);
            double d2 = std::abs(ucell->latvec.e22);
            double d3 = std::abs(ucell->latvec.e33);
            EXPECT_DOUBLE_EQ(ucell->latvec.e11, d1);
            EXPECT_DOUBLE_EQ(ucell->latvec.e12, 0.0);
            EXPECT_DOUBLE_EQ(ucell->latvec.e13, d3);
            EXPECT_DOUBLE_EQ(ucell->latvec.e21, d1);
            EXPECT_DOUBLE_EQ(ucell->latvec.e22, d2);
            EXPECT_DOUBLE_EQ(ucell->latvec.e23, 0.0);
            EXPECT_DOUBLE_EQ(ucell->latvec.e31, 0.0);
            EXPECT_DOUBLE_EQ(ucell->latvec.e32, d2);
            EXPECT_DOUBLE_EQ(ucell->latvec.e33, d3);
        }
        else if (latname_in[i] == "bco")
        {
            double d1 = std::abs(ucell->latvec.e11);
            double d2 = std::abs(ucell->latvec.e22);
            double d3 = std::abs(ucell->latvec.e33);
            EXPECT_DOUBLE_EQ(ucell->latvec.e11, d1);
            EXPECT_DOUBLE_EQ(ucell->latvec.e12, d2);
            EXPECT_DOUBLE_EQ(ucell->latvec.e13, d3);
            EXPECT_DOUBLE_EQ(ucell->latvec.e21, -d1);
            EXPECT_DOUBLE_EQ(ucell->latvec.e22, d2);
            EXPECT_DOUBLE_EQ(ucell->latvec.e23, d3);
            EXPECT_DOUBLE_EQ(ucell->latvec.e31, -d1);
            EXPECT_DOUBLE_EQ(ucell->latvec.e32, -d2);
            EXPECT_DOUBLE_EQ(ucell->latvec.e33, d3);
        }
        else if (latname_in[i] == "sm")
        {
            double a1 = std::sqrt(pow(ucell->latvec.e11, 2) + pow(ucell->latvec.e12, 2) + pow(ucell->latvec.e13, 2));
            double a2 = std::sqrt(pow(ucell->latvec.e21, 2) + pow(ucell->latvec.e22, 2) + pow(ucell->latvec.e23, 2));
            double a3 = std::sqrt(pow(ucell->latvec.e31, 2) + pow(ucell->latvec.e32, 2) + pow(ucell->latvec.e33, 2));
            double a1da2 = (ucell->latvec.e11 * ucell->latvec.e21 + ucell->latvec.e12 * ucell->latvec.e22
                            + ucell->latvec.e13 * ucell->latvec.e23);
            double cosgamma = a1da2 / (a1 * a2);
            double d1 = a2 * cosgamma;
            double d2 = a2 * std::sqrt(1.0 - cosgamma * cosgamma);
            EXPECT_DOUBLE_EQ(ucell->latvec.e11, a1);
            EXPECT_DOUBLE_EQ(ucell->latvec.e12, 0.0);
            EXPECT_DOUBLE_EQ(ucell->latvec.e13, 0.0);
            EXPECT_DOUBLE_EQ(ucell->latvec.e21, d1);
            EXPECT_DOUBLE_EQ(ucell->latvec.e22, d2);
            EXPECT_DOUBLE_EQ(ucell->latvec.e23, 0.0);
            EXPECT_DOUBLE_EQ(ucell->latvec.e31, 0.0);
            EXPECT_DOUBLE_EQ(ucell->latvec.e32, 0.0);
            EXPECT_DOUBLE_EQ(ucell->latvec.e33, a3);
        }
        else if (latname_in[i] == "bacm")
        {
            double d1 = std::abs(ucell->latvec.e11);
            double a2 = std::sqrt(pow(ucell->latvec.e21, 2) + pow(ucell->latvec.e22, 2) + pow(ucell->latvec.e23, 2));
            double d3 = std::abs(ucell->latvec.e13);
            double cosgamma = ucell->latvec.e21 / a2;
            double f1 = a2 * cosgamma;
            double f2 = a2 * std::sqrt(1.0 - cosgamma * cosgamma);
            EXPECT_DOUBLE_EQ(ucell->latvec.e11, d1);
            EXPECT_DOUBLE_EQ(ucell->latvec.e12, 0.0);
            EXPECT_DOUBLE_EQ(ucell->latvec.e13, -d3);
            EXPECT_DOUBLE_EQ(ucell->latvec.e21, f1);
            EXPECT_DOUBLE_EQ(ucell->latvec.e22, f2);
            EXPECT_DOUBLE_EQ(ucell->latvec.e23, 0.0);
            EXPECT_DOUBLE_EQ(ucell->latvec.e31, d1);
            EXPECT_DOUBLE_EQ(ucell->latvec.e32, 0.0);
            EXPECT_DOUBLE_EQ(ucell->latvec.e33, d3);
        }
        else if (latname_in[i] == "triclinic")
        {
            double a1 = std::sqrt(pow(ucell->latvec.e11, 2) + pow(ucell->latvec.e12, 2) + pow(ucell->latvec.e13, 2));
            double a2 = std::sqrt(pow(ucell->latvec.e21, 2) + pow(ucell->latvec.e22, 2) + pow(ucell->latvec.e23, 2));
            double a3 = std::sqrt(pow(ucell->latvec.e31, 2) + pow(ucell->latvec.e32, 2) + pow(ucell->latvec.e33, 2));
            double a1da2 = (ucell->latvec.e11 * ucell->latvec.e21 + ucell->latvec.e12 * ucell->latvec.e22
                            + ucell->latvec.e13 * ucell->latvec.e23);
            double a1da3 = (ucell->latvec.e11 * ucell->latvec.e31 + ucell->latvec.e12 * ucell->latvec.e32
                            + ucell->latvec.e13 * ucell->latvec.e33);
            double a2da3 = (ucell->latvec.e21 * ucell->latvec.e31 + ucell->latvec.e22 * ucell->latvec.e32
                            + ucell->latvec.e23 * ucell->latvec.e33);
            double cosgamma = a1da2 / a1 / a2;
            double singamma = std::sqrt(1.0 - cosgamma * cosgamma);
            double cosbeta = a1da3 / a1 / a3;
            double cosalpha = a2da3 / a2 / a3;
            double d1 = std::sqrt(1.0 + 2.0 * cosgamma * cosbeta * cosalpha - cosgamma * cosgamma - cosbeta * cosbeta
                                  - cosalpha * cosalpha)
                        / singamma;
            EXPECT_DOUBLE_EQ(ucell->latvec.e11, a1);
            EXPECT_DOUBLE_EQ(ucell->latvec.e12, 0.0);
            EXPECT_DOUBLE_EQ(ucell->latvec.e13, 0.0);
            EXPECT_DOUBLE_EQ(ucell->latvec.e21, a2 * cosgamma);
            EXPECT_DOUBLE_EQ(ucell->latvec.e22, a2 * singamma);
            EXPECT_DOUBLE_EQ(ucell->latvec.e23, 0.0);
            EXPECT_DOUBLE_EQ(ucell->latvec.e31, a3 * cosbeta);
            EXPECT_DOUBLE_EQ(ucell->latvec.e32, a3 * (cosalpha - cosbeta * cosgamma) / singamma);
            EXPECT_DOUBLE_EQ(ucell->latvec.e33, a3 * d1);
        }
    }
}

TEST_F(UpdateCellDeathTest, RemakeCellWarnings)
{
    std::vector<std::string> latname_in = {"user_defined_lattice", "trigonal", "bacm", "triclinic", "arbitrary"};
    for (int i = 0; i < latname_in.size(); ++i)
    {
        ucell->latvec.e11 = 10.0;
        ucell->latvec.e12 = 0.00;
        ucell->latvec.e13 = 0.00;
        ucell->latvec.e21 = 10.0;
        ucell->latvec.e22 = 0.00;
        ucell->latvec.e23 = 0.00;
        ucell->latvec.e31 = 0.00;
        ucell->latvec.e32 = 0.00;
        ucell->latvec.e33 = 10.0;
        ucell->latName = latname_in[i];
        testing::internal::CaptureStdout();
        EXPECT_EXIT(unitcell::remake_cell(ucell->lat), ::testing::ExitedWithCode(1), "");
        std::string output = testing::internal::GetCapturedStdout();
        if (latname_in[i] == "user_defined_lattice")
        {
            EXPECT_THAT(output, testing::HasSubstr("to use fixed_ibrav, latname must be provided"));
        }
        else if (latname_in[i] == "trigonal" || latname_in[i] == "bacm" || latname_in[i] == "triclinic")
        {
            EXPECT_THAT(output, testing::HasSubstr("wrong cos12!"));
        }
        else
        {
            EXPECT_THAT(output, testing::HasSubstr("latname type not supported!"));
        }
    }
}

// mohan comment out 2025-07-14
/*
TEST_F(UpdateCellDeathTest, PeriodicBoundaryAdjustment1)
{
    UcellTestPrepare utp = UcellTestLib["C1H2-PBA"];
    ucell = utp.SetUcellInfo();
    testing::internal::CaptureStdout();
    EXPECT_EXIT(unitcell::periodic_boundary_adjustment(
                ucell->atoms,ucell->latvec,ucell->ntype),
                ::testing::ExitedWithCode(1), "");
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("the movement of atom is larger than the length of cell"));
}
*/

TEST_F(UpdateCellTest, PeriodicBoundaryAdjustment2)
{
    UcellTestPrepare utp = UcellTestLib["C1H2-Index"];
    ucell = utp.SetUcellInfo();
    EXPECT_NO_THROW(unitcell::periodic_boundary_adjustment(
                    ucell->atoms,ucell->latvec,ucell->ntype));
}

TEST_F(UpdateCellTest, UpdateVel)
{
    UcellTestPrepare utp = UcellTestLib["C1H2-Index"];
    ucell = utp.SetUcellInfo();
    ModuleBase::Vector3<double>* vel_in = new ModuleBase::Vector3<double>[ucell->nat];
    for (int iat = 0; iat < ucell->nat; ++iat)
    {
        vel_in[iat].set(iat * 0.1, iat * 0.1, iat * 0.1);
    }
    unitcell::update_vel(vel_in,ucell->ntype,ucell->nat,ucell->atoms);
    for (int iat = 0; iat < ucell->nat; ++iat)
    {
        EXPECT_DOUBLE_EQ(vel_in[iat].x, 0.1 * iat);
        EXPECT_DOUBLE_EQ(vel_in[iat].y, 0.1 * iat);
        EXPECT_DOUBLE_EQ(vel_in[iat].z, 0.1 * iat);
    }
    delete[] vel_in;
}

