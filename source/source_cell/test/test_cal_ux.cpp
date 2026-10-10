#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <memory>
#include <valarray>
#include <vector>

#include "source_base/global_variable.h"
#include "source_base/mathzone.h"
#include "source_cell/unitcell.h"
#include "source_cell/cal_ux.h"
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
 *  unit test of cal_ux.cpp
 ***********************************************/

/**
 * - Tested Functions:
 *   - JudgeParallel
 *     - judge_parallel: judge if two vectors a[3] and Vector3<double> b are parallel
 *   - CalUx1
 *     - cal_ux: non-parallel initial moments, no common direction
 *   - CalUx2
 *     - cal_ux: parallel moments, normalized common direction (1,1,1)/sqrt(3)
 */

class CalUxTest : public ::testing::Test
{
  protected:
    std::unique_ptr<UnitCell> ucell;
};

TEST_F(CalUxTest, JudgeParallel)
{
    ModuleBase::Vector3<double> b(1.0, 1.0, 1.0);
    double a[3] = {1.0, 1.0, 1.0};
    EXPECT_TRUE(unitcell::judge_parallel(a, b));

    // the negative case, moved here from MagnetismTest.JudgeParallel when the
    // duplicate Magnetism::judge_parallel was deleted
    double c[3] = {1.0, 0.0, 0.0};
    ModuleBase::Vector3<double> d(0.0, 1.0, 0.0);
    EXPECT_FALSE(unitcell::judge_parallel(c, d));
}

TEST_F(CalUxTest, CalUx1)
{
    UcellTestPrepare utp = UcellTestLib["C1H2-Read"];
    ucell = utp.SetUcellInfo();
    ucell->atoms[0].m_loc_[0].set(0, -1, 0);
    ucell->atoms[1].m_loc_[0].set(1, 1, 1);
    ucell->atoms[1].m_loc_[1].set(0, 0, 0);
    const int nspin = 4;
    unitcell::cal_ux(*ucell, nspin);
    EXPECT_FALSE(ucell->magnet.lsign_);
    EXPECT_DOUBLE_EQ(ucell->magnet.ux_[0], 0);
    EXPECT_DOUBLE_EQ(ucell->magnet.ux_[1], -1);
    EXPECT_DOUBLE_EQ(ucell->magnet.ux_[2], 0);
}

TEST_F(CalUxTest, CalUx2)
{
    UcellTestPrepare utp = UcellTestLib["C1H2-Read"];
    ucell = utp.SetUcellInfo();
    ucell->atoms[0].m_loc_[0].set(0, 0, 0);
    ucell->atoms[1].m_loc_[0].set(1, 1, 1);
    ucell->atoms[1].m_loc_[1].set(0, 0, 0);
    //(0,0,0) is also parallel to (1,1,1)
    const int nspin = 4;
    unitcell::cal_ux(*ucell, nspin);
    EXPECT_TRUE(ucell->magnet.lsign_);
    EXPECT_NEAR(ucell->magnet.ux_[0], 0.57735, 1e-5);
    EXPECT_NEAR(ucell->magnet.ux_[1], 0.57735, 1e-5);
    EXPECT_NEAR(ucell->magnet.ux_[2], 0.57735, 1e-5);
}
