#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cstdio>
#include <fstream>
#include <memory>
#include <string>
#include <valarray>
#include <vector>

#include "source_base/mathzone.h"
#include "source_cell/unitcell.h"
#include "source_cell/read_orb.h"
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
 *  unit test of read_orb.cpp
 ***********************************************/

/**
 * - Tested Functions:
 *   - ReadOrbFile
 *     - read_orb_file(): read header part of orbital file
 *   - ReadOrbFileWarning
 *     - read_orb_file(): ABACUS cannot find the ORBITAL file
 */

class ReadOrbTest : public ::testing::Test
{
  protected:
    std::unique_ptr<UnitCell> ucell;
    std::string output;
};

#ifdef __LCAO
TEST_F(ReadOrbTest, ReadOrbFile)
{
    UcellTestPrepare utp = UcellTestLib["C1H2-Read"];
    ucell = utp.SetUcellInfo();
    std::string orb_file = "./support/C.orb";
    std::ofstream ofs_running;
    ofs_running.open("tmp_readorbfile");
    bool result = unitcell::read_orb_file(0, orb_file, ofs_running, &(ucell->atoms[0]));
    ofs_running << " result=" << result << std::endl;
    EXPECT_TRUE(result);
    ofs_running.close();
    EXPECT_EQ(ucell->atoms[0].nw, 25);
    remove("tmp_readorbfile");
}
#endif

TEST_F(ReadOrbTest, ReadOrbFileWarning)
{
    UcellTestPrepare utp = UcellTestLib["C1H2-Read"];
    ucell = utp.SetUcellInfo();
    std::string orb_file = "./support/CC.orb";
    std::ofstream ofs_running;
    ofs_running.open("tmp_readorbfilewarning");
    testing::internal::CaptureStdout();
    bool result = unitcell::read_orb_file(0, orb_file, ofs_running, &(ucell->atoms[0]));
    output = testing::internal::GetCapturedStdout();
    ofs_running << output << std::endl;
    EXPECT_FALSE(result);
    EXPECT_THAT(output, testing::HasSubstr("Element index 1"));
    EXPECT_THAT(output, testing::HasSubstr("orbital file: ./support/CC.orb"));
    ofs_running.close();
    remove("tmp_readorbfilewarning");
}
