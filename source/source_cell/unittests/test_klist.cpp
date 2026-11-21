#include "gtest/gtest.h"

#include <cstdio>
#include <fstream>
#include <string>

#include "source_cell/klist.h"
#include "source_cell/qlist.h"
#include "source_cell/unitcell.h"

/************************************************
 *  unit test of klist.cpp / qlist.cpp input validation
 ***********************************************/

/**
 * - Tested Functions:
 *   - NegativeTotalCountRejected
 *     - K_Vectors::parse_kfile(): a negative k-point count must be rejected
 *       instead of falling through both the nkstot==0 and nkstot>0 branches
 *   - ZeroMpDimensionRejected / NonNumericMpDimensionRejected
 *     / NegativeMpDimensionRejected / TwoNegativeMpDimensionsRejected
 *     / MpGridProductOverflowRejected / MpGridProductOverflowInt64Rejected
 *     - K_Vectors::read_mp_mesh(): malformed Monkhorst-Pack grid dimensions
 *       (zero, non-numeric, negative, or overflowing int) must be rejected;
 *       a negative dimension previously crashed in renew() with an uncaught
 *       std::length_error, and a negative pair has a positive product;
 *       2097152^3 = 2^63 overflows even a 64-bit product, so the guard must
 *       compare before multiplying
 *   - ValidMpGridAccepted
 *     - a well-formed 2x2x2 Gamma grid still works (8 k-points)
 *   - QptNegativeMpDimensionRejected / QptMpGridProductOverflowInt64Rejected
 *     - QList::read_from_file(): a malformed MP grid resets the q-point count
 *       following its existing warning-and-reset error style
 */

class KptValidationTest : public ::testing::Test
{
  protected:
    UnitCell ucell;
    std::ofstream ofs_running;
    std::ofstream ofs_warning;
    std::string log_path;
    std::string warn_path;
    std::string kpt_path;

    void SetUp() override
    {
        log_path = testing::TempDir() + "kpt_validation.log";
        warn_path = testing::TempDir() + "kpt_validation.warn";
        ofs_running.open(log_path);
        ofs_warning.open(warn_path);
    }

    void TearDown() override
    {
        ofs_running.close();
        ofs_warning.close();
        std::remove(log_path.c_str());
        std::remove(warn_path.c_str());
        if (!kpt_path.empty())
        {
            std::remove(kpt_path.c_str());
        }
    }

    std::string write_kpt_file(const std::string& name, const std::string& content)
    {
        kpt_path = testing::TempDir() + name;
        std::ofstream ofs(kpt_path);
        ofs << content;
        ofs.close();
        ofs_running << "created test KPT file: " << kpt_path << std::endl;
        return kpt_path;
    }

    bool read_kpt(const std::string& content)
    {
        std::string kpt = write_kpt_file("case.kpt", content);
        K_Vectors kv;
        const double kspacing[3] = {-1.0, -1.0, -1.0};
        const double koffset[3] = {0.0, 0.0, 0.0};
        return kv.read_kpoints_for_testing(ucell,
                                           kpt,
                                           false,
                                           kspacing,
                                           "gamma",
                                           koffset,
                                           ofs_running,
                                           ofs_warning,
                                           0);
    }
};

TEST_F(KptValidationTest, NegativeTotalCountRejected)
{
    const bool result = read_kpt("K_POINTS\n-4\nGamma\n4 4 4 0 0 0\n");
    EXPECT_FALSE(result);
}

TEST_F(KptValidationTest, ZeroMpDimensionRejected)
{
    const bool result = read_kpt("K_POINTS\n0\nGamma\n0 4 4 0 0 0\n");
    EXPECT_FALSE(result);
}

TEST_F(KptValidationTest, NonNumericMpDimensionRejected)
{
    const bool result = read_kpt("K_POINTS\n0\nGamma\nx 4 4 0 0 0\n");
    EXPECT_FALSE(result);
}

TEST_F(KptValidationTest, NegativeMpDimensionRejected)
{
    const bool result = read_kpt("K_POINTS\n0\nGamma\n-1 4 4 0 0 0\n");
    EXPECT_FALSE(result);
}

TEST_F(KptValidationTest, TwoNegativeMpDimensionsRejected)
{
    // a negative pair has a positive product; each dimension must be checked
    const bool result = read_kpt("K_POINTS\n0\nGamma\n-1 -1 4 0 0 0\n");
    EXPECT_FALSE(result);
}

TEST_F(KptValidationTest, MpGridProductOverflowRejected)
{
    const bool result = read_kpt("K_POINTS\n0\nGamma\n2000 2000 2000 0 0 0\n");
    EXPECT_FALSE(result);
}

TEST_F(KptValidationTest, MpGridProductOverflowInt64Rejected)
{
    // 2097152^3 = 2^63 overflows even a 64-bit signed product, so the guard
    // must reject by comparing before multiplying
    const bool result = read_kpt("K_POINTS\n0\nGamma\n2097152 2097152 2097152 0 0 0\n");
    EXPECT_FALSE(result);
}

TEST_F(KptValidationTest, ValidMpGridAccepted)
{
    const bool result = read_kpt("K_POINTS\n0\nGamma\n2 2 2 0 0 0\n");
    EXPECT_TRUE(result);
}

TEST_F(KptValidationTest, QptNegativeMpDimensionRejected)
{
    // discriminating case: the old code crashed on a negative dimension
    // (renew(-4) threw), the new one resets the q-point count instead
    std::string qpt = write_kpt_file("case.qpt", "Q_POINTS\n0\nGamma\n-1 4 4 0 0 0\n");
    ModuleCell::QList ql;
    ql.read_from_file(qpt, ucell);
    EXPECT_EQ(ql.nkstot, 0);
}

TEST_F(KptValidationTest, QptMpGridProductOverflowInt64Rejected)
{
    // 2097152^3 = 2^63 overflows even a 64-bit signed product, so the guard
    // must reject by comparing before multiplying
    std::string qpt = write_kpt_file("case.qpt", "Q_POINTS\n0\nGamma\n2097152 2097152 2097152 0 0 0\n");
    ModuleCell::QList ql;
    ql.read_from_file(qpt, ucell);
    EXPECT_EQ(ql.nkstot, 0);
}
