#include "gtest/gtest.h"

#include <cstdio>
#include <fstream>
#include <string>

#include "source_cell/atom_spec.h"
#include "source_cell/read_orb.h"

/************************************************
 *  unit test of read_orb.cpp boundary handling
 ***********************************************/

/**
 * - Tested Functions:
 *   - ZeroCountBeyondLmaxTolerated
 *     - read_orb_file(): orbital generators conventionally keep zero-count
 *       declaration lines for channels beyond Lmax (e.g. an s-only minimal
 *       basis still lists "Number of Porbital--> 0"); these must be tolerated
 *       without accessing l_nchi out of bounds
 *   - NormalMultiChannel
 *     - read_orb_file(): channels within Lmax are counted as before
 *   - NonzeroCountBeyondLmaxRejected
 *     - read_orb_file(): a nonzero count beyond Lmax is a malformed file
 *   - ChannelBeforeLmaxRejected
 *     - read_orb_file(): channel lines before the Lmax declaration are
 *       malformed, because l_nchi is not allocated yet
 *   - NonNumericCountBeyondLmaxRejected / MissingCountBeyondLmaxRejected
 *     - read_orb_file(): a count that fails numeric extraction (garbage or
 *       missing) is malformed; READ_VALUE cannot be used for it because it
 *       hides the failure
 */

class ReadOrbBoundaryTest : public ::testing::Test
{
  protected:
    Atom atom;
    std::ofstream ofs_running;
    std::string log_path;
    std::string orb_path;

    void SetUp() override
    {
        log_path = testing::TempDir() + "read_orb_boundary.log";
        ofs_running.open(log_path);
    }

    void TearDown() override
    {
        ofs_running.close();
        std::remove(log_path.c_str());
        if (!orb_path.empty())
        {
            std::remove(orb_path.c_str());
        }
    }

    std::string write_orb_file(const std::string& name, const std::string& content)
    {
        orb_path = testing::TempDir() + name;
        std::ofstream ofs(orb_path);
        ofs << content;
        ofs.close();
        ofs_running << "created test ORBITAL file: " << orb_path << std::endl;
        return orb_path;
    }
};

TEST_F(ReadOrbBoundaryTest, ZeroCountBeyondLmaxTolerated)
{
    std::string orb = write_orb_file("H_1s_minimal.orb",
                                           "Element                     H\n"
                                           "Energy Cutoff(Ry)           100\n"
                                           "Radius Cutoff(a.u.)         8\n"
                                           "Lmax                        0\n"
                                           "Number of Sorbital-->       1\n"
                                           "Number of Porbital-->       0\n");
    const bool result = unitcell::read_orb_file(0, orb, ofs_running, &atom);
    EXPECT_TRUE(result);
    ASSERT_EQ(atom.l_nchi.size(), 1);
    EXPECT_EQ(atom.l_nchi[0], 1);
    EXPECT_EQ(atom.nw, 1);
}

TEST_F(ReadOrbBoundaryTest, NormalMultiChannel)
{
    std::string orb = write_orb_file("C_spd.orb",
                                           "Element                     C\n"
                                           "Lmax                        2\n"
                                           "Cutoff(a.u.)                8\n"
                                           "Number of Sorbital-->       2\n"
                                           "Number of Porbital-->       1\n"
                                           "Number of Dorbital-->       1\n"
                                           "SUMMARY  END\n");
    const bool result = unitcell::read_orb_file(0, orb, ofs_running, &atom);
    EXPECT_TRUE(result);
    ASSERT_EQ(atom.l_nchi.size(), 3);
    EXPECT_EQ(atom.l_nchi[0], 2);
    EXPECT_EQ(atom.l_nchi[1], 1);
    EXPECT_EQ(atom.l_nchi[2], 1);
    EXPECT_EQ(atom.nw, 2 * 1 + 1 * 3 + 1 * 5);
}

TEST_F(ReadOrbBoundaryTest, NonzeroCountBeyondLmaxRejected)
{
    std::string orb = write_orb_file("H_bad_count.orb",
                                           "Element                     H\n"
                                           "Lmax                        0\n"
                                           "Number of Sorbital-->       1\n"
                                           "Number of Porbital-->       1\n");
    const bool result = unitcell::read_orb_file(0, orb, ofs_running, &atom);
    EXPECT_FALSE(result);
}

TEST_F(ReadOrbBoundaryTest, ChannelBeforeLmaxRejected)
{
    std::string orb = write_orb_file("H_no_lmax.orb",
                                           "Element                     H\n"
                                           "Number of Sorbital-->       1\n"
                                           "Lmax                        0\n");
    const bool result = unitcell::read_orb_file(0, orb, ofs_running, &atom);
    EXPECT_FALSE(result);
}

TEST_F(ReadOrbBoundaryTest, NonNumericCountBeyondLmaxRejected)
{
    std::string orb = write_orb_file("H_bad_token.orb",
                                           "Element                     H\n"
                                           "Lmax                        0\n"
                                           "Number of Sorbital-->       1\n"
                                           "Number of Porbital-->       garbage\n");
    const bool result = unitcell::read_orb_file(0, orb, ofs_running, &atom);
    EXPECT_FALSE(result);
}

TEST_F(ReadOrbBoundaryTest, MissingCountBeyondLmaxRejected)
{
    std::string orb = write_orb_file("H_missing_count.orb",
                                           "Element                     H\n"
                                           "Lmax                        0\n"
                                           "Number of Sorbital-->       1\n"
                                           "Number of Porbital-->\n"
                                           "SUMMARY  END\n");
    const bool result = unitcell::read_orb_file(0, orb, ofs_running, &atom);
    EXPECT_FALSE(result);
}
