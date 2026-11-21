#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cstdio>
#include <fstream>
#include <string>

#include "source_cell/basecell.h"
#include "source_cell/read_stru.h"

/************************************************
 *  unit test of read_atom_species.cpp lattice validation
 ***********************************************/

/**
 * - Tested Functions:
 *   - HexagonalMissingParamsRejected / TrigonalMissingParamsRejected
 *     - read_lattice_constant(): non-cubic named lattices require
 *       LATTICE_PARAMETERS; without it the shape variables stayed zero and
 *       produced a zero-volume or NaN cell (issue #8354)
 *   - TrigonalMalformedParamRejected
 *     - a failed cos(alpha) extraction must not silently become the
 *       simple-cubic (alpha = 90 degree) special case
 *   - ValidHexagonalWithParams
 *     - a hexagonal lattice with c/a given still reads normally
 *   - CubicScNeedsNoParams
 *     - cubic lattice names do not need LATTICE_PARAMETERS
 *   - DegenerateLatticeVectorsRejected
 *     - explicit LATTICE_VECTORS forming a zero-volume cell are rejected
 */

class ReadLatticeValidationTest : public ::testing::Test
{
  protected:
    Lattice lat;
    std::ofstream ofs_running;
    std::string log_path;
    std::string stru_path;

    void SetUp() override
    {
        log_path = testing::TempDir() + "read_lattice_validation.log";
        ofs_running.open(log_path);
    }

    void TearDown() override
    {
        ofs_running.close();
        std::remove(log_path.c_str());
        if (!stru_path.empty())
        {
            std::remove(stru_path.c_str());
        }
    }

    std::string write_stru_file(const std::string& name, const std::string& content)
    {
        stru_path = testing::TempDir() + name;
        std::ofstream ofs(stru_path);
        ofs << content;
        ofs.close();
        ofs_running << "created test STRU file: " << stru_path << std::endl;
        return stru_path;
    }
};

TEST_F(ReadLatticeValidationTest, HexagonalMissingParamsRejected)
{
    const std::string stru = write_stru_file("hex_no_param.stru", "LATTICE_CONSTANT\n4.0\n");
    lat.latName = "hexagonal";
    std::ifstream ifa(stru);
    testing::internal::CaptureStdout();
    EXPECT_EXIT(unitcell::read_lattice_constant(ifa, ofs_running, lat), ::testing::ExitedWithCode(1), "");
    const std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("LATTICE_PARAMETERS is required for lattice type hexagonal"));
    ifa.close();
}

TEST_F(ReadLatticeValidationTest, TrigonalMissingParamsRejected)
{
    const std::string stru = write_stru_file("trigonal_no_param.stru", "LATTICE_CONSTANT\n4.0\n");
    lat.latName = "trigonal";
    std::ifstream ifa(stru);
    testing::internal::CaptureStdout();
    EXPECT_EXIT(unitcell::read_lattice_constant(ifa, ofs_running, lat), ::testing::ExitedWithCode(1), "");
    const std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("LATTICE_PARAMETERS is required for lattice type trigonal"));
    ifa.close();
}

TEST_F(ReadLatticeValidationTest, TrigonalMalformedParamRejected)
{
    // a failed cos(alpha) read must not silently become the simple-cubic
    // (alpha = 90 degree) special case
    const std::string stru =
        write_stru_file("trigonal_bad_param.stru", "LATTICE_CONSTANT\n4.0\nLATTICE_PARAMETERS\nabc\n");
    lat.latName = "trigonal";
    std::ifstream ifa(stru);
    testing::internal::CaptureStdout();
    EXPECT_EXIT(unitcell::read_lattice_constant(ifa, ofs_running, lat), ::testing::ExitedWithCode(1), "");
    const std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("invalid cos(alpha)"));
    ifa.close();
}

TEST_F(ReadLatticeValidationTest, ValidHexagonalWithParams)
{
    const std::string stru =
        write_stru_file("hex_with_param.stru", "LATTICE_CONSTANT\n4.0\nLATTICE_PARAMETERS\n1.633\n");
    lat.latName = "hexagonal";
    std::ifstream ifa(stru);
    EXPECT_TRUE(unitcell::read_lattice_constant(ifa, ofs_running, lat));
    EXPECT_DOUBLE_EQ(lat.latvec.e33, 1.633);
    EXPECT_GT(lat.latvec.Det(), 0.0);
    ifa.close();
}

TEST_F(ReadLatticeValidationTest, CubicScNeedsNoParams)
{
    const std::string stru = write_stru_file("sc_no_param.stru", "LATTICE_CONSTANT\n4.0\n");
    lat.latName = "sc";
    std::ifstream ifa(stru);
    EXPECT_TRUE(unitcell::read_lattice_constant(ifa, ofs_running, lat));
    EXPECT_DOUBLE_EQ(lat.latvec.Det(), 1.0);
    ifa.close();
}

TEST_F(ReadLatticeValidationTest, DegenerateLatticeVectorsRejected)
{
    const std::string stru = write_stru_file("degenerate.stru",
                                             "LATTICE_CONSTANT\n4.0\n"
                                             "LATTICE_VECTORS\n"
                                             "1 0 0\n"
                                             "0 1 0\n"
                                             "0 0 0\n");
    lat.latName = "user_defined_lattice";
    std::ifstream ifa(stru);
    testing::internal::CaptureStdout();
    EXPECT_EXIT(unitcell::read_lattice_constant(ifa, ofs_running, lat),
                ::testing::ExitedWithCode(1),
                "");
    const std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("zero-volume"));
    ifa.close();
}
