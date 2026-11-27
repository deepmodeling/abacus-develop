#include "source_io/module_output/cube_io.h"

#include <fstream>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "source_base/global_variable.h"
#include "source_base/parallel_grid.h"
#include "source_io/module_parameter/parameter.h"
#include "prepare_unitcell.h"

Magnetism::Magnetism()
{
    this->tot_mag = 0.0;
    this->abs_mag = 0.0;
}

Magnetism::~Magnetism()
{
}

/***************************************************************
 *  unit test of malformed-cube-file validation in read_cube
 *  and read_vdata_palgrid (issue #7563)
 ***************************************************************/

/**
 * - Tested Functions:
 *   - read_cube()
 *     - returns false for a truncated file (missing grid data)
 *     - returns false for invalid (zero) grid dimensions
 *   - read_vdata_palgrid()
 *     - aborts the run instead of reading invalid data or hanging other ranks
 */

struct ReadCubeInvalidTest : public ::testing::Test
{
    std::vector<std::string> comment;
    int natom = 0;
    std::vector<double> origin;
    int nx_read = 0;
    int ny_read = 0;
    int nz_read = 0;
    std::vector<double> dx;
    std::vector<double> dy;
    std::vector<double> dz;
    std::vector<int> atom_type;
    std::vector<double> atom_charge;
    std::vector<std::vector<double>> atom_pos;
    std::vector<double> data_read;

    bool call_read_cube(const std::string& fn)
    {
        return ModuleIO::read_cube(fn, comment, natom, origin,
                                   nx_read, ny_read, nz_read,
                                   dx, dy, dz,
                                   atom_type, atom_charge, atom_pos, data_read);
    }
};

TEST_F(ReadCubeInvalidTest, ValidFile)
{
    EXPECT_TRUE(call_read_cube("./support/chg.cube"));
    EXPECT_EQ(natom, 2);
    EXPECT_EQ(nx_read, 36);
    EXPECT_EQ(ny_read, 36);
    EXPECT_EQ(nz_read, 36);
    EXPECT_EQ(data_read.size(), 36 * 36 * 36);
}

TEST_F(ReadCubeInvalidTest, TruncatedData)
{
    // valid header but the grid data is incomplete
    const std::string fn = "test_cube_truncated.cube";
    std::ofstream ofs(fn);
    ofs << "comment line 1\n";
    ofs << "comment line 2\n";
    ofs << "1 0.0 0.0 0.0\n";
    ofs << "2 1.0 0.0 0.0\n";
    ofs << "2 0.0 1.0 0.0\n";
    ofs << "2 0.0 0.0 1.0\n";
    ofs << "14 4.0 0.0 0.0 0.0\n";
    ofs << "1.0 2.0 3.0\n"; // only 3 of the 8 expected values
    ofs.close();

    EXPECT_FALSE(call_read_cube(fn));
    std::remove(fn.c_str());
}

TEST_F(ReadCubeInvalidTest, InvalidDimensions)
{
    // zero grid dimension is invalid
    const std::string fn = "test_cube_bad_dim.cube";
    std::ofstream ofs(fn);
    ofs << "comment line 1\n";
    ofs << "comment line 2\n";
    ofs << "1 0.0 0.0 0.0\n";
    ofs << "0 1.0 0.0 0.0\n";
    ofs << "2 0.0 1.0 0.0\n";
    ofs << "2 0.0 0.0 1.0\n";
    ofs << "14 4.0 0.0 0.0 0.0\n";
    ofs.close();

    EXPECT_FALSE(call_read_cube(fn));
    std::remove(fn.c_str());
}

TEST_F(ReadCubeInvalidTest, NegativeNatom)
{
    const std::string fn = "test_cube_neg_natom.cube";
    std::ofstream ofs(fn);
    ofs << "comment line 1\n";
    ofs << "comment line 2\n";
    ofs << "-1 0.0 0.0 0.0\n";
    ofs << "2 1.0 0.0 0.0\n";
    ofs << "2 0.0 1.0 0.0\n";
    ofs << "2 0.0 0.0 1.0\n";
    ofs.close();

    EXPECT_FALSE(call_read_cube(fn));
    std::remove(fn.c_str());
}

TEST_F(ReadCubeInvalidTest, ReadVdataPalgridFails)
{
    const std::string fn = "test_cube_truncated_palgrid.cube";
    std::ofstream ofs(fn);
    ofs << "comment line 1\n";
    ofs << "comment line 2\n";
    ofs << "1 0.0 0.0 0.0\n";
    ofs << "2 1.0 0.0 0.0\n";
    ofs << "2 0.0 1.0 0.0\n";
    ofs << "2 0.0 0.0 1.0\n";
    ofs << "14 4.0 0.0 0.0 0.0\n";
    ofs << "1.0 2.0 3.0\n";
    ofs.close();

    const int nx = 2;
    const int ny = 2;
    const int nz = 2;
    const int nrxx = nx * ny * nz;
    std::vector<double> data(nrxx, 0.0);
    Parallel_Grid pgrid(nx, ny, nz, nz, nrxx, nz, 1);
    std::ofstream ofs_running("unittest_read_cube.log");

    // read_vdata_palgrid now aborts the run on a malformed cube file so that
    // no rank is left waiting in the grid distribution below.
    EXPECT_DEATH(ModuleIO::read_vdata_palgrid(pgrid, 0, ofs_running, fn, data.data(), 1), "");
    std::remove(fn.c_str());
}
