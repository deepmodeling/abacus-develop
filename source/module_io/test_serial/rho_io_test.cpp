#include "module_io/cube_io.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "module_base/global_variable.h"
#include "module_io/cube_io.h"
#include "prepare_unitcell.h"
#include "module_hamilt_pw/hamilt_pwdft/parallel_grid.h"

#ifdef __LCAO
InfoNonlocal::InfoNonlocal()
{
}
InfoNonlocal::~InfoNonlocal()
{
}
LCAO_Orbitals::LCAO_Orbitals()
{
}
LCAO_Orbitals::~LCAO_Orbitals()
{
}
#endif
Magnetism::Magnetism()
{
    this->tot_magnetization = 0.0;
    this->abs_magnetization = 0.0;
    this->start_magnetization = nullptr;
}
Magnetism::~Magnetism()
{
    delete[] this->start_magnetization;
}
Parallel_Grid::~Parallel_Grid() {}

#define private public
#include "module_parameter/parameter.h"
#undef private

/***************************************************************
 *  unit test of read_rho, write_rho and trilinear_interpolate
 ***************************************************************/

/**
 * - Tested Functions:
 *   - read_rho()
 *     - the function to read_rho from file
 *     - the serial version without MPI
 *   - trilinear_interpolate()
 *     - the trilinear interpolation method
 *     - the serial version without MPI
 */

class RhoIOTest : public ::testing::Test
{
  protected:
    int nspin = 1;
    int nrxx = 36 * 36 * 36;
    int prenspin = 1;
    double** rho;
    UnitCell* ucell;

    int my_rank = 0;
    std::ofstream ofs_running = std::ofstream("unittest.log");

    void SetUp()
    {
        rho = new double*[nspin];
        ucell = new UnitCell;
        for (int is = 0; is < nspin; ++is)
        {
            rho[is] = new double[nrxx];
        }
    }
    void TearDown()
    {
        for (int is = 0; is < nspin; ++is)
        {
            delete[] rho[is];
        }
        delete[] rho;
        delete ucell;
    }
};

TEST_F(RhoIOTest, Read)
{
    int is = 0;
    std::string fn = "./support/SPIN1_CHG.cube";
    int nx = 36;
    int ny = 36;
    int nz = 36;
    double ef;
    UcellTestPrepare utp = UcellTestLib["Si"];
    ucell = utp.SetUcellInfo();
    Parallel_Grid pgrid(nx, ny, nz, nz, nrxx, nz, 1);
    ModuleIO::read_vdata_palgrid(pgrid, my_rank, ofs_running, fn, rho[is], ucell->nat);
    EXPECT_DOUBLE_EQ(rho[0][0], 1.27020863940e-03);
    EXPECT_DOUBLE_EQ(rho[0][46655], 1.33581335706e-02);
}

TEST_F(RhoIOTest, Write)
{
    int nx = 36;
    int ny = 36;
    int nz = 36;
    UcellTestPrepare utp = UcellTestLib["Si"];
    ucell = utp.SetUcellInfo();
    ucell->lat0 = 10.2;
    ucell->latvec = { -0.5,0,0.5,0,0.5,0.5,-0.5,0.5,0 };
    ucell->atoms[0].tau[0] = ModuleBase::Vector3<double>(0.0, 0.0, 0.0);
    ucell->atoms[0].tau[1] = ModuleBase::Vector3<double>(-0.75, 0.75, 0.75);
    ucell->atoms[0].ncpp.zv = 4;
    Parallel_Grid pgrid(nx, ny, nz, nz, nrxx, nz, 1);
    ModuleIO::read_vdata_palgrid(pgrid, my_rank, ofs_running, "support/SPIN1_CHG.cube", rho[0], ucell->nat);
    ModuleIO::write_vdata_palgrid(pgrid, rho[0], 0, nspin, 0, "test_write_vdata_palgrid.cube", 0.461002, ucell, 11, 1);
    EXPECT_EQ(system("diff -q test_write_vdata_palgrid.cube support/SPIN1_CHG.cube"), 0);
}

TEST_F(RhoIOTest, TrilinearInterpolate)
{
    int nx = 36;
    int ny = 40;
    int nz = 44;
    int nx_read = 36;
    int ny_read = 36;
    int nz_read = 36;
    std::ifstream ifs("./support/SPIN1_CHG.cube");
    for (int i = 0; i < 8; ++i)
    {
        ifs.ignore(300, '\n');
    }
    std::vector<double> data_read(nx_read * ny_read * nz_read);
    for (int ix = 0; ix < nx_read; ix++)
    {
        for (int iy = 0; iy < ny_read; iy++)
        {
            for (int iz = 0; iz < nz_read; iz++)
            {
                ifs >> data_read[(ix * ny_read + iy) * nz_read + iz];
            }
        }
    }

    // The old implementation is inconsistent: ifdef MPI, [x][y][z]; else, [z][x][y].
    // Now we use [x][y][z] for both MPI and non-MPI, so here we need to chage the index order.
    auto permute_xyz2zxy = [&](const double* const xyz, double* const zxy) -> void
        {
            for (int ix = 0; ix < nx; ix++)
            {
                for (int iy = 0; iy < ny; iy++)
                {
                    for (int iz = 0; iz < nz; iz++)
                    {
                        zxy[(iz * nx + ix) * ny + iy] = xyz[(ix * ny + iy) * nz + iz];
                    }
                }
            }
        };
    const int nxyz = nx * ny * nz;
    std::vector<double> data_xyz(nxyz);
    std::vector<double> data(nxyz); // z > x > y
    ModuleIO::trilinear_interpolate(data_read.data(), nx_read, ny_read, nz_read, nx, ny, nz, data_xyz.data());
    permute_xyz2zxy(data_xyz.data(), data.data());
    EXPECT_DOUBLE_EQ(data[0], 0.0010824725010374092);
    EXPECT_DOUBLE_EQ(data[10], 0.058649850374240906);
    EXPECT_DOUBLE_EQ(data[100], 0.018931708073604996);
}

struct CubeIOTest : public ::testing::Test
{
    std::vector<std::string> comment;
    int natom = 0;
    std::vector<double> origin;
    std::vector<int> nvoxel;
    int nx_read = 0;
    int ny_read = 0;
    int nz_read = 0;
    std::vector<double> dx;
    std::vector<double> dy;
    std::vector<double> dz;
    std::vector<std::vector<double>> axis_vecs;
    std::vector<int> atom_type;
    std::vector<double> atom_charge;
    std::vector<std::vector<double>> atom_pos;
    std::vector<double> data_read;
    const std::string fn = "./support/SPIN1_CHG.cube";
};


TEST_F(CubeIOTest, ReadCube)
{
    ModuleIO::read_cube(fn, comment, natom, origin, nx_read, ny_read, nz_read, dx, dy, dz, atom_type, atom_charge, atom_pos, data_read);
    EXPECT_EQ(comment[0], "STEP: 0  Cubefile created from ABACUS. Inner loop is z, followed by y and x");
    EXPECT_EQ(comment[1], "1 (nspin) 0.461002 (fermi energy, in Ry)");
    EXPECT_EQ(natom, 2);
    for (auto& o : origin) { EXPECT_EQ(o, 0.0); }
    EXPECT_EQ(nx_read, 36);
    EXPECT_EQ(ny_read, 36);
    EXPECT_EQ(nz_read, 36);
    EXPECT_DOUBLE_EQ(dx[0], -0.141667);
    EXPECT_DOUBLE_EQ(dy[2], 0.141667);
    EXPECT_DOUBLE_EQ(dz[1], 0.141667);
    EXPECT_EQ(atom_type.size(), natom);
    EXPECT_EQ(atom_charge.size(), natom);
    EXPECT_EQ(atom_pos.size(), natom);
    for (auto& t : atom_type) { EXPECT_EQ(t, 14); }
    for (auto& c : atom_charge) { EXPECT_DOUBLE_EQ(c, 4.0); }
    EXPECT_DOUBLE_EQ(atom_pos[1][1], 7.65);
    const int nxyz = nx_read * ny_read * nz_read;
    EXPECT_EQ(data_read.size(), nxyz);
    EXPECT_EQ(data_read[1], 2.64004483879e-03);
    EXPECT_EQ(data_read[nxyz - 1], 1.33581335706e-02);
}

TEST_F(CubeIOTest, WriteCube)
{
    ModuleIO::read_cube(fn, comment, natom, origin, nx_read, ny_read, nz_read, dx, dy, dz, atom_type, atom_charge, atom_pos, data_read);
    ModuleIO::write_cube("test_write.cube", comment, natom, origin, nx_read, ny_read, nz_read, dx, dy, dz, atom_type, atom_charge, atom_pos, data_read, 11);
    EXPECT_EQ(system("diff -q test_write.cube ./support/SPIN1_CHG.cube"), 0);
}