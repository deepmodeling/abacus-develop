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