#include "gtest/gtest.h"
#include <fstream>
#include <iostream>

#include "../paw_cell.h"

/*

Unit Test for the following subroutines, which are used to pass information
from main ABACUS program to LibPAW:

1. set_libpaw_ecut, which sets kinetic energy cutoff
2. set_libpaw_cell, which sets quantities related to cell parameters
3. set_libpaw_fft, which sets the real-space FFT grid
4. set_libpaw_atom, which sets information of atoms in the unit cell
5. set_libpaw_files, which sets the names of PAW xml files

*/

class Test_Libpaw_Cell : public testing::Test
{
    protected:

    Paw_Cell paw_cell;
};

TEST_F(Test_Libpaw_Cell, test_paw)
{
    ModuleBase::Matrix3 latvec;

    double ecut = 30.0;
    paw_cell.set_libpaw_ecut(ecut,ecut);

    latvec.e11 = 0.5; latvec.e12 = 0.5; latvec.e13 = 0.0;
    latvec.e21 = 0.0; latvec.e22 = 0.5; latvec.e23 = 0.5;
    latvec.e31 = 0.5; latvec.e32 = 0.0; latvec.e33 = 0.5;

    double lat0 = 10.2;

    paw_cell.set_libpaw_cell(latvec, lat0);

    int nx = 30;
    int ny = 30;
    int nz = 30;

    paw_cell.set_libpaw_fft(nx, ny, nz, nx, ny, nz);

    int natom = 2;
    int ntypat = 2;
    int typat[2] = {1,2};
    double xred[6] = {0.0, 0.0, 0.0, 0.5, 0.5, 0.5};

    paw_cell.set_libpaw_atom(natom, ntypat, typat, xred);

    paw_cell.set_libpaw_files();

    EXPECT_NEAR(paw_cell.get_libpaw_ecut(),30.0,1e-10);
    EXPECT_NEAR(paw_cell.get_libpaw_ecutpaw(),30.0,1e-10);
    EXPECT_NEAR(paw_cell.get_libpaw_ucvol(),265.302,1e-10);
    
    EXPECT_NEAR(paw_cell.get_libpaw_rprimd()[0],5.1,1e-10);
    EXPECT_NEAR(paw_cell.get_libpaw_rprimd()[1],5.1,1e-10);
    EXPECT_NEAR(paw_cell.get_libpaw_rprimd()[2],0.0,1e-10);
    EXPECT_NEAR(paw_cell.get_libpaw_rprimd()[3],0.0,1e-10);
    EXPECT_NEAR(paw_cell.get_libpaw_rprimd()[4],5.1,1e-10);
    EXPECT_NEAR(paw_cell.get_libpaw_rprimd()[5],5.1,1e-10);
    EXPECT_NEAR(paw_cell.get_libpaw_rprimd()[6],5.1,1e-10);
    EXPECT_NEAR(paw_cell.get_libpaw_rprimd()[7],0.0,1e-10);
    EXPECT_NEAR(paw_cell.get_libpaw_rprimd()[8],5.1,1e-10);

    EXPECT_NEAR(paw_cell.get_libpaw_gprimd()[0],0.09803921568,1e-10);
    EXPECT_NEAR(paw_cell.get_libpaw_gprimd()[1],0.09803921568,1e-10);
    EXPECT_NEAR(paw_cell.get_libpaw_gprimd()[2],-0.09803921568,1e-10);
    EXPECT_NEAR(paw_cell.get_libpaw_gprimd()[3],-0.09803921568,1e-10);
    EXPECT_NEAR(paw_cell.get_libpaw_gprimd()[4],0.09803921568,1e-10);
    EXPECT_NEAR(paw_cell.get_libpaw_gprimd()[5],0.09803921568,1e-10);
    EXPECT_NEAR(paw_cell.get_libpaw_gprimd()[6],0.09803921568,1e-10);
    EXPECT_NEAR(paw_cell.get_libpaw_gprimd()[7],-0.09803921568,1e-10);
    EXPECT_NEAR(paw_cell.get_libpaw_gprimd()[8],0.09803921568,1e-10);

    EXPECT_NEAR(paw_cell.get_libpaw_gmet()[0],0.02883506343,1e-10);
    EXPECT_NEAR(paw_cell.get_libpaw_gmet()[1],-0.00961168781,1e-10);
    EXPECT_NEAR(paw_cell.get_libpaw_gmet()[2],-0.00961168781,1e-10);
    EXPECT_NEAR(paw_cell.get_libpaw_gmet()[3],-0.00961168781,1e-10);
    EXPECT_NEAR(paw_cell.get_libpaw_gmet()[4],0.02883506343,1e-10);
    EXPECT_NEAR(paw_cell.get_libpaw_gmet()[5],-0.00961168781,1e-10);
    EXPECT_NEAR(paw_cell.get_libpaw_gmet()[6],-0.00961168781,1e-10);
    EXPECT_NEAR(paw_cell.get_libpaw_gmet()[7],-0.00961168781,1e-10);
    EXPECT_NEAR(paw_cell.get_libpaw_gmet()[8],0.02883506343,1e-10);

    EXPECT_EQ(paw_cell.get_libpaw_ngfft()[0],30);
    EXPECT_EQ(paw_cell.get_libpaw_ngfft()[1],30);
    EXPECT_EQ(paw_cell.get_libpaw_ngfft()[2],30);

    EXPECT_EQ(paw_cell.get_libpaw_ngfftdg()[0],30);
    EXPECT_EQ(paw_cell.get_libpaw_ngfftdg()[1],30);
    EXPECT_EQ(paw_cell.get_libpaw_ngfftdg()[2],30);

    EXPECT_EQ(paw_cell.get_libpaw_natom(),2);
    EXPECT_EQ(paw_cell.get_libpaw_ntypat(),2);

    EXPECT_EQ(paw_cell.get_libpaw_typat()[0],1);
    EXPECT_EQ(paw_cell.get_libpaw_typat()[1],2);

    EXPECT_EQ(paw_cell.get_libpaw_filename_list()[0],'C');
    EXPECT_EQ(paw_cell.get_libpaw_filename_list()[1],'.');
    EXPECT_EQ(paw_cell.get_libpaw_filename_list()[2],'x');
    EXPECT_EQ(paw_cell.get_libpaw_filename_list()[3],'m');
    EXPECT_EQ(paw_cell.get_libpaw_filename_list()[4],'l');
    EXPECT_EQ(paw_cell.get_libpaw_filename_list()[264],'H');
    EXPECT_EQ(paw_cell.get_libpaw_filename_list()[265],'.');
    EXPECT_EQ(paw_cell.get_libpaw_filename_list()[266],'x');
    EXPECT_EQ(paw_cell.get_libpaw_filename_list()[267],'m');
    EXPECT_EQ(paw_cell.get_libpaw_filename_list()[268],'l');
}