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

    double ecut = 10.0;
    double ecutpaw = 40.0;
    paw_cell.set_libpaw_ecut(ecut,ecutpaw);

    latvec.e11 = 1.0; latvec.e12 = 0.0; latvec.e13 = 0.0;
    latvec.e21 = 0.0; latvec.e22 = 1.0; latvec.e23 = 0.0;
    latvec.e31 = 0.0; latvec.e32 = 0.0; latvec.e33 = 1.0;

    double lat0 = 10.0;

    paw_cell.set_libpaw_cell(latvec, lat0);

    int nx = 30;
    int ny = 30;
    int nz = 30;
    int nx_dg = 60;
    int ny_dg = 60;
    int nz_dg = 60;

    paw_cell.set_libpaw_fft(nx, ny, nz, nx_dg, ny_dg, nz_dg);

    int natom = 5;
    int ntypat = 2;
    int typat[5] = {2,1,1,1,1};
    double xred[15] = {-0.279789547400000, 7.109405980000000E-002, 0.000000000000000E+000,
        -0.212391634800000,      -0.119543389500000,       0.000000000000000E+000,
        -0.212388153900000,       0.166411496700000,       0.165096194500000,
        -0.212388153900000,       0.166411496700000,      -0.165096194500000,
        -0.481990228200000,       7.109655050000001E-002,  0.000000000000000E+000,};

    paw_cell.set_libpaw_atom(natom, ntypat, typat, xred);

    paw_cell.set_libpaw_files();

    paw_cell.set_libpaw_xc(1,7);

    paw_cell.set_nspin(1);

    EXPECT_NEAR(paw_cell.get_libpaw_ecut(),10.0,1e-10);
    EXPECT_NEAR(paw_cell.get_libpaw_ecutpaw(),40.0,1e-10);
    EXPECT_NEAR(paw_cell.get_libpaw_ucvol(),1000.0,1e-10);
    
    EXPECT_NEAR(paw_cell.get_libpaw_rprimd()[0],10.0,1e-10);
    EXPECT_NEAR(paw_cell.get_libpaw_rprimd()[1], 0.0,1e-10);
    EXPECT_NEAR(paw_cell.get_libpaw_rprimd()[2], 0.0,1e-10);
    EXPECT_NEAR(paw_cell.get_libpaw_rprimd()[3], 0.0,1e-10);
    EXPECT_NEAR(paw_cell.get_libpaw_rprimd()[4],10.0,1e-10);
    EXPECT_NEAR(paw_cell.get_libpaw_rprimd()[5], 0.0,1e-10);
    EXPECT_NEAR(paw_cell.get_libpaw_rprimd()[6], 0.0,1e-10);
    EXPECT_NEAR(paw_cell.get_libpaw_rprimd()[7], 0.0,1e-10);
    EXPECT_NEAR(paw_cell.get_libpaw_rprimd()[8],10.0,1e-10);

    EXPECT_NEAR(paw_cell.get_libpaw_gprimd()[0],0.1,1e-10);
    EXPECT_NEAR(paw_cell.get_libpaw_gprimd()[1],0.0,1e-10);
    EXPECT_NEAR(paw_cell.get_libpaw_gprimd()[2],0.0,1e-10);
    EXPECT_NEAR(paw_cell.get_libpaw_gprimd()[3],0.0,1e-10);
    EXPECT_NEAR(paw_cell.get_libpaw_gprimd()[4],0.1,1e-10);
    EXPECT_NEAR(paw_cell.get_libpaw_gprimd()[5],0.0,1e-10);
    EXPECT_NEAR(paw_cell.get_libpaw_gprimd()[6],0.0,1e-10);
    EXPECT_NEAR(paw_cell.get_libpaw_gprimd()[7],0.0,1e-10);
    EXPECT_NEAR(paw_cell.get_libpaw_gprimd()[8],0.1,1e-10);

    EXPECT_NEAR(paw_cell.get_libpaw_gmet()[0],0.01,1e-10);
    EXPECT_NEAR(paw_cell.get_libpaw_gmet()[1],0.0 ,1e-10);
    EXPECT_NEAR(paw_cell.get_libpaw_gmet()[2],0.0 ,1e-10);
    EXPECT_NEAR(paw_cell.get_libpaw_gmet()[3],0.0 ,1e-10);
    EXPECT_NEAR(paw_cell.get_libpaw_gmet()[4],0.01,1e-10);
    EXPECT_NEAR(paw_cell.get_libpaw_gmet()[5],0.0 ,1e-10);
    EXPECT_NEAR(paw_cell.get_libpaw_gmet()[6],0.0 ,1e-10);
    EXPECT_NEAR(paw_cell.get_libpaw_gmet()[7],0.0 ,1e-10);
    EXPECT_NEAR(paw_cell.get_libpaw_gmet()[8],0.01,1e-10);

    EXPECT_EQ(paw_cell.get_libpaw_ngfft()[0],30);
    EXPECT_EQ(paw_cell.get_libpaw_ngfft()[1],30);
    EXPECT_EQ(paw_cell.get_libpaw_ngfft()[2],30);

    EXPECT_EQ(paw_cell.get_libpaw_ngfftdg()[0],60);
    EXPECT_EQ(paw_cell.get_libpaw_ngfftdg()[1],60);
    EXPECT_EQ(paw_cell.get_libpaw_ngfftdg()[2],60);

    EXPECT_EQ(paw_cell.get_libpaw_natom(),5);
    EXPECT_EQ(paw_cell.get_libpaw_ntypat(),2);

    EXPECT_EQ(paw_cell.get_libpaw_typat()[0],2);
    EXPECT_EQ(paw_cell.get_libpaw_typat()[1],1);
    EXPECT_EQ(paw_cell.get_libpaw_typat()[2],1);
    EXPECT_EQ(paw_cell.get_libpaw_typat()[3],1);
    EXPECT_EQ(paw_cell.get_libpaw_typat()[4],1);

    EXPECT_EQ(paw_cell.get_libpaw_filename_list()[0],'H');
    EXPECT_EQ(paw_cell.get_libpaw_filename_list()[1],'.');
    EXPECT_EQ(paw_cell.get_libpaw_filename_list()[2],'x');
    EXPECT_EQ(paw_cell.get_libpaw_filename_list()[3],'m');
    EXPECT_EQ(paw_cell.get_libpaw_filename_list()[4],'l');
    EXPECT_EQ(paw_cell.get_libpaw_filename_list()[264],'C');
    EXPECT_EQ(paw_cell.get_libpaw_filename_list()[265],'.');
    EXPECT_EQ(paw_cell.get_libpaw_filename_list()[266],'x');
    EXPECT_EQ(paw_cell.get_libpaw_filename_list()[267],'m');
    EXPECT_EQ(paw_cell.get_libpaw_filename_list()[268],'l');

    EXPECT_EQ(paw_cell.get_libpaw_ixc(),7);
    EXPECT_EQ(paw_cell.get_libpaw_xclevel(),1);

    EXPECT_EQ(paw_cell.get_nspin(),1);

#ifdef USE_PAW
    paw_cell.prepare_paw();

    double *vloc, *ncoret;
    int nfft = nx_dg * ny_dg * nz_dg;
    vloc = new double[nfft];
    ncoret = new double[nfft];

    paw_cell.get_vloc_ncoret(vloc, ncoret);
    std::ifstream ifs("fort.26");

    for(int i = 0; i < nfft; i++)
    {
        double tmp;
        ifs >> tmp;
        EXPECT_NEAR(tmp,vloc[i],1e-10);
    }

    for(int i = 0; i < nfft; i++)
    {
        double tmp;
        ifs >> tmp;
        EXPECT_NEAR(tmp,ncoret[i],1e-10);
    }

    delete[] vloc;
    delete[] ncoret;

    std::ifstream ifs_rhoij("rhoij");
    int nrhoijsel, *rhoijselect;
    double *rhoijp;
    int size_rhoij[5] = {36,15,15,15,15};
    for(int iat = 0; iat < natom; iat ++)
    {
        rhoijselect = new int[size_rhoij[iat]];
        rhoijp = new double[size_rhoij[iat]];

        ifs_rhoij >> nrhoijsel;
        for(int i = 0; i < size_rhoij[iat]; i++)
        {
            ifs_rhoij >> rhoijselect[i];
        }
        for(int i = 0; i < size_rhoij[iat]; i++)
        {
            ifs_rhoij >> rhoijp[i];
        }
        paw_cell.set_rhoij(iat,nrhoijsel,size_rhoij[iat],rhoijselect,rhoijp);

        delete[] rhoijselect;
        delete[] rhoijp;
    }

    double *nhat, *nhatgr;
    nhat = new double[nfft];
    nhatgr = new double[nfft*3];
    paw_cell.get_nhat(nhat,nhatgr);

    std::ifstream ifs_nhat("fort.19");
    for(int i=0; i<nfft; i++)
    {
        double tmp;
        ifs_nhat >> tmp;
        EXPECT_NEAR(tmp,nhat[i],1e-10);
    }

    delete[] nhat;
    delete[] nhatgr;

    std::ifstream ifs_veff("veff");
    double *vks, *vxc;
    vks = new double[nfft];
    vxc = new double[nfft];
    for(int i=0; i<nfft; i++)
    {
        ifs_veff >> vks[i];
    }
    for(int i=0; i<nfft; i++)
    {
        ifs_veff >> vxc[i];
    }

    paw_cell.calculate_dij(vks,vxc);
    delete[] vks;
    delete[] vxc;

    std::ifstream ifs_dij("dij_ref");
    double *dij;
    for(int iat = 0; iat < natom; iat ++)
    {
        dij = new double[size_rhoij[iat]];

        paw_cell.get_dij(iat,size_rhoij[iat],dij);

        for(int i=0; i<size_rhoij[iat]; i++)
        {
            double tmp;
            ifs_dij >> tmp;
            EXPECT_NEAR(tmp,dij[i],1e-10);
        }

        delete[] dij;
    }
#endif
}