#include "gtest/gtest.h"
#include <fstream>
#include <iostream>

#include "../paw_cell.h"

/*
Unit Test for paw_force, used to calculate PAW on-site contribution
to force and stress
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
    double ecutpaw = 10.0;
    paw_cell.set_libpaw_ecut(ecut,ecutpaw);

    latvec.e11 = 1.0; latvec.e12 = 0.0; latvec.e13 = 0.0;
    latvec.e21 = 0.0; latvec.e22 = 1.0; latvec.e23 = 0.0;
    latvec.e31 = 0.0; latvec.e32 = 0.0; latvec.e33 = 1.0;

    double lat0 = 10.0;

    paw_cell.set_libpaw_cell(latvec, lat0);

    int nx = 30;
    int ny = 30;
    int nz = 30;

    GlobalV::NPROC = 1;
    std::vector<int> start_z = {0};
    std::vector<int> num_z = {nz};
    paw_cell.set_libpaw_fft(nx, ny, nz, nx, ny, nz, start_z.data(), num_z.data());

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

    paw_cell.prepare_paw();

    int nfft = nx * ny * nz;

    std::ifstream ifs_rhoij("rhoij1");
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

    double* vks = new double[nfft];
    double* vxc = new double[nfft];
    std::ifstream ifs_vks("vks.dat");
    std::ifstream ifs_vxc("vxc.dat");

    for(int i = 0; i < nfft; i++)
    {
        ifs_vks >> vks[i];
        ifs_vxc >> vxc[i];
    }

    double * force = new double[15];
    paw_cell.calculate_force(vks,vxc,force);

    double force_ref[15] = {0.00250618058688 ,-0.00704773416561 ,1.05723313446e-09 ,
        0.00234697314963 ,0.00363698488924 ,0.00615245758785 ,
        0.00234696689977 ,0.00363697886513 ,-0.00615241947266 ,
        -0.0076578601745 ,0.000104586221829 ,2.03032806302e-09 ,
        0.000144056652824 ,-8.1352863205e-05 ,-2.42047466188e-07 };

    for(int i = 0; i < 15; i ++)
    {
        EXPECT_NEAR(force[i],force_ref[i],1e-10);
    }

    delete[] vks;
    delete[] vxc;
    delete[] force;
}