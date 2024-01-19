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
    int typat[5] = {1,2,2,2,2};
    double xred[15] = {
        -0.279789547400000, 7.109405980000000E-002, 0.000000000000000E+000,
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
    double* rhor = new double[nfft];
    std::ifstream ifs_vks("vks.dat");
    std::ifstream ifs_vxc("vxc.dat");
    std::ifstream ifs_rhor("rhor.dat");

    for(int i = 0; i < nfft; i++)
    {
        ifs_vks >> vks[i];
        ifs_vxc >> vxc[i];
        ifs_rhor >> rhor[i];
    }

    double * force = new double[15];
    paw_cell.calculate_force(vks,vxc,rhor,force);

    double force_ref[15] = {
        0.00577967131142355, 0.00205616575556812,-4.88151255062632e-07,
        -0.363338769176905, 1.02333816217455, -3.50206921722388e-07,
        -0.36452161929387, -0.515139168331754, -0.888193793442404,
        -0.364521875200586, -0.51513940804039, 0.888192716951071,
        1.08375115425131, 0.001417226939304, -3.89259689050867e-07 };

    for(int i = 0; i < 15; i ++)
    {
        EXPECT_NEAR(force[i],force_ref[i],1e-10);
    }

    delete[] vks;
    delete[] vxc;
    delete[] force;
}