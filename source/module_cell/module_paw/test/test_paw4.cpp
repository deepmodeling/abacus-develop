#include "gtest/gtest.h"
#include <fstream>
#include <iostream>

#include "../paw_cell.h"

/*

Unit Test for PAW subroutines; is a combination of test3 and test1

1. set_libpaw_ecut, which sets kinetic energy cutoff
2. set_libpaw_cell, which sets quantities related to cell parameters
3. set_libpaw_fft, which sets the real-space FFT grid
4. set_libpaw_atom, which sets information of atoms in the unit cell
5. set_libpaw_files, which sets the names of PAW xml files

6. init_paw_cell and set_paw_k, which collects necessary information
7.  accumulate_rhoij, which calculates sum_n <psi_n|p_i><p_j|psi_n>
*/

class Test_PAW : public testing::Test
{
    protected:

    Paw_Cell paw_cell;

    double ecut = 20.0, cell_factor = 1.2, omega = 1000, tpiba = 0.62831853071;
    int nat = 5, ntyp = 2;
    int atom_type[5] = {1,0,0,0,0}; // Si, Si
    std::vector<std::string> filename_list;
    int nx = 30, ny = 30, nz = 30;
    double ** atom_coord;
    std::complex<double> *eigts1_in, *eigts2_in, *eigts3_in;
};

TEST_F(Test_PAW, test_paw)
{

// I'm setting up lots of things, to be worse I need to do it twice ..
// one for ABACUS side, the other for libpaw side
// might be confusing, will clean up later

    atom_coord = new double * [nat];
    for(int ia = 0; ia < nat; ia ++)
    {
        atom_coord[ia] = new double[3];
    }
    atom_coord[0][0] = -0.2797895474; atom_coord[0][1] = 7.10940598E-002; atom_coord[0][2] = 0.0;
    atom_coord[1][0] = -0.2123916348; atom_coord[1][1] = -0.1195433895; atom_coord[1][2] = 0.0;
    atom_coord[2][0] = -0.2123881539; atom_coord[2][1] = 0.1664114967; atom_coord[2][2] = 0.1650961945;
    atom_coord[3][0] = -0.2123881539; atom_coord[3][1] = 0.1664114967; atom_coord[3][2] = -0.1650961945;
    atom_coord[4][0] = -0.4819902282; atom_coord[4][1] = 7.10965505E-002; atom_coord[4][2] = 0.0;

    eigts1_in = new std::complex<double> [nat * (2 * nx + 1)];
    eigts2_in = new std::complex<double> [nat * (2 * ny + 1)];
    eigts3_in = new std::complex<double> [nat * (2 * nz + 1)];

    filename_list.resize(ntyp);
    filename_list[0] = "H.xml";
    filename_list[1] = "C.xml";

    ModuleBase::Matrix3 latvec;
    latvec.e11 = 1.0; latvec.e12 = 0.0; latvec.e13 = 0.0;
    latvec.e21 = 0.0; latvec.e22 = 1.0; latvec.e23 = 0.0;
    latvec.e31 = 0.0; latvec.e32 = 0.0; latvec.e33 = 1.0;

    double lat0 = 10.0;

    std::ifstream ifs_eigts("test4_eigts.dat");

    for(int i = 0; i < nat*(2*nx+1); i ++)
    {
        ifs_eigts >> eigts1_in[i];
    }
    for(int i = 0; i < nat*(2*ny+1); i ++)
    {
        ifs_eigts >> eigts2_in[i];
    }
    for(int i = 0; i < nat*(2*nz+1); i ++)
    {
        ifs_eigts >> eigts3_in[i];
    }

    paw_cell.init_paw_cell(ecut, cell_factor, omega, nat, ntyp,
        atom_type, (const double **) atom_coord, filename_list);
    paw_cell.set_eigts(nx, ny, nz,eigts1_in, eigts2_in, eigts3_in);

    delete[] eigts1_in;
    delete[] eigts2_in;
    delete[] eigts3_in;

    //=========================================

    int npw = 1503;
    int *ig_to_ix, *ig_to_iy, *ig_to_iz;
    double ** kpg = new double * [npw];    
    double kpt[3] = {0.0,0.0,0.0};

    ig_to_ix = new int[npw];
    ig_to_iy = new int[npw];
    ig_to_iz = new int[npw];

    std::ifstream ifs_igxyz("test4_igxyz.dat");
    std::ifstream ifs_kpg("test4_kpg.dat");

    int ig;
    for(int i = 0; i < npw; i++)
    {
        ifs_igxyz >> ig >> ig_to_ix[i] >> ig_to_iy[i] >> ig_to_iz[i];
        kpg[i] = new double[3];
        ifs_kpg >> ig >> kpg[i][0] >> kpg[i][1] >> kpg[i][2];
    }

    paw_cell.set_paw_k(npw, kpt, ig_to_ix, ig_to_iy, ig_to_iz, (const double **) kpg, tpiba);

    delete[] ig_to_ix;
    delete[] ig_to_iy;
    delete[] ig_to_iz;
    for(int i = 0; i < npw; i++)
    {
        delete[] kpg[i];
    }
    delete[] kpg;

    double *xred;
    int    *typat;

    xred = new double[nat * 3];
    typat = new int[nat];

    for(int iat = 0; iat < nat; iat ++)
    {
        typat[iat] = atom_type[iat]+1;
        for(int i = 0; i < 3; i ++)
        {
            xred[iat*3+i] = atom_coord[iat][i];
        }
    }

    for(int ia = 0; ia < nat; ia ++)
    {
        delete[] atom_coord[ia];
    }
    delete[] atom_coord;

    paw_cell.set_libpaw_ecut(ecut,ecut);

    paw_cell.set_libpaw_cell(latvec, lat0);

    paw_cell.set_libpaw_fft(nx, ny, nz, nx, ny, nz);

    paw_cell.set_libpaw_atom(nat, ntyp, typat, xred);

    paw_cell.set_libpaw_files();

    paw_cell.set_libpaw_xc(1,7);
    
    paw_cell.set_nspin(1);

    paw_cell.prepare_paw();

    paw_cell.get_vkb();

    std::complex<double> *psi;
    psi = new std::complex<double>[npw];
    const int nband = 6;
    std::vector<double> weight={2,2,2,2,0,0};

    std::ifstream ifs_psi("test4_psi.dat");

    paw_cell.reset_rhoij();
    for(int iband = 0; iband < nband; iband ++)
    {
        for(int ipw = 0; ipw < npw; ipw ++)
        {
            ifs_psi >> psi[ipw];
        }
        paw_cell.accumulate_rhoij(psi,weight[iband]);
    }

    delete[] psi;

    std::vector<std::vector<double>> rhoijp;
    std::vector<std::vector<int>> rhoijselect;
    std::vector<int> nrhoijsel;

    paw_cell.get_rhoijp(rhoijp, rhoijselect, nrhoijsel);

    for(int iat = 0; iat < nat; iat ++)
    {
        paw_cell.set_rhoij(iat,nrhoijsel[iat],rhoijp[iat].size(),rhoijselect[iat].data(),rhoijp[iat].data());
    }

    int nfft = nx * ny * nz;
    double *nhat, *nhatgr;
    nhat = new double[nfft];
    nhatgr = new double[nfft*3];
    paw_cell.get_nhat(nhat,nhatgr);
}