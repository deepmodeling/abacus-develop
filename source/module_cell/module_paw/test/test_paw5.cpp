#include "gtest/gtest.h"
#include <fstream>
#include <iostream>

#include "../paw_cell.h"

class Test_PAW_Cell_k : public testing::Test
{
    protected:

    Paw_Cell paw_cell;
};

TEST_F(Test_PAW_Cell_k, test_paw)
{

    GlobalV::CAL_FORCE = 1;

    int natom = 5;
    int ntypat = 2;

    double ecut = 20.0;
    double cell_factor = 1.2;
    double omega = 1000;
    int nat = 5, ntyp = 2;
    int atom_type[5] = {1,0,0,0,0}; // C,H,H,H,H

    std::vector<std::string> filename_list;
    filename_list.resize(2);
    filename_list[0] = "H.xml";
    filename_list[1] = "C.xml";

    double ** atom_coord;
    atom_coord = new double * [nat];
    for(int ia = 0; ia < nat; ia ++)
    {
        atom_coord[ia] = new double[3];
    }

    atom_coord[0][0] = 0.7202104526;
    atom_coord[0][1] = 0.0710940598;
    atom_coord[0][2] = 0;
    atom_coord[1][0] = 0.7876083652;
    atom_coord[1][1] = 0.8804566105;
    atom_coord[1][2] = 0;
    atom_coord[2][0] = 0.7876118461;
    atom_coord[2][1] = 0.1664114967;
    atom_coord[2][2] = 0.1650961945;
    atom_coord[3][0] = 0.7876118461;
    atom_coord[3][1] = 0.1664114967;
    atom_coord[3][2] = 0.8349038055;
    atom_coord[4][0] = 0.5180097718;
    atom_coord[4][1] = 0.0710965505;
    atom_coord[4][2] = 0;

    paw_cell.init_paw_cell(ecut, cell_factor, omega, nat, ntyp,
        atom_type, (const double **) atom_coord, filename_list);

    for(int ia = 0; ia < nat; ia ++)
    {
        delete[] atom_coord[ia];
    }
    delete[] atom_coord;

    int nx = 30;
    int ny = 30;
    int nz = 30;

    std::ifstream ifs_eigts("eigts.dat");
    std::complex<double> *eigts1_in, *eigts2_in, *eigts3_in;
    eigts1_in = new std::complex<double> [nat * (2 * nx + 1)];
    eigts2_in = new std::complex<double> [nat * (2 * ny + 1)];
    eigts3_in = new std::complex<double> [nat * (2 * nz + 1)];

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

    paw_cell.set_eigts(nx,ny,nz,eigts1_in, eigts2_in, eigts3_in);

    delete[] eigts1_in;
    delete[] eigts2_in;
    delete[] eigts3_in;

    //=========================================

    int npw = 1503;
    double tpiba = 6.28318530718;
    int *ig_to_ix, *ig_to_iy, *ig_to_iz;
    double ** kpg = new double * [npw];    
    double kpt[3] = {0.0,0.0,0.0};

    ig_to_ix = new int[npw];
    ig_to_iy = new int[npw];
    ig_to_iz = new int[npw];

    std::ifstream ifs_igxyz("igxyz.dat");
    std::ifstream ifs_kpg("kpg1.dat");

    int ig;
    for(int i = 0; i < npw; i++)
    {
        ifs_igxyz >> ig >> ig_to_ix[i] >> ig_to_iy[i] >> ig_to_iz[i];
        kpg[i] = new double[3];
        ifs_kpg >> ig >> kpg[i][0] >> kpg[i][1] >> kpg[i][2];
    }

    int* isk_in = new int[1];
    isk_in[0] = 0;
    paw_cell.set_isk(1,isk_in);
    paw_cell.set_currentk(0);
    delete[] isk_in;

    // this is gamma point so gcar = kpg
    paw_cell.set_paw_k(npw, kpt, ig_to_ix, ig_to_iy, ig_to_iz, (const double **) kpg, tpiba, (const double**) kpg);

    delete[] ig_to_ix;
    delete[] ig_to_iy;
    delete[] ig_to_iz;
    for(int i = 0; i < npw; i++)
    {
        delete[] kpg[i];
    }
    delete[] kpg;

    paw_cell.get_vkb();
    auto vkb = paw_cell.output_vkb();

    EXPECT_EQ(vkb.size(),28);
    EXPECT_EQ(vkb[0].size(),npw);

    std::ifstream ifs_vkb("vkb_ref.dat");
    for(int iproj = 0; iproj < 28; iproj ++)
    {
        for(int ipw = 0; ipw < npw; ipw ++)
        {
            std::complex<double> tmp;
            ifs_vkb >> tmp;
            EXPECT_NEAR(tmp.real(),vkb[iproj][ipw].real(),1e-8);
            EXPECT_NEAR(tmp.imag(),vkb[iproj][ipw].imag(),1e-8);
        }
    }

    std::complex<double> *psi;
    const int nband = 6;
    psi = new std::complex<double>[npw*nband];
    std::vector<double> weight={2,2,2,2,0,0};

    std::ifstream ifs_psi("psi.dat");

    paw_cell.reset_rhoij();
    for(int iband = 0; iband < nband; iband ++)
    {
        for(int ipw = 0; ipw < npw; ipw ++)
        {
            ifs_psi >> psi[iband*npw+ipw];
        }
        paw_cell.accumulate_rhoij(&psi[iband*npw],weight[iband]);
    }

    std::vector<std::vector<double>> rhoijp;
    std::vector<std::vector<int>> rhoijselect;
    std::vector<int> nrhoijsel;

    paw_cell.get_rhoijp(rhoijp, rhoijselect, nrhoijsel);

    EXPECT_EQ(rhoijp.size(),nat);
    EXPECT_EQ(rhoijselect.size(),nat);
    EXPECT_EQ(nrhoijsel.size(),nat);

    std::ifstream ifs_rhoij("rhoij1.dat");
    for(int iat = 0; iat < nat; iat ++)
    {
        // As all entries are larger than 1e-10, rhoijp is the same as rhoij
        int n;
        ifs_rhoij >> n;
        EXPECT_EQ(nrhoijsel[iat],n);
        for(int i = 0; i < n; i ++)
        {
            int tmp;
            ifs_rhoij >> tmp;
            EXPECT_EQ(rhoijselect[iat][i],tmp);
        }
        for(int i = 0; i < n; i ++)
        {
            double tmp;
            ifs_rhoij >> tmp;
            EXPECT_NEAR(rhoijp[iat][i],tmp,1e-8);
        }
    }

    //set sij
    std::vector<double> sij;

    // C atom
    sij.resize(64);
    for(int i = 0; i < sij.size(); i ++)
    {
        sij[i] = 0.0;
    }
    sij[0] = -0.0589997365289;
    sij[1] = 0.747236491416;
    sij[8] = 0.747236491416;
    sij[9] = -9.30610950738;
    sij[18] = 0.0699425732126;
    sij[21] = -0.483805550235;
    sij[27] = 0.0699425732126;
    sij[30] = -0.483805550235;
    sij[36] = 0.0699425732126;
    sij[39] = -0.483805550235;
    sij[42] = -0.483805550235;
    sij[45] = 3.20343380497;
    sij[51] = -0.483805550235;
    sij[54] = 3.20343380497;
    sij[60] = -0.483805550235;
    sij[63] = 3.20343380497;

    paw_cell.set_sij(0,sij.data());

    // H atom
    sij.resize(25);
    for(int i = 0; i < sij.size(); i ++)
    {
        sij[i] = 0.0;
    }

    sij[0] = 0.00772422455533;
    sij[1] = 0.0123172347802;
    sij[5] = 0.0123172347802;
    sij[6] = 0.0196386056473;
    sij[12] = 0.000908274554872;
    sij[18] = 0.000908274554872;
    sij[24] = 0.000908274554872;

    for(int iat = 1; iat < nat; iat ++)
    {
        paw_cell.set_sij(iat,sij.data());
    }

    // set dij
    std::ifstream ifs_dij("dij_in.dat");
    double** dij;
    dij = new double*[1];

    // C atom
    dij[0]=new double[64];
    for(int i = 0; i < 64; i ++)
    {
        ifs_dij >> dij[0][i];
    }
    paw_cell.set_dij(0,dij);
    delete[] dij[0];

    dij[0] = new double[25];
    for(int iat = 1; iat < nat; iat ++)
    {
        for(int i = 0; i < 25; i ++)
        {
            ifs_dij >> dij[0][i];
        }
        paw_cell.set_dij(iat,dij);
    }
    delete[] dij[0];
    delete[] dij;

    double force[15];
    double epsilon[6] = {-1.1820349744334246e+00,-6.1373585531593766e-01,-6.1245076531494447e-01,
        -6.1185767211855080e-01,-8.1899553911957745e-02,1.7727505277597955e-01};
    
    paw_cell.paw_nl_force(psi,epsilon,weight.data(),6,force);

    double force_ref[15] = {
        0.00585756099172754,0.00503158883480039,2.0116281335185e-08,
        -0.0154573057564919,0.0441062361839232,7.31842461588338e-11,
        -0.0154715252567464,-0.0219321805387927,-0.0380705126948971,
        -0.0154715261725948,-0.021932181037128,0.0380705141839775,
        0.0467646730276522,4.17425941637841e-05,6.61390458383683e-11
    };
    
    for(int i = 0; i < 15; i ++)
    {
        EXPECT_NEAR(force[i],force_ref[i],1e-8);
    }
    delete[] psi;
}