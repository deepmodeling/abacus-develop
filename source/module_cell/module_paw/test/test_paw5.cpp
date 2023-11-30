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

    paw_cell.set_paw_k(npw, kpt, ig_to_ix, ig_to_iy, ig_to_iz, (const double **) kpg, tpiba);

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
    psi = new std::complex<double>[npw];
    const int nband = 6;
    std::vector<double> weight={2,2,2,2,0,0};

    std::ifstream ifs_psi("psi.dat");

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

/*
    const int nproj = 8;
    double** dij;
    dij = new double*[1];
    dij[0]=new double[nproj * nproj];
    for(int i = 0; i < nproj * nproj; i ++)
    {
        dij[0][i] = 0.0;
    }
    dij[0][0] = 13.407893;
    dij[0][9] = 0.8201733412;
    dij[0][18] = 5.491609854;
    dij[0][27] = 5.491609854;
    dij[0][36] = 5.491609854;
    dij[0][45] = 0.59649632;
    dij[0][54] = 0.59649632;
    dij[0][63] = 0.59649632;

    std::vector<double> sij;
    sij.resize(nproj * nproj);
    for(int i = 0; i < sij.size(); i ++)
    {
        sij[i] = 0.0;
    }
    sij[0] = -4.902127221589223E-002;
    sij[9] = -9.18672607663861;
    sij[18] = -6.319002149104143E-003;
    sij[27] = -6.319002149104143E-003;
    sij[36] = -6.319002149104143E-003;
    sij[45] = -2.38515151080165;
    sij[54] = -2.38515151080165;
    sij[63] = -2.38515151080165;
    
    sij[1] = 0.726604973599628;
    sij[8] = 0.726604973599628;
    sij[21] = 0.156822922280989;
    sij[30] = 0.156822922280989;
    sij[39] = 0.156822922280989;
    sij[42] = 0.156822922280989;
    sij[51] = 0.156822922280989;
    sij[60] = 0.156822922280989;

    for(int iat = 0; iat < nat; iat ++)
    {
        paw_cell.set_dij(iat,dij);
        paw_cell.set_sij(iat,sij.data());
    }
    delete[] dij[0];
    delete[] dij;

    psi = new std::complex<double>[npw];
    std::complex<double> *vnlpsi, *snlpsi;
    vnlpsi = new std::complex<double>[npw];
    snlpsi = new std::complex<double>[npw];

    ifs_psi.clear();
    ifs_psi.seekg (0, std::ios::beg);

    std::ifstream ifs_vnlpsi("vnlpsi_ref.dat");
    std::ifstream ifs_snlpsi("snlpsi_ref.dat");

    std::cout << std::setprecision(10);
    for(int iband = 0; iband < nband; iband ++)
    {
        for(int ipw = 0; ipw < npw; ipw ++)
        {
            ifs_psi >> psi[ipw];
            vnlpsi[ipw] = 0.0;
        }

        paw_cell.paw_nl_psi(0, psi, vnlpsi);
        paw_cell.paw_nl_psi(1, psi, snlpsi);

        for(int ipw = 0; ipw < npw; ipw ++)
        {
            std::complex<double> tmp;
            ifs_vnlpsi >> tmp;
            EXPECT_NEAR(tmp.real(),vnlpsi[ipw].real(),1e-8);
            EXPECT_NEAR(tmp.imag(),vnlpsi[ipw].imag(),1e-8);

            ifs_snlpsi >> tmp;
            EXPECT_NEAR(tmp.real(),snlpsi[ipw].real(),1e-8);
            EXPECT_NEAR(tmp.imag(),snlpsi[ipw].imag(),1e-8);
        }
    }

    delete[] psi;
    delete[] vnlpsi;
    delete[] snlpsi;
*/
}