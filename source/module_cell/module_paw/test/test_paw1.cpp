#include "gtest/gtest.h"
#include <fstream>

#include "../paw_cell.h"

class Test_Paw_Cell : public testing::Test
{
    protected:

    Paw_Cell paw_cell;

    double ecut = 50.0, cell_factor = 1.2, omega = 1.0;
    int nat = 5, ntyp = 3;
    int atom_type[5] = {0,1,2,1,2}; // Fe,O,H,O,H
    std::vector<std::string> filename_list;
    int nx = 1, ny = 1, nz = 1;
    double ** atom_coord;
    std::complex<double> *eigts1_in, *eigts2_in, *eigts3_in;
};

TEST_F(Test_Paw_Cell, test_paw)
{

    filename_list.resize(3);
    filename_list[0] = "Fe.GGA_PBE-JTH.xml";
    filename_list[1] =  "O.GGA_PBE-JTH.xml";
    filename_list[2] =  "H.LDA_PW-JTH.xml";

    /*
    Fe : mstate = 1+1+3+3+5+5 = 18
    <valence_states>
    <state n=" 3" l="0" f=" 2.0000000E+00" rc=" 2.0126126692" e="-3.2883573E+00" id= "Fe1"/>
    <state n=" 4" l="0" f=" 1.0000000E+00" rc=" 2.0126126692" e="-1.5638272E-01" id= "Fe2"/>
    <state n=" 3" l="1" f=" 6.0000000E+00" rc=" 1.8068708301" e="-2.0449140E+00" id= "Fe3"/>
    <state        l="1"                    rc=" 1.8068708301" e=" 1.7500000E+00" id= "Fe4"/>
    <state n=" 3" l="2" f=" 7.0000000E+00" rc=" 2.0126126692" e="-1.4018767E-01" id= "Fe5"/>
    <state        l="2"                    rc=" 2.0126126692" e=" 1.0000000E+00" id= "Fe6"/>
    </valence_states>    
    O : mstate = 1+1+3+3 = 8
    <valence_states>
    <state n=" 2" l="0" f=" 2.0000000E+00" rc=" 1.4146523028" e="-8.8057208E-01" id=  "O1"/>
    <state        l="0"                    rc=" 1.4146523028" e=" 1.0000000E+00" id=  "O2"/>
    <state n=" 2" l="1" f=" 4.0000000E+00" rc=" 1.4146523028" e="-3.3186932E-01" id=  "O3"/>
    <state        l="1"                    rc=" 1.4146523028" e=" 1.0000000E+00" id=  "O4"/>
    </valence_states>    
    H : mstate = 1+1+3 = 5
    <valence_states>
    <state n=" 1" l="0" f=" 1.0000000E+00" rc=" 0.9949503343" e="-2.3345876E-01" id=  "H1"/>
    <state        l="0"                    rc=" 0.9949503343" e=" 6.0000000E+00" id=  "H2"/>
    <state        l="1"                    rc=" 0.9949503343" e=" 1.2500000E+00" id=  "H3"/>
    </valence_states>    
    */

    atom_coord = new double * [5];
    for(int ia = 0; ia < nat; ia ++)
    {
        atom_coord[ia] = new double [3];
    }

    eigts1_in = new std::complex<double> [nx];
    eigts2_in = new std::complex<double> [ny];
    eigts3_in = new std::complex<double> [nz];

    paw_cell.init_paw_cell(ecut, cell_factor, omega, nat, ntyp, 
        atom_type, (const double **) atom_coord, filename_list,
        nx, ny, nz, eigts1_in, eigts2_in, eigts3_in);

    int nproj_tot = paw_cell.get_nproj_tot();
    EXPECT_EQ(nproj_tot,44);// 18 + 2 * 8 + 2 * 5 = 44

    int lmax = paw_cell.get_lmax();
    EXPECT_EQ(lmax,2);

    std::vector<int> iprj_to_ia = paw_cell.get_iprj_to_ia();
    EXPECT_EQ(iprj_to_ia.size(),44);
    for(int ip = 0; ip < 18; ip ++) EXPECT_EQ(iprj_to_ia[ip],0);
    for(int ip = 0; ip <  8; ip ++) EXPECT_EQ(iprj_to_ia[18+ip],1);
    for(int ip = 0; ip <  5; ip ++) EXPECT_EQ(iprj_to_ia[26+ip],2);
    for(int ip = 0; ip <  8; ip ++) EXPECT_EQ(iprj_to_ia[31+ip],3);
    for(int ip = 0; ip <  5; ip ++) EXPECT_EQ(iprj_to_ia[39+ip],4);

    std::vector<int> iprj_to_im = paw_cell.get_iprj_to_im();
    EXPECT_EQ(iprj_to_im.size(),44);
    for(int ip = 0; ip < 18; ip ++) EXPECT_EQ(iprj_to_im[ip],ip);
    for(int ip = 0; ip <  8; ip ++) EXPECT_EQ(iprj_to_im[18+ip],ip);
    for(int ip = 0; ip <  5; ip ++) EXPECT_EQ(iprj_to_im[26+ip],ip);
    for(int ip = 0; ip <  8; ip ++) EXPECT_EQ(iprj_to_im[31+ip],ip);
    for(int ip = 0; ip <  5; ip ++) EXPECT_EQ(iprj_to_im[39+ip],ip);

    std::vector<int> iprj_to_il = paw_cell.get_iprj_to_il();
    EXPECT_EQ(iprj_to_il.size(),44);
    std::vector<int> iprj_to_il_ref = {0,1,2,2,2,3,3,3,4,4,4,4,4,5,5,5,5,5,
        0,1,2,2,2,3,3,3,0,1,2,2,2,0,1,2,2,2,3,3,3,0,1,2,2,2};
    for(int ip = 0; ip < 44; ip ++)
    {
        EXPECT_EQ(iprj_to_il[ip],iprj_to_il_ref[ip]);
    }

    std::vector<int> iprj_to_l = paw_cell.get_iprj_to_l();
    EXPECT_EQ(iprj_to_l.size(),44);
    std::vector<int> iprj_to_l_ref = {0,0,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,
        0,0,1,1,1,1,1,1,0,0,1,1,1,0,0,1,1,1,1,1,1,0,0,1,1,1};
    for(int ip = 0; ip < 44; ip ++)
    {
        EXPECT_EQ(iprj_to_l[ip],iprj_to_l_ref[ip]);
    }

    std::vector<int> iprj_to_m = paw_cell.get_iprj_to_m();
    EXPECT_EQ(iprj_to_m.size(),44);
    std::vector<int> iprj_to_m_ref = {0,0,-1,0,1,-1,0,1,-2,-1,0,1,2,-2,-1,0,1,2,
        0,0,-1,0,1,-1,0,1,0,0,-1,0,1,0,0,-1,0,1,-1,0,1,0,0,-1,0,1};
    for(int ip = 0; ip < 44; ip ++)
    {
        EXPECT_EQ(iprj_to_m[ip],iprj_to_m_ref[ip]);
    }

    std::vector<int> start_iprj = paw_cell.get_start_iprj();
    EXPECT_EQ(start_iprj.size(),5);
    std::vector<int> start_iprj_ref = {0, 18, 26, 31, 39};
    for(int ia = 0; ia < 5; ia ++)
    {
        EXPECT_EQ(start_iprj[ia],start_iprj_ref[ia]);
    }

    for(int ia = 0; ia < nat; ia ++)
    {
        delete[] atom_coord[ia];
    }
    delete[] atom_coord;

    delete[] eigts1_in;
    delete[] eigts2_in;
    delete[] eigts3_in;
}

class Test_Ylm : public testing::Test
{
    protected:
};

TEST_F(Test_Ylm, test_paw)
{

    std::ifstream ifs_kpg("kpg.dat");
    std::ifstream ifs_ylm("ylm_ref.dat");

    int nkpg = 1491;

    for(int ikpg = 0; ikpg < nkpg; ikpg ++)
    {
        double r[3];
        ifs_kpg >> r[0] >> r[1] >> r[2];

        int lmax = 5;

        std::vector<double> ylm = Paw_Cell::calc_ylm(lmax, r);
        EXPECT_EQ(ylm.size(), 36);

        for(int i = 0; i < 36; i ++)
        {
            double ylm_ref;
            ifs_ylm >> ylm_ref;
            EXPECT_NEAR(ylm_ref,ylm[i],1e-8);
        }
    }
}

class Test_Leg_Pol : public testing::Test
{
    protected:
};

TEST_F(Test_Leg_Pol, test_paw)
{

    std::ifstream ifs("leg_pol.dat");

    int count = 29820;
    for(int i = 0; i < count; i ++)
    {
        int l,m;
        double arg,ref;
        ifs >> l >> m >> arg >> ref;
        double result = Paw_Cell::ass_leg_pol(l,m,arg);
        EXPECT_NEAR(ref,result,1e-8);
    }
}

class Test_PAW_Cell_k : public testing::Test
{
    protected:

    Paw_Cell paw_cell;

    double ecut = 20.0, cell_factor = 1.2, omega = 265.302, tpiba = 0.6159985595;
    int nat = 2, ntyp = 1;
    int atom_type[2] = {0,0}; // Si, Si
    std::vector<std::string> filename_list;
    int nx = 24, ny = 24, nz = 24;
    double ** atom_coord;
    std::complex<double> *eigts1_in, *eigts2_in, *eigts3_in;
};

TEST_F(Test_PAW_Cell_k, test_paw)
{
    // The subroutines are executed but the results are not checked yet;
    // it is because the results are only intermediate and are private members of the class
    // I will check the final results later once I finish them

    atom_coord = new double * [2];
    atom_coord[0] = new double [3];
    atom_coord[0][0] = 0.0; atom_coord[0][1] = 0.0; atom_coord[0][2] = 0.0;
    atom_coord[1] = new double [3];
    atom_coord[1][0] = 0.25; atom_coord[1][1] = 0.25; atom_coord[1][2] = 0.25;

    eigts1_in = new std::complex<double> [nat * (2 * nx + 1)];
    eigts2_in = new std::complex<double> [nat * (2 * ny + 1)];
    eigts3_in = new std::complex<double> [nat * (2 * nz + 1)];

    filename_list.resize(1);
    filename_list[0] = "Si_test.xml";

    std::ifstream ifs_eigts("eigts.dat");

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
        atom_type, (const double **) atom_coord, filename_list, nx, ny, nz,
        eigts1_in, eigts2_in, eigts3_in);

    for(int ia = 0; ia < nat; ia ++)
    {
        delete[] atom_coord[ia];
    }
    delete[] atom_coord;

    delete[] eigts1_in;
    delete[] eigts2_in;
    delete[] eigts3_in;

    //=========================================

    int npw = 411;
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

    EXPECT_EQ(vkb.size(),16);
    EXPECT_EQ(vkb[0].size(),npw);

    std::ifstream ifs_vkb("vkb_ref.dat");
    for(int iproj = 0; iproj < 16; iproj ++)
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

    for(int iband = 0; iband < nband; iband ++)
    {
        for(int ipw = 0; ipw < npw; ipw ++)
        {
            ifs_psi >> psi[ipw];
        }
        paw_cell.accumulate_rhoij(psi,weight[iband]);
    }

    delete[] psi;

    std::ifstream ifs_rhoij("rhoij1.dat");

    std::vector<std::vector<double>> rhoij = paw_cell.get_rhoij();
    EXPECT_EQ(rhoij.size(),nat);
    for(int iat = 0; iat < nat; iat ++)
    {
        EXPECT_EQ(rhoij[iat].size(),36);
        for(int i = 0; i < 36; i ++)
        {
            double tmp;
            ifs_rhoij >> tmp;

            EXPECT_NEAR(tmp,rhoij[iat][i],1e-8);
        }
    }

    std::vector<std::vector<double>> rhoijp;
    std::vector<std::vector<int>> rhoijselect;
    std::vector<int> nrhoijsel;

    paw_cell.get_rhoijp(rhoijp, rhoijselect, nrhoijsel);

    EXPECT_EQ(rhoijp.size(),nat);
    EXPECT_EQ(rhoijselect.size(),nat);
    EXPECT_EQ(nrhoijsel.size(),nat);

    for(int iat = 0; iat < nat; iat ++)
    {
        // As all entries are larger than 1e-10, rhoijp is the same as rhoij
        EXPECT_EQ(rhoijp[iat].size(),36);
        EXPECT_EQ(rhoijselect[iat].size(),36);
        EXPECT_EQ(nrhoijsel[iat],36);
        for(int i = 0; i < 36; i ++)
        {
            EXPECT_EQ(rhoijselect[iat][i],i);
            EXPECT_NEAR(rhoijp[iat][i],rhoij[iat][i],1e-8);
        }
    }
}