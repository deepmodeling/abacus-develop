#include "gtest/gtest.h"
#include <fstream>
#include <iostream>

#include "../paw_atom.h"
/*

Unit Test for subroutines related to the processing of rhoij,
the on-site density matrix to be used for PAW:
1. set_ca, which passes <psi|ptilde> from outside and saves it
2. accumulate_rhoij, which accumulates the contribution of one band
3. convert_rhoij, which converts rhoij to format required by PAW

*/

class Test_Paw_Atom : public testing::Test
{
    protected:

    Paw_Atom paw_atom;
    const int nproj = 8;
};

TEST_F(Test_Paw_Atom, test_paw)
{
    paw_atom.init_paw_atom(nproj);
    paw_atom.reset_rhoij();

    std::vector<std::complex<double>> ca;
    ca.resize(nproj);
    std::ifstream ifs_ca("ca.dat");
    std::ifstream ifs_rhoij("rhoij.dat");

    // there are altogether 4 bands
    const int nband = 4;
    for(int iband = 0; iband < nband; iband ++)
    {
        for(int ip = 0; ip < nproj; ip ++)
        {
            double re,im;
            ifs_ca >> re >> im;
            ca[ip] = std::complex<double>(re,im);
        }

        paw_atom.set_ca(ca, 2.0); //pass coefficient into the class
        paw_atom.accumulate_rhoij(); //accumulate the contribution of current band
    }

    std::vector<double> rhoij = paw_atom.get_rhoij();

    for(int i=0; i<rhoij.size();i++)
    {
        double tmp;
        ifs_rhoij >> tmp;
        EXPECT_NEAR(tmp,rhoij[i],1e-8);        
    }

    paw_atom.convert_rhoij();

    int nrhoijsel = paw_atom.get_nrhoijsel();
    EXPECT_EQ(nrhoijsel,21);

    std::vector<int> rhoijselect = paw_atom.get_rhoijselect();
    std::vector<int> rhoijselect_ref = {1,2,3,7,8,10,11,12,14,15,22,23,25,26,
        28,29,30,32,33,35,36,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    for(int i = 0; i < rhoijselect.size(); i ++)
    {
        EXPECT_EQ(rhoijselect[i],rhoijselect_ref[i]-1); // -1 because Fortran index starts from 1
    }

    std::vector<double> rhoijp = paw_atom.get_rhoijp();
    for(int i=0; i<rhoijp.size();i++)
    {
        double tmp;
        ifs_rhoij >> tmp;
        EXPECT_NEAR(tmp,rhoijp[i],1e-8);        
    }
}