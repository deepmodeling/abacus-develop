#include "gtest/gtest.h"
#include <fstream>

#include "../paw_element.h"

class Test_Read_Paw : public testing::Test
{
    protected:

    Paw_Element paw_element;
};

TEST_F(Test_Read_Paw, test_paw)
{
    double ecutpaw = 50;
    double cellfac = 1.2;
    
    paw_element.init_paw(ecutpaw, cellfac);
    paw_element.read_paw_xml("H.LDA_PW-JTH.xml");
    // I have not done any checks for now; will add later
}

class Test_SphBes_Func : public testing::Test
{
    protected:

};

// Test the first 6 spherical bessel functions
TEST_F(Test_SphBes_Func, test_paw)
{
    double dr = 0.01;
    int    nr = 1000;
    double sph_bes_ref;

    std::ifstream ifs("sphbes_ref.dat");

    for(int l = 0; l < 6; l ++)
    {
        for(int ir = 0; ir < nr; ir ++)
        {
            double r = double(ir) * dr;
            double sph_bes = Paw_Element::spherical_bessel_function(l,r);
            ifs >> sph_bes_ref;
            EXPECT_NEAR(sph_bes,sph_bes_ref,1e-8);
        }
    }
}

class Test_SphBes_Transform : public testing::Test
{
    protected:

    Paw_Element paw_element;
};

// This is the first projector
TEST_F(Test_SphBes_Transform, test_paw)
{
    paw_element.read_paw_xml("H.LDA_PW-JTH.xml");
    
    std::ifstream ifs_fr("func.dat");
    std::ifstream ifs_q("qlist.dat");
    std::ifstream ifs_fq("fq_ref.dat");
    
    std::vector<double> fr;
    std::vector<double> qlist;

    fr.resize(1500);
    for(int ir=0; ir<1500; ir++)
    {
        ifs_fr >> fr[ir];
    }

    qlist.resize(2999);
    for(int iq=0; iq<2999; iq++)
    {
        double q;
        ifs_q >> qlist[iq];
    }

    for(int iq=0; iq<2999; iq++)
    {
        double fq = paw_element.spherical_bessel_transform(0, fr, qlist[iq]);
        double fq_ref;
        ifs_fq >> fq_ref;

        EXPECT_NEAR(fq,fq_ref,1e-8);
    }

}