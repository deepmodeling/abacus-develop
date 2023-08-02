#include "gtest/gtest.h"
#include "../hcontainer_funcs.h"
#include "../hcontainer.h"

// test_size is the number of atoms in the unitcell
// modify test_size to test different size of unitcell
int test_size = 10;
int test_nw = 10;

/**
 * Unit test of folding_HR function
 * folding_HR is a function to calculate the Hk matrix with specific k vector
 * in this test, we test the correctness and time consuming of the function.
*/
class FoldingTest : public ::testing::Test
{
  protected:
    void SetUp() override
    {
        // set up a unitcell, with one element and three atoms, each atom has 2 orbitals
        ucell.ntype = 1;
        ucell.nat = test_size;
        ucell.atoms = new Atom[ucell.ntype];
        ucell.iat2it = new int[ucell.nat];
        for (int iat = 0; iat < ucell.nat; iat++)
        {
            ucell.iat2it[iat] = 0;
        }
        ucell.atoms[0].nw = test_nw;

        // set up a HContainer with ucell
        HR = new hamilt::HContainer<std::complex<double>>(ucell);
    }

    void TearDown() override
    {
        delete HR;
        delete[] ucell.atoms;
        delete[] ucell.iat2it;
    }

    UnitCell ucell;
    hamilt::HContainer<std::complex<double>>* HR;
};

// using TEST_F to test folding_HR
TEST_F(FoldingTest, folding_HR_cd2cd)
{
    // fill HR with constant value
    for (int i = 0; i < HR->size_atom_pairs(); i++)
    {
        std::complex<double>* ptr1 = HR->get_atom_pair(i).get_HR_values(1, 1, 1).get_pointer();
        std::complex<double>* ptr = HR->get_atom_pair(i).get_HR_values(0,0,0).get_pointer();
        for (int j = 0; j < HR->get_atom_pair(i).get_size(); j++)
        {
            ptr[j] = std::complex<double>(1.0, 1.0);
            ptr1[j] = std::complex<double>(2.0, 2.0);
        }
    }
    std::vector<std::complex<double>> hk(test_size * test_nw * test_size * test_nw, std::complex<double>(0.0, 0.0));
    ModuleBase::Vector3<double> kvec_d_in(0.1, 0.2, 0.3);
    hamilt::folding_HR(*HR, hk.data(), kvec_d_in, test_size * test_nw, 0);
    for (int i = 0; i < test_size * test_nw * test_size * test_nw; i++)
    {
        EXPECT_NEAR(hk[i].real(), 0.55753651583505226, 1e-10);
        EXPECT_NEAR(hk[i].imag(), -1.7936044933348412, 1e-10);
        //std::cout<<__FILE__<<__LINE__<<" hk["<<i<<"] = "<<hk[i]<<std::endl;
    }
    EXPECT_EQ(1, 1);
}