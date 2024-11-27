#include "../ekinetic_new.h"

#include "gtest/gtest.h"

//---------------------------------------
// Unit test of EkineticNew class
// EkineticNew is a derivative class of Operator, it is used to calculate the kinetic matrix
// It use HContainer to store the real space HR matrix
// In this test, we test the correctness and time consuming of 3 functions in EkineticNew class
// - initialize_HR() called in constructor
// - contributeHR()
// - contributeHk()
// - HR(double) and SK(complex<double>) are tested in constructHRd2cd
// - HR(double) and SK(double) are tested in constructHRd2d
//---------------------------------------

// test_size is the number of atoms in the unitcell
// modify test_size to test different size of unitcell
int test_size = 10;
int test_nw = 10;
class EkineticNewTest : public ::testing::Test
{
  protected:
    void SetUp() override
    {
#ifdef __MPI
        // MPI parallel settings
        MPI_Comm_size(MPI_COMM_WORLD, &dsize);
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#endif

        // set up a unitcell, with one element and test_size atoms, each atom has test_nw orbitals
        ucell.ntype = 1;
        ucell.nat = test_size;
        ucell.atoms = new Atom[ucell.ntype];
        ucell.iat2it = new int[ucell.nat];
        ucell.iat2ia = new int[ucell.nat];
        ucell.atoms[0].tau.resize(ucell.nat);
        ucell.itia2iat.create(ucell.ntype, ucell.nat);
        for (int iat = 0; iat < ucell.nat; iat++)
        {
            ucell.iat2it[iat] = 0;
            ucell.iat2ia[iat] = iat;
            ucell.atoms[0].tau[iat] = ModuleBase::Vector3<double>(0.0, 0.0, 0.0);
            ucell.itia2iat(0, iat) = iat;
        }
        ucell.atoms[0].na = test_size;
        ucell.atoms[0].nw = test_nw;
        ucell.atoms[0].iw2l.resize(test_nw);
        ucell.atoms[0].iw2m.resize(test_nw);
        ucell.atoms[0].iw2n.resize(test_nw);
        for (int iw = 0; iw < test_nw; ++iw)
        {
            ucell.atoms[0].iw2l[iw] = 0;
            ucell.atoms[0].iw2m[iw] = 0;
            ucell.atoms[0].iw2n[iw] = 0;
        }
        ucell.set_iat2iwt(1);
        init_parav();
        // set up a HContainer with ucell
        HR = new hamilt::HContainer<double>(paraV);
    }

    void TearDown() override
    {
        delete HR;
        delete paraV;
        delete[] ucell.atoms;
        delete[] ucell.iat2it;
        delete[] ucell.iat2ia;
    }

#ifdef __MPI
    void init_parav()
    {
        int nb = 10;
        int global_row = test_size * test_nw;
        int global_col = test_size * test_nw;
        std::ofstream ofs_running;
        paraV = new Parallel_Orbitals();
        paraV->init(global_row, global_col, nb, MPI_COMM_WORLD);
        paraV->set_atomic_trace(ucell.get_iat2iwt(), test_size, global_row);
    }
#else
    void init_parav()
    {
    }
#endif

    UnitCell ucell;
    hamilt::HContainer<double>* HR;
    Parallel_Orbitals* paraV;
    TwoCenterIntegrator intor_;

    int dsize;
    int my_rank = 0;
};

// using TEST_F to test EkineticNew
TEST_F(EkineticNewTest, constructHRd2d)
{
    std::vector<ModuleBase::Vector3<double>> kvec_d_in(1, ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
    hamilt::HS_Matrix_K<double> hsk(paraV, true);
    hsk.set_zero_hk();
    Grid_Driver gd(0, 0);
    hamilt::EkineticNew<hamilt::OperatorLCAO<double, double>>
        op(&hsk, kvec_d_in, HR, &ucell, {1.0}, &gd, &intor_);
    op.contributeHR();
    // check the value of HR
    for (int iap = 0; iap < HR->size_atom_pairs(); ++iap)
    {
        hamilt::AtomPair<double>& tmp = HR->get_atom_pair(iap);
        int iat1 = tmp.get_atom_i();
        int iat2 = tmp.get_atom_j();
        auto indexes1 = paraV->get_indexes_row(iat1);
        auto indexes2 = paraV->get_indexes_col(iat2);
        int nwt = indexes1.size() * indexes2.size();
        for (int i = 0; i < nwt; ++i)
        {
            EXPECT_EQ(tmp.get_pointer(0)[i], 1.0);
        }
    }
    // calculate HK
    op.contributeHk(0);
    // check the value of HK
    double* hk = hsk.get_hk();
    for (int i = 0; i < paraV->get_row_size() * paraV->get_col_size(); ++i)
    {
        EXPECT_EQ(hk[i], 1.0);
    }
    // calculate HR again
    op.contributeHR();
    // check the value of HR
    for (int iap = 0; iap < HR->size_atom_pairs(); ++iap)
    {
        hamilt::AtomPair<double>& tmp = HR->get_atom_pair(iap);
        int iat1 = tmp.get_atom_i();
        int iat2 = tmp.get_atom_j();
        auto indexes1 = paraV->get_indexes_row(iat1);
        auto indexes2 = paraV->get_indexes_col(iat2);
        int nwt = indexes1.size() * indexes2.size();
        for (int i = 0; i < nwt; ++i)
        {
            EXPECT_EQ(tmp.get_pointer(0)[i], 2.0);
        }
    }
}

TEST_F(EkineticNewTest, constructHRd2cd)
{
    std::vector<ModuleBase::Vector3<double>> kvec_d_in(2, ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
    kvec_d_in[1] = ModuleBase::Vector3<double>(0.1, 0.2, 0.3);
    hamilt::HS_Matrix_K<std::complex<double>> hsk(paraV, true);
    hsk.set_zero_hk();
    Grid_Driver gd(0, 0);
    hamilt::EkineticNew<hamilt::OperatorLCAO<std::complex<double>, double>>
        op(&hsk, kvec_d_in, HR, &ucell, {1.0}, &gd, &intor_);
    op.contributeHR();
    // check the value of HR
    for (int iap = 0; iap < HR->size_atom_pairs(); ++iap)
    {
        hamilt::AtomPair<double>& tmp = HR->get_atom_pair(iap);
        int iat1 = tmp.get_atom_i();
        int iat2 = tmp.get_atom_j();
        auto indexes1 = paraV->get_indexes_row(iat1);
        auto indexes2 = paraV->get_indexes_col(iat2);
        int nwt = indexes1.size() * indexes2.size();
        for (int i = 0; i < nwt; ++i)
        {
            EXPECT_EQ(tmp.get_pointer(0)[i], 1.0);
        }
    }
    // calculate HK for gamma point
    op.contributeHk(0);
    auto* hk = hsk.get_hk();
    // check the value of HK of gamma point
    for (int i = 0; i < paraV->get_row_size() * paraV->get_col_size(); ++i)
    {
        EXPECT_EQ(hk[i].real(), 1.0);
        EXPECT_EQ(hk[i].imag(), 0.0);
    }
    // calculate HK for k point
    hsk.set_zero_hk();
    op.contributeHk(1);
    // check the value of HK
    for (int i = 0; i < paraV->get_row_size() * paraV->get_col_size(); ++i)
    {
        EXPECT_NEAR(hk[i].real(), -1.6180339887498945 / 2, 1e-10);
        EXPECT_NEAR(hk[i].imag(), -1.1755705045849467 / 2, 1e-10);
    }
}

int main(int argc, char** argv)
{
#ifdef __MPI
    MPI_Init(&argc, &argv);
#endif
    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();
#ifdef __MPI
    MPI_Finalize();
#endif
    return result;
}
