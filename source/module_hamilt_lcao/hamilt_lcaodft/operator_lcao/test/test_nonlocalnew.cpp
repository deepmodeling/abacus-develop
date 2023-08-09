#include "gtest/gtest.h"
#include "../nonlocal_new.h"


//---------------------------------------
// Unit test of NonlocalNew class
// NonlocalNew is a derivative class of Operator, it is used to calculate the kinetic matrix
// It use HContainer to store the real space HR matrix
// In this test, we test the correctness and time consuming of 3 functions in NonlocalNew class
// - initialize_HR() called in constructor
// - contributeHR()
// - contributeHk()
// - HR(double) and SK(complex<double>) are tested in constructHRd2cd
// - HR(double) and SK(double) are tested in constructHRd2d
//---------------------------------------

// test_size is the number of atoms in the unitcell
// modify test_size to test different size of unitcell
int test_size = 2;
int test_nw = 5;
class NonlocalNewTest : public ::testing::Test
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
        ucell.infoNL.Beta = new Numerical_Nonlocal[ucell.ntype];
        ucell.nat = test_size;
        ucell.atoms = new Atom[ucell.ntype];
        ucell.iat2it = new int[ucell.nat];
        ucell.iat2ia = new int[ucell.nat];
        ucell.atoms[0].tau = new ModuleBase::Vector3<double>[ucell.nat];
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
        ucell.atoms[0].iw2l = new int[test_nw];
        ucell.atoms[0].iw2m = new int[test_nw];
        ucell.atoms[0].iw2n = new int[test_nw];
        for (int iw = 0; iw < test_nw; ++iw)
        {
            ucell.atoms[0].iw2l[iw] = 0;
            ucell.atoms[0].iw2m[iw] = 0;
            ucell.atoms[0].iw2n[iw] = 0;
        }
        ucell.atoms[0].ncpp.dion.create(5, 5);
        ucell.atoms[0].ncpp.dion.zero_out();
        ucell.atoms[0].ncpp.d_so.create(4, 5, 5);
        ucell.atoms[0].ncpp.d_so.zero_out();
        ucell.atoms[0].ncpp.non_zero_count_soc[0] = 5;
        ucell.atoms[0].ncpp.non_zero_count_soc[3] = 5;
        ucell.atoms[0].ncpp.index1_soc[0] = new int[5];
        ucell.atoms[0].ncpp.index2_soc[0] = new int[5];
        ucell.atoms[0].ncpp.index1_soc[3] = new int[5];
        ucell.atoms[0].ncpp.index2_soc[3] = new int[5];
        for(int i = 0; i < 5; ++i)
        {
            ucell.atoms[0].ncpp.dion(i, i) = 1.0;
            ucell.atoms[0].ncpp.d_so(0, i, i) = std::complex<double>(2.0, 0.0);
            ucell.atoms[0].ncpp.d_so(3, i, i) = std::complex<double>(2.0, 0.0);
            ucell.atoms[0].ncpp.index1_soc[0][i] = i;
            ucell.atoms[0].ncpp.index2_soc[0][i] = i;
            ucell.atoms[0].ncpp.index1_soc[3][i] = i;
            ucell.atoms[0].ncpp.index2_soc[3][i] = i;
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
        delete[] ucell.atoms[0].tau;
        delete[] ucell.atoms[0].iw2l;
        delete[] ucell.atoms[0].iw2m;
        delete[] ucell.atoms[0].iw2n;
        delete[] ucell.atoms;
        delete[] ucell.iat2it;
        delete[] ucell.iat2ia;
    }

#ifdef __MPI
    void init_parav()
    {
        int global_row = test_size * test_nw;
        int global_col = test_size * test_nw;
        std::ofstream ofs_running;
        paraV = new Parallel_Orbitals();
        paraV->set_block_size(2/* nb_2d set to be 2*/);
        paraV->set_proc_dim(dsize, 0);
        paraV->mpi_create_cart(MPI_COMM_WORLD);
        paraV->set_local2global(global_row, global_col, ofs_running, ofs_running);
        int lr = paraV->get_row_size();
        int lc = paraV->get_col_size();
        paraV->set_desc(global_row, global_col, lr);
        paraV->set_global2local(global_row, global_col, true, ofs_running);
        paraV->set_atomic_trace(ucell.get_iat2iwt(), test_size, global_row);
    }
#else
    void init_parav()
    {}
#endif

    UnitCell ucell;
    hamilt::HContainer<double>* HR;
    Parallel_Orbitals *paraV;

    int dsize;
    int my_rank = 0;
};

// using TEST_F to test NonlocalNew
TEST_F(NonlocalNewTest, constructHRd2d)
{
    std::vector<ModuleBase::Vector3<double>> kvec_d_in(1, ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
    std::vector<double> hk(paraV->get_row_size() * paraV->get_col_size(), 0.0);
    Grid_Driver gd(0,0,0);
    // check some input values
    EXPECT_EQ(ucell.infoNL.Beta[0].get_rcut_max(), 1.0);
    EXPECT_EQ(LCAO_Orbitals::get_const_instance().Phi[0].getRcut(), 1.0);
    hamilt::NonlocalNew<hamilt::OperatorLCAO<double>, double> op(
        nullptr, 
        kvec_d_in, 
        HR, 
        hk.data(), 
        &ucell, 
        &gd,
        paraV
    );
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
            EXPECT_EQ(tmp.get_pointer(0)[i], 10.0);
        }
    }
    // calculate SK
    op.contributeHk(0);
    // check the value of SK
    for (int i = 0; i < paraV->get_row_size() * paraV->get_col_size(); ++i)
    {
        EXPECT_EQ(hk[i], 10.0);
    }
}

TEST_F(NonlocalNewTest, constructHRd2cd)
{
    std::vector<ModuleBase::Vector3<double>> kvec_d_in(2, ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
    kvec_d_in[1] = ModuleBase::Vector3<double>(0.1, 0.2, 0.3);
    std::vector<std::complex<double>> hk(paraV->get_row_size() * paraV->get_col_size(), std::complex<double>(0.0, 0.0));
    Grid_Driver gd(0,0,0);
    hamilt::NonlocalNew<hamilt::OperatorLCAO<std::complex<double>>, double> op(
        nullptr, 
        kvec_d_in, 
        HR, 
        hk.data(), 
        &ucell, 
        &gd,
        paraV
    );
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
            EXPECT_EQ(tmp.get_pointer(0)[i], 10.0);
        }
    }
    // calculate SK for gamma point
    op.contributeHk(0);
    // check the value of SK of gamma point
    for (int i = 0; i < paraV->get_row_size() * paraV->get_col_size(); ++i)
    {
        EXPECT_EQ(hk[i].real(), 10.0);
        EXPECT_EQ(hk[i].imag(), 0.0);
    }
    // calculate HK for k point
    hk.assign(paraV->get_row_size() * paraV->get_col_size(), std::complex<double>(0.0, 0.0) );
    op.contributeHk(1);
    // check the value of HK
    for (int i = 0; i < paraV->get_row_size() * paraV->get_col_size(); ++i)
    {
        EXPECT_NEAR(hk[i].real(), 10.0, 1e-10);
        EXPECT_NEAR(hk[i].imag(), 0.0, 1e-10);
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