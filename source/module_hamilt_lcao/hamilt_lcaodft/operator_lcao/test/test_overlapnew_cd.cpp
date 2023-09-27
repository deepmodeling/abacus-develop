#include "gtest/gtest.h"
#include "../overlap_new.h"


//---------------------------------------
// Unit test of OverlapNew class
// OverlapNew is a derivative class of Operator, it is used to calculate the overlap matrix
// It use HContainer to store the real space SR matrix
// In this test, we test the correctness and time consuming of 3 functions in OverlapNew class
// - initialize_SR() called in constructor
// - contributeHR()
// - contributeHk()
// - SR(complex<double>) and SK(complex<double>) are tested in constructHRd2d
//---------------------------------------

// test_size is the number of atoms in the unitcell
// modify test_size to test different size of unitcell
int test_size = 10;
int test_nw = 10;
class OverlapNewTest : public ::testing::Test
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
        ucell.set_iat2iwt(2);
        init_parav();
        // set up a HContainer with ucell
        SR = new hamilt::HContainer<std::complex<double>>(paraV);
    }

    void TearDown() override
    {
        delete SR;
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
        int global_row = test_size * test_nw * 2;
        int global_col = test_size * test_nw * 2;
        std::ofstream ofs_running;
        paraV = new Parallel_Orbitals();
        paraV->set_block_size(10/* nb_2d set to be 2*/);
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
    hamilt::HContainer<std::complex<double>>* SR;
    Parallel_Orbitals *paraV;

    int dsize;
    int my_rank = 0;
};

TEST_F(OverlapNewTest, constructHRcd2cd)
{
    int npol = ucell.get_npol();
    std::vector<ModuleBase::Vector3<double>> kvec_d_in(2, ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
    kvec_d_in[1] = ModuleBase::Vector3<double>(0.1, 0.2, 0.3);
    std::vector<std::complex<double>> hk(paraV->get_row_size() * paraV->get_col_size(), std::complex<double>(0.0, 0.0));
    Grid_Driver gd(0,0,0);
    hamilt::OverlapNew<hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>> op(
        nullptr, 
        kvec_d_in, 
        nullptr,
        nullptr,
        SR, 
        &hk, 
        &ucell, 
        &gd,
        paraV
    );
    op.contributeHR();
    // check the value of SR
    for (int iap = 0; iap < SR->size_atom_pairs(); ++iap)
    {
        hamilt::AtomPair<std::complex<double>>& tmp = SR->get_atom_pair(iap);
        int iat1 = tmp.get_atom_i();
        int iat2 = tmp.get_atom_j();
        auto indexes1 = paraV->get_indexes_row(iat1);
        auto indexes2 = paraV->get_indexes_col(iat2);
        int i = 0;
        for (int mu = 0; mu < indexes1.size(); ++mu)
        {
            for(int nu = 0; nu < indexes2.size(); ++nu)
            {
                if(mu % npol == nu % npol)
                {
                    EXPECT_EQ(tmp.get_pointer(0)[i].real(), 1.0);
                    EXPECT_EQ(tmp.get_pointer(0)[i].imag(), 0.0);
                }
                else
                {
                    EXPECT_EQ(tmp.get_pointer(0)[i].real(), 0.0);
                    EXPECT_EQ(tmp.get_pointer(0)[i].imag(), 0.0);
                }
                ++i;
            }
        }
    }
    // calculate SK for gamma point
    op.contributeHk(0);
    // check the value of SK of gamma point
    int i = 0;
    for ( int irow = 0; irow < paraV->get_row_size(); ++irow)
    {
        for ( int icol = 0; icol < paraV->get_col_size(); ++icol)
        {
            if (irow%npol == icol%npol)
            {
                EXPECT_NEAR(hk[i].real(), 1.0, 1e-10);
                EXPECT_NEAR(hk[i].imag(), 0.0, 1e-10);
            }
            else
            {
                EXPECT_NEAR(hk[i].real(), 0.0, 1e-10);
                EXPECT_NEAR(hk[i].imag(), 0.0, 1e-10);
            }
            ++i;
        }
    }
    // calculate SK for k point
    hk.assign(paraV->get_row_size() * paraV->get_col_size(), std::complex<double>(0.0, 0.0) );
    op.contributeHk(1);
    // check the value of SK
    i = 0;
    for ( int irow = 0; irow < paraV->get_row_size(); ++irow)
    {
        for ( int icol = 0; icol < paraV->get_col_size(); ++icol)
        {
            if (irow%npol == icol%npol)
            {
                EXPECT_NEAR(hk[i].real(), -0.80901699437494723, 1e-10);
                EXPECT_NEAR(hk[i].imag(), -0.58778525229247336, 1e-10);
            }
            else
            {
                EXPECT_NEAR(hk[i].real(), 0.0, 1e-10);
                EXPECT_NEAR(hk[i].imag(), 0.0, 1e-10);
            }
            ++i;
        }
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