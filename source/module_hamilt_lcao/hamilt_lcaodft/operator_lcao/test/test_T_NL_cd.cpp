#include "gtest/gtest.h"
#include "../ekinetic_new.h"
#include "../nonlocal_new.h"
#include <chrono>

//---------------------------------------
// Unit test of EkineticNew + NonlocalNew class
// EkineticNew and NonlocalNew are derivative classes of Operator, used to calculate the T+VNL matrix
// It use HContainer to store the real space HR matrix
// In this test, we test the correctness and time consuming of 6 functions in T+VNL class
// - initialize_HR() called in two constructors
// - contributeHR() in two Operators
// - contributeHk() in two Operators
// - HR(complex<double>) and HK(complex<double>) are tested in constructHRd2d
//---------------------------------------

// test_size is the number of atoms in the unitcell
// modify test_size to test different size of unitcell
int test_size = 10;
int test_nw = 10;
class TNLTest : public ::testing::Test
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
        // for NonlocalNew
        ucell.infoNL.Beta = new Numerical_Nonlocal[ucell.ntype];
        ucell.atoms[0].ncpp.d_real.create(5, 5);
        ucell.atoms[0].ncpp.d_real.zero_out();
        ucell.atoms[0].ncpp.d_so.create(4, 5, 5);
        ucell.atoms[0].ncpp.d_so.zero_out();
        ucell.atoms[0].ncpp.non_zero_count_soc[0] = 5;
        ucell.atoms[0].ncpp.non_zero_count_soc[1] = 0;
        ucell.atoms[0].ncpp.non_zero_count_soc[2] = 0;
        ucell.atoms[0].ncpp.non_zero_count_soc[3] = 5;
        ucell.atoms[0].ncpp.index1_soc[0] = new int[5];
        ucell.atoms[0].ncpp.index2_soc[0] = new int[5];
        ucell.atoms[0].ncpp.index1_soc[3] = new int[5];
        ucell.atoms[0].ncpp.index2_soc[3] = new int[5];
        for(int i = 0; i < 5; ++i)
        {
            ucell.atoms[0].ncpp.d_real(i, i) = 1.0;
            ucell.atoms[0].ncpp.d_so(0, i, i) = std::complex<double>(2.0, 0.0);
            ucell.atoms[0].ncpp.d_so(3, i, i) = std::complex<double>(2.0, 0.0);
            ucell.atoms[0].ncpp.index1_soc[0][i] = i;
            ucell.atoms[0].ncpp.index2_soc[0][i] = i;
            ucell.atoms[0].ncpp.index1_soc[3][i] = i;
            ucell.atoms[0].ncpp.index2_soc[3][i] = i;
        }
        // end of set up a unitcell
        init_parav();
        // set up a HContainer with ucell
        HR = new hamilt::HContainer<std::complex<double>>(paraV);
    }

    void TearDown() override
    {
        delete HR;
        delete paraV;
        delete[] ucell.atoms[0].tau;
        delete[] ucell.atoms[0].iw2l;
        delete[] ucell.atoms[0].iw2m;
        delete[] ucell.atoms[0].iw2n;
        delete[] ucell.atoms[0].ncpp.index1_soc[0];
        delete[] ucell.atoms[0].ncpp.index2_soc[0];
        delete[] ucell.atoms[0].ncpp.index1_soc[3];
        delete[] ucell.atoms[0].ncpp.index2_soc[3];
        delete[] ucell.atoms;
        delete[] ucell.iat2it;
        delete[] ucell.iat2ia;
        delete[] ucell.infoNL.Beta;

    }

#ifdef __MPI
    void init_parav()
    {
        int global_row = test_size * test_nw * 2;
        int global_col = test_size * test_nw * 2;
        std::ofstream ofs_running;
        paraV = new Parallel_Orbitals();
        paraV->set_block_size(20/* nb_2d set to be 2*/);
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
    hamilt::HContainer<std::complex<double>>* HR;
    Parallel_Orbitals *paraV;

    int dsize;
    int my_rank = 0;
};

TEST_F(TNLTest, testTVNLcd2cd)
{
    int npol = ucell.get_npol();
    std::vector<ModuleBase::Vector3<double>> kvec_d_in(2, ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
    kvec_d_in[1] = ModuleBase::Vector3<double>(0.1, 0.2, 0.3);
    std::vector<std::complex<double>> hk(paraV->get_row_size() * paraV->get_col_size(), std::complex<double>(0.0, 0.0));
    Grid_Driver gd(0,0,0);
    std::chrono::high_resolution_clock::time_point start_time = std::chrono::high_resolution_clock::now();
    hamilt::Operator<std::complex<double>> *op = new hamilt::EkineticNew<hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>>(
        nullptr, 
        kvec_d_in, 
        HR, 
        &hk, 
        &ucell, 
        &gd,
        paraV
    );
    hamilt::Operator<std::complex<double>> *op1 = new hamilt::NonlocalNew<hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>>(
        nullptr, 
        kvec_d_in, 
        HR, 
        &hk, 
        &ucell, 
        &gd,
        paraV
    );
    // merge two Operators to a chain
    op->add(op1);
    std::chrono::high_resolution_clock::time_point end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time0 = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    start_time = std::chrono::high_resolution_clock::now();
    // calculate HR and folding HK for gamma point
    op->init(0);
    end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time1 = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);

    // check the value of HR
    double result_ref = test_size * 10;
    for (int iap = 0; iap < HR->size_atom_pairs(); ++iap)
    {
        hamilt::AtomPair<std::complex<double>>& tmp = HR->get_atom_pair(iap);
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
                    EXPECT_EQ(tmp.get_pointer(0)[i].real(), 2.0);
                    EXPECT_EQ(tmp.get_pointer(0)[i].imag(), 0.0);
                    EXPECT_EQ(tmp.get_pointer(1)[i].real(), result_ref);
                    EXPECT_EQ(tmp.get_pointer(1)[i].imag(), 0.0);
                }
                else
                {
                    EXPECT_EQ(tmp.get_pointer(0)[i].real(), 0.0);
                    EXPECT_EQ(tmp.get_pointer(0)[i].imag(), 0.0);
                    EXPECT_EQ(tmp.get_pointer(1)[i].real(), 0.0);
                    EXPECT_EQ(tmp.get_pointer(1)[i].imag(), 0.0);
                }
                ++i;
            }
        }
    }
    // check the value of HK of gamma point
    result_ref += 2.0;
    int i = 0;
    for ( int irow = 0; irow < paraV->get_row_size(); ++irow)
    {
        for ( int icol = 0; icol < paraV->get_col_size(); ++icol)
        {
            if (irow%npol == icol%npol)
            {
                EXPECT_NEAR(hk[i].real(), result_ref, 1e-10);
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
    // calculate HK for k point
    start_time = std::chrono::high_resolution_clock::now();
    hk.assign(paraV->get_row_size() * paraV->get_col_size(), std::complex<double>(0.0, 0.0) );
    op->init(1);
    end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time2 = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    std::cout << "Test terms:   " <<std::setw(15)<< "constructor" <<std::setw(15)<< "init(HR+HK)" <<std::setw(15)<< "2nd-init(HK)" << std::endl;
    std::cout << "Elapsed time: " <<std::setw(15)<< elapsed_time0.count()<<std::setw(15)<<elapsed_time1.count()<<std::setw(15)<<elapsed_time2.count() << " seconds." << std::endl;
    // check the value of HK
    double result_ref1 = -1.6180339887498931 + test_size * 10;
    double result_ref2 = -1.1755705045849467;
    i = 0;
    for ( int irow = 0; irow < paraV->get_row_size(); ++irow)
    {
        for ( int icol = 0; icol < paraV->get_col_size(); ++icol)
        {
            if (irow%npol == icol%npol)
            {
                EXPECT_NEAR(hk[i].real(), result_ref1, 1e-10);
                EXPECT_NEAR(hk[i].imag(), result_ref2, 1e-10);
            }
            else
            {
                EXPECT_NEAR(hk[i].real(), 0.0, 1e-10);
                EXPECT_NEAR(hk[i].imag(), 0.0, 1e-10);
            }
            ++i;
        }
    }
    delete op;
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