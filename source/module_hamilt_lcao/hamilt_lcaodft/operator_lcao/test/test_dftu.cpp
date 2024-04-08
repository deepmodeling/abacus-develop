#include "gtest/gtest.h"
#include "../dftu_lcao.h"
#include <chrono> 

// mock of DFTU
#include "module_hamilt_lcao/module_dftu/dftu.h"
ModuleDFTU::DFTU::DFTU(){};
ModuleDFTU::DFTU::~DFTU(){};
namespace GlobalC
{
    ModuleDFTU::DFTU dftu;
}
const hamilt::HContainer<double>* tmp_DMR;
const hamilt::HContainer<double>* ModuleDFTU::DFTU::get_dmr(int ispin) const
{
    return tmp_DMR;
}



//---------------------------------------
// Unit test of DFTU class
// DFTU is a derivative class of Operator, it is used to calculate the kinetic matrix
// It use HContainer to store the real space HR matrix
// In this test, we test the correctness and time consuming of 3 functions in DFTU class
// - initialize_HR() called in constructor
// - contributeHR()
// - contributeHk()
// - HR(double) and SK(complex<double>) are tested in constructHRd2cd
// - HR(double) and SK(double) are tested in constructHRd2d
//---------------------------------------

// test_size is the number of atoms in the unitcell
// modify test_size to test different size of unitcell
int test_size = 10; 
int test_nw = 10;   // please larger than 5
class DFTUTest : public ::testing::Test
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
        ucell.lat0 = 1.0;
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
            ucell.atoms[0].iw2l[iw] = 2;
            ucell.atoms[0].iw2m[iw] = 0;
            ucell.atoms[0].iw2n[iw] = 0;
        }
        ucell.set_iat2iwt(1);
        init_parav();
        // set up a HContainer with ucell
        HR = new hamilt::HContainer<double>(ucell, paraV);
        // initialize DMR and set default values of 1.0
        DMR = new hamilt::HContainer<double>(*HR);
        tmp_DMR = DMR;

        //setting of DFTU
        GlobalC::dftu.locale.resize(test_size);
        for(int iat=0;iat<test_size;iat++)
        {
            GlobalC::dftu.locale[iat].resize(3);
            for(int l=0;l<3;l++)
            {
                GlobalC::dftu.locale[iat][l].resize(1);
                GlobalC::dftu.locale[iat][l][0].resize(2);
                GlobalC::dftu.locale[iat][l][0][0].create(2*l+1, 2*l+1);
                GlobalC::dftu.locale[iat][l][0][1].create(2*l+1, 2*l+1);
            }
        }
        GlobalC::dftu.U = &U_test;
        GlobalC::dftu.orbital_corr = &orbital_c_test;

        GlobalV::onsite_radius = 1.0;

    }

    void TearDown() override
    {
        delete HR;
        delete DMR;
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
    hamilt::HContainer<double>* HR;
    hamilt::HContainer<double>* DMR;
    Parallel_Orbitals *paraV;

    int dsize;
    int my_rank = 0;
    double U_test = 1.0;
    int orbital_c_test = 2;
};

// using TEST_F to test DFTU
TEST_F(DFTUTest, constructHRd2d)
{
    //test for nspin=1
    GlobalV::NSPIN = 1;
    std::vector<ModuleBase::Vector3<double>> kvec_d_in(1, ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
    std::vector<double> hk(paraV->get_row_size() * paraV->get_col_size(), 0.0);
    Grid_Driver gd(0,0,0);
    // check some input values
    EXPECT_EQ(LCAO_Orbitals::get_const_instance().Phi[0].getRcut(), 1.0);
    // reset HR and DMR
    const double factor = 1.0 / test_nw / test_nw / test_size / test_size;
    for(int i=0;i<DMR->get_nnr();i++)
    {
        DMR->get_wrapper()[i] = factor;
        HR->get_wrapper()[i] = 0.0;
    }
    std::chrono::high_resolution_clock::time_point start_time = std::chrono::high_resolution_clock::now();
    hamilt::DFTU<hamilt::OperatorLCAO<double, double>> op(
        nullptr, 
        kvec_d_in, 
        HR, 
        &hk, 
        ucell, 
        &gd,
        &GlobalC::dftu,
        *paraV
    );
    std::chrono::high_resolution_clock::time_point end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    start_time = std::chrono::high_resolution_clock::now();
    op.contributeHR();
    end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time1 = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    // check the occupations of dftu
    for(int iat=0;iat<test_size;iat++)
    {
        for(int icc=0;icc<25;icc++)
        {
            EXPECT_NEAR(GlobalC::dftu.locale[iat][2][0][0].c[icc], 0.5, 1e-10);
        }
    }
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
            EXPECT_NEAR(tmp.get_pointer(0)[i], -10.0*test_size, 1e-10);
        }
    }
    // calculate SK
    start_time = std::chrono::high_resolution_clock::now();
    op.contributeHk(0);
    end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time2 = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    // check the value of SK
    for (int i = 0; i < paraV->get_row_size() * paraV->get_col_size(); ++i)
    {
        EXPECT_NEAR(hk[i], -10.0*test_size, 1e-10);
    }
    std::cout << "Test terms:   " <<std::setw(15)<< "initialize_HR" <<std::setw(15)<< "contributeHR" <<std::setw(15)<< "contributeHk" << std::endl;
    std::cout << "Elapsed time: " <<std::setw(15)<< elapsed_time.count()<<std::setw(15)<<elapsed_time1.count()<<std::setw(15)<<elapsed_time2.count() << " seconds." << std::endl;
}

TEST_F(DFTUTest, constructHRd2cd)
{
    // test for nspin=2
    GlobalV::NSPIN = 2;
    std::vector<ModuleBase::Vector3<double>> kvec_d_in(2, ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
    std::vector<std::complex<double>> hk(paraV->get_row_size() * paraV->get_col_size(), std::complex<double>(0.0, 0.0));
    Grid_Driver gd(0,0,0);
    EXPECT_EQ(LCAO_Orbitals::get_const_instance().Phi[0].getRcut(), 1.0);
    // reset HR and DMR
    const double factor = 0.5 / test_nw / test_nw / test_size / test_size;
    for(int i=0;i<DMR->get_nnr();i++)
    {
        DMR->get_wrapper()[i] = factor;
        HR->get_wrapper()[i] = 0.0;
    }
    hamilt::DFTU<hamilt::OperatorLCAO<std::complex<double>, double>> op(
        nullptr, 
        kvec_d_in, 
        HR, 
        &hk, 
        ucell, 
        &gd,
        &GlobalC::dftu,
        *paraV
    );
    op.contributeHR();
    // check the occupations of dftu for spin-up
    for(int iat=0;iat<test_size;iat++)
    {
        for(int icc=0;icc<25;icc++)
        {
            EXPECT_NEAR(GlobalC::dftu.locale[iat][2][0][0].c[icc], 0.5, 1e-10);
        }
    }
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
            EXPECT_NEAR(tmp.get_pointer(0)[i], -10.0*test_size, 1e-10);
        }
    }
    // calculate HK for gamma point
    op.contributeHk(0);
    // check the value of HK of gamma point
    for (int i = 0; i < paraV->get_row_size() * paraV->get_col_size(); ++i)
    {
        EXPECT_NEAR(hk[i].real(), -10.0*test_size, 1e-10);
        EXPECT_NEAR(hk[i].imag(), 0.0, 1e-10);
    }
    // calculate spin-down hamiltonian
    op.contributeHR();
    // check the occupations of dftu for spin-down
    for(int iat=0;iat<test_size;iat++)
    {
        for(int icc=0;icc<25;icc++)
        {
            EXPECT_NEAR(GlobalC::dftu.locale[iat][2][0][1].c[icc], 0.5, 1e-10);
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