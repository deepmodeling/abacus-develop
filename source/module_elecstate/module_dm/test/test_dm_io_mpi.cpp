#include <chrono>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "module_elecstate/module_dm/density_matrix.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "prepare_unitcell.h"

// mock functions
#ifdef __LCAO
InfoNonlocal::InfoNonlocal()
{
}
InfoNonlocal::~InfoNonlocal()
{
}
LCAO_Orbitals::LCAO_Orbitals()
{
}
LCAO_Orbitals::~LCAO_Orbitals()
{
}
#endif
Magnetism::Magnetism()
{
    this->tot_magnetization = 0.0;
    this->abs_magnetization = 0.0;
    this->start_magnetization = nullptr;
}
Magnetism::~Magnetism()
{
    delete[] this->start_magnetization;
}
// mocke functions

/************************************************
 *  unit test of DensityMatrix constructor
 ***********************************************/

/**
 * This unit test construct a DensityMatrix object
 */

// test_size is the number of atoms in the unitcell
// modify test_size to test different size of unitcell
int test_size = 2;
int test_nw = 13;

class DMTest : public testing::Test
{
  protected:
    Parallel_Orbitals* paraV;
    int dsize;
    int my_rank = 0;
    UnitCell* ucell;
    UcellTestPrepare utp = UcellTestLib["Si"];
    std::vector<ModuleBase::ComplexMatrix> DMK;
    K_Vectors* kv = nullptr;
    // nw is the number of orbitals of each atom
    // it should container ucell.nat elements
    std::vector<int> nw = {13};
    int nks = 2;
    int nlocal = 0;
    void SetUp() override
    {
#ifdef __MPI
        // MPI parallel settings
        MPI_Comm_size(MPI_COMM_WORLD, &dsize);
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#endif
        // initalize a unitcell
        ucell = utp.SetUcellInfo(nw, nlocal);
        ucell->set_iat2iwt(1);
        // initalize a DMK
        DMK.reserve(nks);
        ModuleBase::ComplexMatrix zero_dmk(nlocal, nlocal);
        for (int ik = 0; ik < nks; ik++)
        {
            DMK.push_back(zero_dmk);
        }
        // initalize a kvectors
        kv = new K_Vectors;
        kv->nks = nks;
        kv->kvec_d.resize(nks);
        kv->kvec_d[1].x = 0.5;
        // set paraV
        init_parav();
    }

    void TearDown()
    {
        DMK.clear();
        delete kv;
        delete paraV;
    }

#ifdef __MPI
    void init_parav()
    {
        int global_row = test_size * test_nw;
        int global_col = test_size * test_nw;
        std::ofstream ofs_running;
        paraV = new Parallel_Orbitals();
        paraV->set_block_size(2 /* nb_2d set to be 2*/);
        paraV->set_proc_dim(dsize, 0);
        paraV->mpi_create_cart(MPI_COMM_WORLD);
        paraV->set_local2global(global_row, global_col, ofs_running, ofs_running);
        int lr = paraV->get_row_size();
        int lc = paraV->get_col_size();
        paraV->set_desc(global_row, global_col, lr);
        paraV->set_global2local(global_row, global_col, true, ofs_running);
        paraV->set_atomic_trace(ucell->get_iat2iwt(), test_size, global_row);
    }
#else
    void init_parav()
    {
    }
#endif
};

TEST_F(DMTest, DMConstructor1)
{
    // 
    // construct DM
    std::cout << paraV->nrow << paraV->ncol << std::endl;
    elecstate::DensityMatrix<double> DM(kv, paraV);
    // read DMK
    DM.set_DMK_files("./support/");
    // write DMK
    DM.output_DMK("./support/output");
    // construct a new DM
    elecstate::DensityMatrix<double> DM1(kv,paraV);
    DM1.set_DMK_files("./support/output");
    // compare DMK1 with DMK
    EXPECT_NEAR(DM.get_DMK(0,0,0),DM1.get_DMK(0,0,0),1e-6);
	EXPECT_NEAR(DM.get_DMK(1,25,25),DM1.get_DMK(1,25,25),1e-6);
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