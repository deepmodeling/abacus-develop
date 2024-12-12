#include <chrono>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "module_cell/unitcell.h"
#include "module_elecstate/module_dm/density_matrix.h"
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

#include "module_cell/klist.h"

K_Vectors::K_Vectors()
{
}

K_Vectors::~K_Vectors()
{
}

#include "module_cell/module_neighbor/sltk_grid_driver.h"
// mock find_atom() function
void Grid_Driver::Find_atom(const UnitCell& ucell,
                            const ModuleBase::Vector3<double>& tau,
                            const int& T,
                            const int& I,
                            AdjacentAtomInfo* adjs) const
{
}
Grid::Grid(const int& test_grid_in) : test_grid(test_grid_in)
{
}
Grid::~Grid()
{
}
Grid_Driver::Grid_Driver(const int& test_d_in,const int& test_grid_in)
    : Grid(test_grid_in), test_deconstructor(test_d_in){}
Grid_Driver::~Grid_Driver()
{
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
        // initalize a kvectors
        kv = new K_Vectors;
        kv->set_nks(nks);
        kv->kvec_d.resize(nks);
        kv->kvec_d[1].x = 0.5;
        // set paraV
        init_parav();
    }

    void TearDown() override
    {
        DMK.clear();
        delete kv;
        delete paraV;
    }

#ifdef __MPI
    void init_parav()
    {
        int nb = 2;
        int global_row = test_size * test_nw;
        int global_col = test_size * test_nw;
        std::ofstream ofs_running;
        paraV = new Parallel_Orbitals();
        paraV->init(global_row, global_col, nb, MPI_COMM_WORLD);
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
    int nspin = 1;
    // construct DM
    std::cout << paraV->nrow << paraV->ncol << std::endl;
    elecstate::DensityMatrix<double, double> DM(paraV, nspin, kv->kvec_d, kv->get_nks());
    // read DMK
    std::string directory = "./support/";
    for (int is = 1; is <= nspin; ++is)
    {
        for (int ik = 0; ik < kv->get_nks() / nspin; ++ik)
        {
            DM.read_DMK(directory, is, ik);
        }
    }
    // write DMK
    directory = "./support/output";
    for (int is = 1; is <= nspin; ++is)
    {
        for (int ik = 0; ik < kv->get_nks() / nspin; ++ik)
        {
            DM.write_DMK(directory, is, ik);
        }
    }
    // construct a new DM
    elecstate::DensityMatrix<double, double> DM1(paraV, nspin, kv->kvec_d, kv->get_nks());
    directory = "./support/output";
    for (int is = 1; is <= nspin; ++is)
    {
        for (int ik = 0; ik < kv->get_nks() / nspin; ++ik)
        {
            DM1.read_DMK(directory, is, ik);
        }
    }
    // compare DMK1 with DMK
    EXPECT_NEAR(DM.get_DMK(1, 0, 0, 0), DM1.get_DMK(1, 0, 0, 0), 1e-6);
    EXPECT_NEAR(DM.get_DMK(1, 1, 25, 25), DM1.get_DMK(1, 1, 25, 25), 1e-6);
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
