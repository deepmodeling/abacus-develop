#include <chrono>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "module_elecstate/module_dm/density_matrix.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"

/************************************************
 *  unit test of DensityMatrix constructor
 ***********************************************/

/**
 * This unit test construct a DensityMatrix object
 */

// test_size is the number of atoms in the unitcell
// modify test_size to test different size of unitcell
int test_size = 10;
int test_nw = 26;

class DMTest : public testing::Test
{
  protected:
    Parallel_Orbitals* paraV;
    int dsize;
    int my_rank = 0;
    UnitCell ucell;
    K_Vectors* kv = nullptr;
    int nks = 2;
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
        ucell.set_iat2iwt(1);
        init_parav();

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
        delete kv;
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
        paraV->set_block_size(2 /* nb_2d set to be 2*/);
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
    {
    }
#endif
};

TEST_F(DMTest, DMInit1)
{
    // construct DM
    std::cout << "dim0: " << paraV->dim0 << "    dim1:" << paraV->dim1 << std::endl;
    std::cout << "nrow: " << paraV->nrow << "    ncol:" << paraV->ncol << std::endl;
    elecstate::DensityMatrix<double> DM(kv, paraV);
    // initialize this->_DMR
    Grid_Driver gd(0, 0, 0);
    DM.init_DMR(&gd, &ucell);
    // compare
    EXPECT_EQ(DM.get_DMR_pointer()->size_atom_pairs(), test_size * test_size);
    EXPECT_EQ(DM.get_DMR_pointer()->get_atom_pair(2, 2).get_atom_i(), 2);
    EXPECT_EQ(DM.get_DMR_pointer()->get_atom_pair(2, 2).get_atom_j(), 2);
    EXPECT_EQ(DM.get_DMR_pointer()->get_atom_pair(2, 2).get_row_size(), paraV->get_row_size(2));
    EXPECT_EQ(DM.get_DMR_pointer()->get_atom_pair(2, 2).get_col_size(), paraV->get_col_size(2));
}

TEST_F(DMTest, DMInit2)
{
    // construct a DM
    elecstate::DensityMatrix<double> DM(kv, paraV);
    Grid_Driver gd(0, 0, 0);
    DM.init_DMR(&gd, &ucell);
    std::cout << "dim0: " << paraV->dim0 << "    dim1:" << paraV->dim1 << std::endl;
    // construct another DM
    elecstate::DensityMatrix<double> DM1(kv, paraV);
    DM1.init_DMR(*DM.get_DMR_pointer());
    // compare
    EXPECT_EQ(DM1.get_DMR_pointer()->size_atom_pairs(), test_size * test_size);
    EXPECT_EQ(DM1.get_DMR_pointer()->get_atom_pair(2, 2).get_atom_i(), 2);
    EXPECT_EQ(DM1.get_DMR_pointer()->get_atom_pair(2, 2).get_atom_j(), 2);
    EXPECT_EQ(DM1.get_DMR_pointer()->get_atom_pair(2, 2).get_row_size(), paraV->get_row_size(2));
    EXPECT_EQ(DM1.get_DMR_pointer()->get_atom_pair(2, 2).get_col_size(), paraV->get_col_size(2));
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