#include <chrono>

#include <array>
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "source_base/basis_flattened_indexing.h"
#include "source_base/element_basis_index.h"
#include "source_estate/module_dm/density_matrix.h"
#include "source_hamilt/module_hcontainer/hcontainer.h"
#include "source_cell/klist.h"
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
        ucell.iat2it.resize(ucell.nat);
        ucell.iat2ia.resize(ucell.nat);
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
        // set paraV
        init_parav();
    }

    void TearDown() override
    {
        delete paraV;
        delete[] ucell.atoms;
    }

#ifdef __MPI
    void init_parav()
    {
        int nb = 2;
        int global_row = test_size * test_nw;
        int global_col = test_size * test_nw;
        std::ofstream ofs_running;
        paraV = new Parallel_Orbitals();
        paraV->para_2d().init(global_row, global_col, nb, MPI_COMM_WORLD, false);
        paraV->set_atomic_trace(ucell.get_iat2iwt(), test_size, global_row);
    }
#else
    void init_parav()
    {
    }
#endif
};

// test for construct DMR from GridD and UnitCell
TEST_F(DMTest, DMInit1)
{
    // initalize a kvectors
    K_Vectors* kv = nullptr;
    int nspin = 1;
    int nks = 2; // since nspin = 1
    kv = new K_Vectors;
    kv->set_nks(nks);
    kv->kvec_d.resize(nks);
    kv->kvec_d[1].x = 0.5;
    // construct DM
    std::cout << "dim0: " << paraV->para_2d().nrproc() << "    dim1:" << paraV->para_2d().ncproc() << std::endl;
    std::cout << "nrow: " << paraV->para_2d().nrloc() << "    ncol:" << paraV->para_2d().ncloc() << std::endl;
    module_dm::DensityMatrix<double, double> DM(paraV, nspin, kv->kvec_d, nks);
    // initialize this->dmr
    Grid_Driver gd(0,0);
    DM.init_dmr(&gd, &ucell);
    // compare
    EXPECT_EQ(DM.get_dmr_ptr(1)->size_atom_pairs(), test_size * test_size);
    EXPECT_EQ(DM.get_dmr_ptr(1)->get_atom_pair(2, 2).get_atom_i(), 2);
    EXPECT_EQ(DM.get_dmr_ptr(1)->get_atom_pair(2, 2).get_atom_j(), 2);
    EXPECT_EQ(DM.get_dmr_ptr(1)->get_atom_pair(2, 2).get_row_size(), paraV->get_nrow_atom(2));
    EXPECT_EQ(DM.get_dmr_ptr(1)->get_atom_pair(2, 2).get_col_size(), paraV->get_ncol_atom(2));
    delete kv;
}

TEST_F(DMTest, DMInitExplicitLayout)
{
    ModuleBase::Element_Basis_Index::Range range(1);
    range[0].resize(1);
    range[0][0].N = test_nw;
    range[0][0].M = 1;
    const ModuleBase::Element_Basis_Index::IndexLNM hierarchy
        = ModuleBase::Element_Basis_Index::construct_index(range);
    const std::vector<int> atom_types(test_size, 0);
    const BasisFlattenedIndexing basis_index(hierarchy, atom_types, 1);

    std::vector<ModuleBase::Vector3<double>> kvec_d(2);
    kvec_d[1].x = 0.5;
    module_dm::DensityMatrix<double, double> DM(
        basis_index, paraV->para_2d(), 1, kvec_d, 1);
    Grid_Driver gd(0, 0);
    DM.init_dmr(&gd, &ucell);

    EXPECT_EQ(DM.get_dmr_ptr(1)->size_atom_pairs(), test_size * test_size);
    EXPECT_EQ(DM.get_dmr_ptr(1)->get_atom_pair(2, 2).get_row_size(), paraV->get_nrow_atom(2));
    EXPECT_EQ(DM.get_dmr_ptr(1)->get_atom_pair(2, 2).get_col_size(), paraV->get_ncol_atom(2));
}

// test for construct DMR from RA and UnitCell
TEST_F(DMTest, DMInit2)
{
    // initalize a kvectors
    K_Vectors* kv = nullptr;
    int nspin = 1;
    int nks = 2; // since nspin = 1
    kv = new K_Vectors;
    kv->set_nks(nks);
    kv->kvec_d.resize(nks);
    kv->kvec_d[1].x = 0.5;
    // construct DM
    std::cout << "dim0: " << paraV->para_2d().nrproc() << "    dim1:" << paraV->para_2d().ncproc() << std::endl;
    std::cout << "nrow: " << paraV->para_2d().nrloc() << "    ncol:" << paraV->para_2d().ncloc() << std::endl;
    module_dm::DensityMatrix<double, double> DM(paraV, nspin, kv->kvec_d, nks);
    // initialize Record_adj using Grid_Driver
    Grid_Driver gd(0,0);
    Record_adj ra;
    ra.na_each.resize(ucell.nat);
    ra.info_offset.resize(ucell.nat);
    // count adjacent atoms for each atom first to lay out the flat info
    std::vector<AdjacentAtomInfo> all_adjs(ucell.nat);
    int total = 0;
    for (int iat1 = 0; iat1 < ucell.nat; iat1++)
    {
        auto tau1 = ucell.get_tau(iat1);
        int T1, I1;
        ucell.iat2iait(iat1, &I1, &T1);
        gd.Find_atom(ucell, tau1, T1, I1, &all_adjs[iat1]);
        ra.na_each[iat1] = all_adjs[iat1].adj_num + 1;
        ra.info_offset[iat1] = total;
        total += ra.na_each[iat1];
    }
    ra.info.resize(total);
    for (int iat1 = 0; iat1 < ucell.nat; iat1++)
    {
        const AdjacentAtomInfo& adjs = all_adjs[iat1];
        for (int ad = 0; ad < ra.na_each[iat1]; ++ad)
        {
            const int T2 = adjs.ntype[ad];
            const int I2 = adjs.natom[ad];
            const ModuleBase::Vector3<int>& R_index = adjs.box[ad];
            std::array<int, 5>& rec = ra.info[ra.info_offset[iat1] + ad];
            rec[0] = R_index.x;
            rec[1] = R_index.y;
            rec[2] = R_index.z;
            rec[3] = T2;
            rec[4] = I2;
        }
    }
    DM.init_dmr(ra, &ucell);
    // compare
    EXPECT_EQ(DM.get_dmr_ptr(1)->size_atom_pairs(), test_size * test_size);
    EXPECT_EQ(DM.get_dmr_ptr(1)->get_atom_pair(2, 2).get_atom_i(), 2);
    EXPECT_EQ(DM.get_dmr_ptr(1)->get_atom_pair(2, 2).get_atom_j(), 2);
    EXPECT_EQ(DM.get_dmr_ptr(1)->get_atom_pair(2, 2).get_row_size(), paraV->get_nrow_atom(2));
    EXPECT_EQ(DM.get_dmr_ptr(1)->get_atom_pair(2, 2).get_col_size(), paraV->get_ncol_atom(2));
    // release memory
    delete kv;
}

// test for construct DMR from another HContainer<double>
TEST_F(DMTest, DMInit3)
{
    // initalize a kvectors
    K_Vectors* kv = nullptr;
    int nspin = 2;
    int nks = 4; // since nspin = 2
    kv = new K_Vectors;
    kv->set_nks(nks);
    kv->kvec_d.resize(nks);
    kv->kvec_d[1].x = 0.5;
    kv->kvec_d[3].x = 0.5;
    // construct a DM
    module_dm::DensityMatrix<std::complex<double>, double> DM(paraV, nspin, kv->kvec_d, kv->get_nks() / nspin);
    Grid_Driver gd(0, 0);
    DM.init_dmr(&gd, &ucell);
    std::cout << "dim0: " << paraV->para_2d().nrproc() << "    dim1:" << paraV->para_2d().ncproc() << std::endl;
    // construct another DM
    module_dm::DensityMatrix<std::complex<double>, double> DM1(paraV, nspin, kv->kvec_d, kv->get_nks() / nspin);
    DM1.init_dmr(*DM.get_dmr_ptr(1));
    // compare
    EXPECT_EQ(DM1.get_dmr_ptr(2)->size_atom_pairs(), test_size * test_size);
    EXPECT_EQ(DM1.get_dmr_ptr(2)->get_atom_pair(2, 2).get_atom_i(), 2);
    EXPECT_EQ(DM1.get_dmr_ptr(1)->get_atom_pair(2, 2).get_atom_j(), 2);
    EXPECT_EQ(DM1.get_dmr_ptr(1)->get_atom_pair(2, 2).get_row_size(), paraV->get_nrow_atom(2));
    EXPECT_EQ(DM1.get_dmr_ptr(2)->get_atom_pair(2, 2).get_col_size(), paraV->get_ncol_atom(2));
    //
    delete kv;
}

// test for construct DMR from another HContainer <std::complex<double>>
TEST_F(DMTest, DMInit4)
{
    // initalize a kvectors
    K_Vectors* kv = nullptr;
    int nspin = 2;
    int nks = 4; // since nspin = 2
    kv = new K_Vectors;
    kv->set_nks(nks);
    kv->kvec_d.resize(nks);
    kv->kvec_d[1].x = 0.5;
    kv->kvec_d[3].x = 0.5;
    // construct a new HContainer
    Grid_Driver gd(0, 0);
    hamilt::HContainer<std::complex<double>>* tmp_DMR;
    tmp_DMR = new hamilt::HContainer<std::complex<double>>(paraV);
    // set up a HContainer
    for (int iat1 = 0; iat1 < ucell.nat; iat1++)
    {
        auto tau1 = ucell.get_tau(iat1);
        int T1, I1;
        ucell.iat2iait(iat1, &I1, &T1);
        AdjacentAtomInfo adjs;
        gd.Find_atom(ucell, tau1, T1, I1, &adjs);
        // std::cout << "adjs.adj_num: " <<adjs.adj_num << std::endl;
        for (int ad = 0; ad < adjs.adj_num + 1; ++ad)
        {
            const int T2 = adjs.ntype[ad];
            const int I2 = adjs.natom[ad];
            int iat2 = ucell.itia2iat(T2, I2);
            if (paraV->is_invalid_atom_pair(iat1, iat2))
            {
                continue;
            }
            ModuleBase::Vector3<int>& R_index = adjs.box[ad];
            // std::cout << "R_index: " << R_index.x << " " << R_index.y << " " << R_index.z << std::endl;
            hamilt::AtomPair<std::complex<double>> tmp_ap(iat1, iat2, R_index.x, R_index.y, R_index.z, paraV);
            tmp_DMR->insert_pair(tmp_ap);
        }
    }
    // construct a DM from this HContainer
    module_dm::DensityMatrix<std::complex<double>, double> DM(paraV, nspin, kv->kvec_d, kv->get_nks() / nspin);
    DM.init_dmr(*tmp_DMR);
    std::cout << "dim0: " << paraV->para_2d().nrproc() << "    dim1:" << paraV->para_2d().ncproc() << std::endl;
    // compare
    EXPECT_EQ(DM.get_dmr_ptr(2)->size_atom_pairs(), test_size * test_size);
    EXPECT_EQ(DM.get_dmr_ptr(2)->get_atom_pair(2, 2).get_atom_i(), 2);
    EXPECT_EQ(DM.get_dmr_ptr(1)->get_atom_pair(2, 2).get_atom_j(), 2);
    EXPECT_EQ(DM.get_dmr_ptr(1)->get_atom_pair(2, 2).get_row_size(), paraV->get_nrow_atom(2));
    EXPECT_EQ(DM.get_dmr_ptr(2)->get_atom_pair(2, 2).get_col_size(), paraV->get_ncol_atom(2));
    //
    delete kv;
}

TEST_F(DMTest, InitDMRFromOppositeValueType)
{
    K_Vectors kv;
    kv.set_nks(1);
    kv.kvec_d.resize(1);

    hamilt::HContainer<std::complex<double>> source(paraV);
    const ModuleBase::Vector3<int> R_index(1, 0, 0);
    source.insert_pair(2, 3, R_index);
    source.allocate(nullptr, true);

    module_dm::DensityMatrix<std::complex<double>, double> dm(paraV, 1, kv.kvec_d, 1);
    dm.init_dmr(source);

    const hamilt::HContainer<double>* result = dm.get_dmr_ptr(1);
    ASSERT_NE(result, nullptr);
    ASSERT_EQ(result->size_atom_pairs(), 1);
    EXPECT_EQ(result->get_atom_pair(0).get_atom_i(), 2);
    EXPECT_EQ(result->get_atom_pair(0).get_atom_j(), 3);
    EXPECT_EQ(result->get_atom_pair(0).find_R(R_index), 0);
    EXPECT_NE(result->find_matrix(2, 3, R_index), nullptr);
}

// test for save_dmr
TEST_F(DMTest, saveDMR)
{
    // initalize a kvectors
    K_Vectors* kv = nullptr;
    int nspin = 2;
    int nks = 4; // since nspin = 2
    kv = new K_Vectors;
    kv->set_nks(nks);
    kv->kvec_d.resize(nks);
    kv->kvec_d[1].x = 0.5;
    kv->kvec_d[3].x = 0.5;
    // construct a DM
    module_dm::DensityMatrix<std::complex<double>, double> DM(paraV, nspin, kv->kvec_d, kv->get_nks() / nspin);
    Grid_Driver gd(0, 0);
    DM.init_dmr(&gd, &ucell);
    // construct another DM
    module_dm::DensityMatrix<std::complex<double>, double> DM_test(paraV, nspin, kv->kvec_d, kv->get_nks() / nspin);
    DM_test.init_dmr(*DM.get_dmr_ptr(1));
    DM_test.save_dmr();
    EXPECT_EQ(DM_test.get_dmr_ptr(1)->get_nnr(), DM.get_dmr_ptr(1)->get_nnr());
    EXPECT_EQ(DM_test.get_dmr_ptr(1)->get_nnr(), DM_test.get_dmr_save()[0].size());
    // add a new AtomPair, act as a relaxation
    hamilt::AtomPair<double> tmp_ap(9, 9, 1, 0, 0, paraV);
    DM_test.get_dmr_ptr(1)->insert_pair(tmp_ap);
    DM_test.get_dmr_ptr(1)->allocate();
    // update DMR_save
    DM_test.save_dmr();
    EXPECT_EQ(DM_test.get_dmr_ptr(1)->get_nnr(), DM_test.get_dmr_save()[0].size());
    // delete 
    delete kv;   
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
