#include <fstream>
#include <iomanip>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#define private public
#include "source_estate/module_dm/density_matrix.h"
#include "source_lcao/module_hcontainer/hcontainer.h"
#include "source_lcao/module_hcontainer/read_hcontainer.h"
#include "source_lcao/setup_dm.h"
#include "source_cell/klist.h"
#undef private

/************************************************
 *  unit test of init_dm_from_file (nspin=1 & nspin=2)
 ***********************************************/

// Small test system: 2 atoms, 13 orbitals each => nlocal=26
int test_size = 2;
int test_nw = 13;

class InitDMFileTest : public testing::Test
{
  protected:
    Parallel_Orbitals* paraV;
    UnitCell ucell;
    int nlocal;

    void SetUp() override
    {
        nlocal = test_size * test_nw;

        // set up a unitcell
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
        ucell.atoms[0].iw2l.resize(test_nw, 0);
        ucell.atoms[0].iw2m.resize(test_nw, 0);
        ucell.atoms[0].iw2n.resize(test_nw, 0);
        ucell.set_iat2iwt(1);

        // set up parallel orbitals (serial mode)
        paraV = new Parallel_Orbitals();
        paraV->set_serial(nlocal, nlocal);
        paraV->set_atomic_trace(ucell.get_iat2iwt(), ucell.nat, nlocal);
    }

    void TearDown() override
    {
        delete paraV;
        delete[] ucell.atoms;
    }

    /// Write a minimal CSR file with a diagonal matrix at R=(0,0,0)
    void write_test_csr(const std::string& filename, double scale)
    {
        std::ofstream ofs(filename);
        ofs << "IONIC_STEP: 1" << std::endl;
        ofs << "Matrix Dimension of DM(R): " << nlocal << std::endl;
        ofs << "Matrix number of DM(R): 1" << std::endl;

        // R coordinate header: rx ry rz nnz
        int nnz = nlocal; // diagonal
        ofs << "0 0 0 " << nnz << std::endl;

        // values line
        for (int i = 0; i < nlocal; i++)
        {
            if (i > 0) ofs << " ";
            ofs << std::scientific << std::setprecision(8) << scale * (i + 1) * 0.01;
        }
        ofs << std::endl;

        // column indices line
        for (int i = 0; i < nlocal; i++)
        {
            if (i > 0) ofs << " ";
            ofs << i;
        }
        ofs << std::endl;

        // row pointers line (CSR format: each row has exactly 1 element on diagonal)
        for (int i = 0; i <= nlocal; i++)
        {
            if (i > 0) ofs << " ";
            ofs << i;
        }
        ofs << std::endl;

        ofs.close();
    }

    /// Create DensityMatrix with given nspin and initialize DMR from an HContainer template
    elecstate::DensityMatrix<double, double>* create_dm(int nspin)
    {
        K_Vectors kv;
        int nks = (nspin == 2) ? 2 : 1; // gamma_only: nk=1 per spin
        kv.set_nks(nks * (nspin == 2 ? 2 : 1));
        kv.kvec_d.resize(kv.get_nks());

        int nspin_dm = (nspin == 2) ? 2 : 1;
        auto* dm = new elecstate::DensityMatrix<double, double>(
            paraV, nspin_dm, kv.kvec_d, kv.get_nks() / nspin_dm);

        // Create a template HContainer and init DMR from it
        hamilt::HContainer<double> tmp_HR(paraV);
        // Add atom pairs for all atom-atom combinations at R=(0,0,0)
        for (int i = 0; i < ucell.nat; i++)
        {
            for (int j = 0; j < ucell.nat; j++)
            {
                hamilt::AtomPair<double> ap(i, j, 0, 0, 0, paraV);
                tmp_HR.insert_pair(ap);
            }
        }
        tmp_HR.allocate(nullptr, true);
        dm->init_DMR(tmp_HR);
        return dm;
    }
};

TEST_F(InitDMFileTest, Nspin1_ReadSingleFile)
{
    // Create test directory and CSR file
    system("mkdir -p ./test_dm_dir");
    write_test_csr("./test_dm_dir/dmrs1_nao.csr", 1.0);

    // Create DM with nspin=1
    auto* dm = create_dm(1);
    ASSERT_EQ(dm->_DMR.size(), 1);

    // Read from file using Read_HContainer (same as init_dm_from_file does)
    hamilt::HContainer<double>* dmr0 = dm->get_DMR_vector()[0];
    hamilt::Read_HContainer<double> reader(dmr0, "./test_dm_dir/dmrs1_nao.csr", nlocal, &ucell);
    reader.read();

    // Verify DMR[0] has data
    EXPECT_GT(dmr0->size_atom_pairs(), 0);

    // Check diagonal element (0,0) at R=(0,0,0) is non-zero
    auto* ap = dmr0->find_pair(0, 0);
    ASSERT_NE(ap, nullptr);
    bool has_nonzero = false;
    for (int i = 0; i < ap->get_size(); i++)
    {
        if (std::abs(ap->get_pointer()[i]) > 1e-15)
        {
            has_nonzero = true;
            break;
        }
    }
    EXPECT_TRUE(has_nonzero);

    delete dm;
    system("rm -rf ./test_dm_dir");
}

TEST_F(InitDMFileTest, Nspin2_ReadTwoFiles)
{
    // Create test directory and two CSR files with different scale factors
    system("mkdir -p ./test_dm_dir");
    write_test_csr("./test_dm_dir/dmrs1_nao.csr", 1.0);  // spin-up
    write_test_csr("./test_dm_dir/dmrs2_nao.csr", 0.5);  // spin-down

    // Create DM with nspin=2
    auto* dm = create_dm(2);
    ASSERT_EQ(dm->_DMR.size(), 2);

    // Read spin-up
    hamilt::HContainer<double>* dmr0 = dm->get_DMR_vector()[0];
    hamilt::Read_HContainer<double> reader0(dmr0, "./test_dm_dir/dmrs1_nao.csr", nlocal, &ucell);
    reader0.read();

    // Read spin-down
    hamilt::HContainer<double>* dmr1 = dm->get_DMR_vector()[1];
    hamilt::Read_HContainer<double> reader1(dmr1, "./test_dm_dir/dmrs2_nao.csr", nlocal, &ucell);
    reader1.read();

    // Verify both DMR components have data
    EXPECT_GT(dmr0->size_atom_pairs(), 0);
    EXPECT_GT(dmr1->size_atom_pairs(), 0);

    // Verify spin-up and spin-down have different values (scale 1.0 vs 0.5)
    auto* ap0 = dmr0->find_pair(0, 0);
    auto* ap1 = dmr1->find_pair(0, 0);
    ASSERT_NE(ap0, nullptr);
    ASSERT_NE(ap1, nullptr);

    bool values_differ = false;
    int check_size = std::min(ap0->get_size(), ap1->get_size());
    for (int i = 0; i < check_size; i++)
    {
        double v0 = ap0->get_pointer()[i];
        double v1 = ap1->get_pointer()[i];
        if (std::abs(v0) > 1e-15 && std::abs(v0 - v1) > 1e-15)
        {
            values_differ = true;
            // spin-down should be ~half of spin-up
            EXPECT_NEAR(v1 / v0, 0.5, 1e-6);
            break;
        }
    }
    EXPECT_TRUE(values_differ);

    delete dm;
    system("rm -rf ./test_dm_dir");
}

TEST_F(InitDMFileTest, Nspin2_DMRVectorSize)
{
    // Verify that nspin=2 creates exactly 2 DMR components
    auto* dm = create_dm(2);
    EXPECT_EQ(dm->_DMR.size(), 2);
    EXPECT_NE(dm->_DMR[0], nullptr);
    EXPECT_NE(dm->_DMR[1], nullptr);
    delete dm;
}

TEST_F(InitDMFileTest, Nspin1_DMRVectorSize)
{
    // Verify that nspin=1 creates exactly 1 DMR component
    auto* dm = create_dm(1);
    EXPECT_EQ(dm->_DMR.size(), 1);
    EXPECT_NE(dm->_DMR[0], nullptr);
    delete dm;
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
