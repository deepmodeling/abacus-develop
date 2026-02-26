#include <cerrno>
#include <fstream>
#include <iomanip>
#include <sys/stat.h>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#define private public
#include "source_estate/module_dm/density_matrix.h"
#include "source_lcao/module_hcontainer/hcontainer.h"
#include "source_lcao/module_hcontainer/read_hcontainer.h"
#include "source_lcao/setup_dm.h"
#include "source_cell/klist.h"
#undef private
#include "source_io/module_dm/write_dmr.h"

/************************************************
 *  unit test of init_dm_from_file (nspin=1 & nspin=2)
 *
 *  Uses write_dmr_csr to generate CSR files in the
 *  current format, then reads them back via Read_HContainer
 *  to verify the round-trip.
 ***********************************************/

// Small test system: 2 atoms, 4 orbitals each => nlocal=8
int test_size = 2;
int test_nw = 4;

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
        ucell.atoms[0].taud.resize(ucell.nat);
        ucell.itia2iat.create(ucell.ntype, ucell.nat);
        for (int iat = 0; iat < ucell.nat; iat++)
        {
            ucell.iat2it[iat] = 0;
            ucell.iat2ia[iat] = iat;
            ucell.atoms[0].tau[iat] = ModuleBase::Vector3<double>(0.0, 0.0, 0.0);
            ucell.atoms[0].taud[iat] = ModuleBase::Vector3<double>(0.0, 0.0, 0.0);
            ucell.itia2iat(0, iat) = iat;
        }
        ucell.atoms[0].na = test_size;
        ucell.atoms[0].nw = test_nw;
        ucell.atoms[0].iw2l.resize(test_nw, 0);
        ucell.atoms[0].iw2m.resize(test_nw, 0);
        ucell.atoms[0].iw2n.resize(test_nw, 0);
        ucell.atoms[0].label = "Si";
        ucell.latName = "fcc";
        ucell.lat0 = 10.0;
        ucell.latvec.e11 = 1.0; ucell.latvec.e12 = 0.0; ucell.latvec.e13 = 0.0;
        ucell.latvec.e21 = 0.0; ucell.latvec.e22 = 1.0; ucell.latvec.e23 = 0.0;
        ucell.latvec.e31 = 0.0; ucell.latvec.e32 = 0.0; ucell.latvec.e33 = 1.0;
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

    /// Build an HContainer with diagonal values at R=(0,0,0) and write it via write_dmr_csr
    void write_test_csr(const std::string& filename, double scale, int ispin, int nspin)
    {
        hamilt::HContainer<double> hc(paraV);
        for (int i = 0; i < ucell.nat; i++)
        {
            for (int j = 0; j < ucell.nat; j++)
            {
                hamilt::AtomPair<double> ap(i, j, 0, 0, 0, paraV);
                hc.insert_pair(ap);
            }
        }
        hc.allocate(nullptr, true);

        // Fill diagonal elements with scale-dependent values
        for (int i = 0; i < ucell.nat; i++)
        {
            auto* ap = hc.find_pair(i, i);
            if (ap)
            {
                int nw = ucell.atoms[0].nw;
                for (int k = 0; k < nw; k++)
                {
                    ap->get_pointer()[k * nw + k] = scale * (k + 1) * 0.1;
                }
            }
        }

        std::string fname = filename;
        ModuleIO::write_dmr_csr(fname, &ucell, 8, &hc, 0, ispin, nspin);
    }

    /// Create DensityMatrix with given nspin and initialize DMR from an HContainer template
    elecstate::DensityMatrix<double, double>* create_dm(int nspin)
    {
        K_Vectors kv;
        int nks = (nspin == 2) ? 2 : 1;
        kv.set_nks(nks * (nspin == 2 ? 2 : 1));
        kv.kvec_d.resize(kv.get_nks());

        int nspin_dm = (nspin == 2) ? 2 : 1;
        auto* dm = new elecstate::DensityMatrix<double, double>(
            paraV, nspin_dm, kv.kvec_d, kv.get_nks() / nspin_dm);

        // Create a template HContainer and init DMR from it
        hamilt::HContainer<double> tmp_HR(paraV);
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
    mkdir("./test_dm_dir", 0755);
    write_test_csr("./test_dm_dir/dmrs1_nao.csr", 1.0, 0, 1);

    auto* dm = create_dm(1);
    ASSERT_EQ(dm->_DMR.size(), 1);

    hamilt::HContainer<double>* dmr0 = dm->get_DMR_vector()[0];
    hamilt::Read_HContainer<double> reader(dmr0, "./test_dm_dir/dmrs1_nao.csr", nlocal, &ucell);
    reader.read();

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
}

TEST_F(InitDMFileTest, Nspin2_ReadTwoFiles)
{
    mkdir("./test_dm_dir", 0755);
    write_test_csr("./test_dm_dir/dmrs1_nao.csr", 1.0, 0, 2);  // spin-up
    write_test_csr("./test_dm_dir/dmrs2_nao.csr", 0.5, 1, 2);  // spin-down

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
            EXPECT_NEAR(v1 / v0, 0.5, 1e-6);
            break;
        }
    }
    EXPECT_TRUE(values_differ);

    delete dm;
}

TEST_F(InitDMFileTest, Nspin2_DMRVectorSize)
{
    auto* dm = create_dm(2);
    EXPECT_EQ(dm->_DMR.size(), 2);
    EXPECT_NE(dm->_DMR[0], nullptr);
    EXPECT_NE(dm->_DMR[1], nullptr);
    delete dm;
}

TEST_F(InitDMFileTest, Nspin1_DMRVectorSize)
{
    auto* dm = create_dm(1);
    EXPECT_EQ(dm->_DMR.size(), 1);
    EXPECT_NE(dm->_DMR[0], nullptr);
    delete dm;
}

/************************************************
 *  unit test of init_hr_from_file (init_chg=hr)
 *
 *  Tests HR CSR file round-trip via Read_HContainer,
 *  and nspin=2 dual-buffer read into two halves of
 *  a single HContainer (same logic as init_chg_hr).
 ***********************************************/

TEST_F(InitDMFileTest, HR_Nspin1_ReadSingleFile)
{
    // Write an HR CSR file (same format as DM CSR)
    mkdir("./test_hr_dir", 0755);
    write_test_csr("./test_hr_dir/hrs1_nao.csr", 2.0, 0, 1);

    // Create an HContainer to read into
    hamilt::HContainer<double> hR(paraV);
    for (int i = 0; i < ucell.nat; i++)
    {
        for (int j = 0; j < ucell.nat; j++)
        {
            hamilt::AtomPair<double> ap(i, j, 0, 0, 0, paraV);
            hR.insert_pair(ap);
        }
    }
    hR.allocate(nullptr, true);

    // Read HR from file (same as init_hr_from_file does internally)
    hR.set_zero();
    hamilt::Read_HContainer<double> reader(&hR, "./test_hr_dir/hrs1_nao.csr", nlocal, &ucell);
    reader.read();

    // Verify data was loaded
    EXPECT_GT(hR.size_atom_pairs(), 0);
    auto* ap = hR.find_pair(0, 0);
    ASSERT_NE(ap, nullptr);

    // Check diagonal has expected values (scale=2.0, value = 2.0*(k+1)*0.1)
    int nw = ucell.atoms[0].nw;
    for (int k = 0; k < nw; k++)
    {
        double expected = 2.0 * (k + 1) * 0.1;
        EXPECT_NEAR(ap->get_pointer()[k * nw + k], expected, 1e-6);
    }
}

TEST_F(InitDMFileTest, HR_Nspin2_ReadTwoFiles)
{
    // Write two HR CSR files with different scale factors
    mkdir("./test_hr_dir", 0755);
    write_test_csr("./test_hr_dir/hrs1_nao.csr", 1.0, 0, 2);  // spin-up
    write_test_csr("./test_hr_dir/hrs2_nao.csr", 3.0, 1, 2);  // spin-down

    // Create two independent HContainers for spin-up and spin-down
    // (mirrors init_chg_hr reading two separate HR files)
    auto create_hcontainer = [&]() {
        hamilt::HContainer<double>* hR = new hamilt::HContainer<double>(paraV);
        for (int i = 0; i < ucell.nat; i++)
        {
            for (int j = 0; j < ucell.nat; j++)
            {
                hamilt::AtomPair<double> ap(i, j, 0, 0, 0, paraV);
                hR->insert_pair(ap);
            }
        }
        hR->allocate(nullptr, true);
        return hR;
    };

    // Read spin-up
    auto* hR_up = create_hcontainer();
    hR_up->set_zero();
    hamilt::Read_HContainer<double> reader_up(hR_up, "./test_hr_dir/hrs1_nao.csr", nlocal, &ucell);
    reader_up.read();

    // Read spin-down
    auto* hR_down = create_hcontainer();
    hR_down->set_zero();
    hamilt::Read_HContainer<double> reader_down(hR_down, "./test_hr_dir/hrs2_nao.csr", nlocal, &ucell);
    reader_down.read();

    // Verify both have data
    EXPECT_GT(hR_up->size_atom_pairs(), 0);
    EXPECT_GT(hR_down->size_atom_pairs(), 0);

    // Verify spin-up values (scale=1.0)
    auto* ap_up = hR_up->find_pair(0, 0);
    ASSERT_NE(ap_up, nullptr);
    int nw = ucell.atoms[0].nw;
    for (int k = 0; k < nw; k++)
    {
        double expected_up = 1.0 * (k + 1) * 0.1;
        EXPECT_NEAR(ap_up->get_pointer()[k * nw + k], expected_up, 1e-6);
    }

    // Verify spin-down values (scale=3.0)
    auto* ap_down = hR_down->find_pair(0, 0);
    ASSERT_NE(ap_down, nullptr);
    for (int k = 0; k < nw; k++)
    {
        double expected_down = 3.0 * (k + 1) * 0.1;
        EXPECT_NEAR(ap_down->get_pointer()[k * nw + k], expected_down, 1e-6);
    }

    // Verify the two are independent (different values)
    double val_up = ap_up->get_pointer()[0];
    double val_down = ap_down->get_pointer()[0];
    EXPECT_GT(std::abs(val_up), 1e-15);
    EXPECT_NEAR(val_down / val_up, 3.0, 1e-6);

    delete hR_up;
    delete hR_down;
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
