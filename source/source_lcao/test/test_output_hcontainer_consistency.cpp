#ifdef __MPI
#include <mpi.h>
#endif
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include <fstream>
#include <sstream>
#include <sys/stat.h>

#define private public
#include "source_lcao/module_hcontainer/hcontainer.h"
#include "source_lcao/module_hcontainer/output_hcontainer.h"
#include "source_lcao/module_hcontainer/read_hcontainer.h"
#undef private

#include "source_cell/unitcell.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_io/module_dm/write_dmr.h"

/************************************************
 *  Unit test for Output_HContainer write-read
 *  round-trip consistency.
 *
 *  Uses write_dmr_csr (which internally calls
 *  Output_HContainer::write) to produce a full
 *  CSR file, then reads it back via Read_HContainer
 *  (which uses csrFileReader) and compares values.
 ***********************************************/

class OutputHContainerTest : public testing::Test
{
protected:
    Parallel_Orbitals* paraV;
    UnitCell ucell;
    int nlocal;
    int test_size = 2;
    int test_nw = 4;

    void SetUp() override
    {
        nlocal = test_size * test_nw;

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

        paraV = new Parallel_Orbitals();
        paraV->set_serial(nlocal, nlocal);
        paraV->set_atomic_trace(ucell.get_iat2iwt(), ucell.nat, nlocal);
    }

    void TearDown() override
    {
        delete paraV;
        delete[] ucell.atoms;
    }

    /// Build an HContainer with diagonal values at R=(0,0,0)
    hamilt::HContainer<double>* create_hcontainer(double scale)
    {
        auto* hc = new hamilt::HContainer<double>(paraV);
        for (int i = 0; i < ucell.nat; i++)
        {
            for (int j = 0; j < ucell.nat; j++)
            {
                hamilt::AtomPair<double> ap(i, j, 0, 0, 0, paraV);
                hc->insert_pair(ap);
            }
        }
        hc->allocate(nullptr, true);

        for (int i = 0; i < ucell.nat; i++)
        {
            auto* ap = hc->find_pair(i, i);
            if (ap)
            {
                int nw = ucell.atoms[0].nw;
                for (int k = 0; k < nw; k++)
                {
                    ap->get_pointer()[k * nw + k] = scale * (k + 1) * 0.1;
                }
            }
        }
        return hc;
    }

    /// Write a full CSR file via write_dmr_csr (uses Output_HContainer internally)
    void write_full_csr(const std::string& filename, hamilt::HContainer<double>* hc,
                        int ispin, int nspin, int precision = 8)
    {
        std::string fname = filename;
        ModuleIO::write_dmr_csr(fname, &ucell, precision, hc, 0, ispin, nspin);
    }
};

// ---- Test 1: write-read round-trip, verify element-wise consistency ----
TEST_F(OutputHContainerTest, WriteReadConsistency)
{
    mkdir("./test_ohc_dir", 0755);

    auto* hc_write = create_hcontainer(2.5);
    write_full_csr("./test_ohc_dir/hrs1_nao.csr", hc_write, 0, 1);

    // Read back
    auto* hc_read = create_hcontainer(0.0);
    hc_read->set_zero();
    hamilt::Read_HContainer<double> reader(hc_read, "./test_ohc_dir/hrs1_nao.csr", nlocal, &ucell);
    reader.read();

    // Compare every diagonal element of every atom pair
    int nw = ucell.atoms[0].nw;
    for (int i = 0; i < ucell.nat; i++)
    {
        auto* ap_w = hc_write->find_pair(i, i);
        auto* ap_r = hc_read->find_pair(i, i);
        ASSERT_NE(ap_w, nullptr);
        ASSERT_NE(ap_r, nullptr);

        for (int k = 0; k < nw; k++)
        {
            EXPECT_NEAR(ap_w->get_pointer()[k * nw + k],
                        ap_r->get_pointer()[k * nw + k], 1e-8)
                << "Mismatch at atom " << i << " orbital " << k;
        }
    }

    // Off-diagonal atom pairs should be zero in both
    for (int i = 0; i < ucell.nat; i++)
    {
        for (int j = 0; j < ucell.nat; j++)
        {
            if (i == j) continue;
            auto* ap_w = hc_write->find_pair(i, j);
            auto* ap_r = hc_read->find_pair(i, j);
            if (ap_w && ap_r)
            {
                for (int k = 0; k < nw * nw; k++)
                {
                    EXPECT_NEAR(ap_w->get_pointer()[k],
                                ap_r->get_pointer()[k], 1e-10);
                }
            }
        }
    }

    delete hc_write;
    delete hc_read;
    std::remove("./test_ohc_dir/hrs1_nao.csr");
    rmdir("./test_ohc_dir");
}

// ---- Test 2: sparse threshold filters out tiny values ----
TEST_F(OutputHContainerTest, SparseThresholdFiltering)
{
    mkdir("./test_ohc_dir", 0755);

    // scale = 1e-12 => all values < default sparse_threshold (1e-10)
    auto* hc = create_hcontainer(1e-12);
    write_full_csr("./test_ohc_dir/hrs1_nao.csr", hc, 0, 1);

    // Verify the CSR file has no actual data values after the header.
    // When all values are below sparse_threshold, Output_HContainer
    // writes no R-block data (no "Rx Ry Rz nonZero" line appears).
    std::ifstream ifs("./test_ohc_dir/hrs1_nao.csr");
    ASSERT_TRUE(ifs.is_open());
    std::string line;
    std::string file_content;
    while (std::getline(ifs, line))
    {
        file_content += line + "\n";
    }
    ifs.close();

    // After the fix, nonZero=0 R-blocks are also written to keep
    // nR consistent. Verify the file contains "0 0 0 0" (R=(0,0,0) with 0 nonzero).
    EXPECT_NE(file_content.find("0 0 0 0"), std::string::npos)
        << "Expected R=(0,0,0) block with 0 nonzero elements";

    delete hc;
    std::remove("./test_ohc_dir/hrs1_nao.csr");
    rmdir("./test_ohc_dir");
}

// ---- Test 3: precision parameter affects round-trip accuracy ----
TEST_F(OutputHContainerTest, PrecisionParameter)
{
    mkdir("./test_ohc_dir", 0755);

    auto* hc = create_hcontainer(1.23456789);

    // Write with low precision (4 digits)
    write_full_csr("./test_ohc_dir/hrs1_nao.csr", hc, 0, 1, 4);

    auto* hc_read = create_hcontainer(0.0);
    hc_read->set_zero();
    hamilt::Read_HContainer<double> reader(hc_read, "./test_ohc_dir/hrs1_nao.csr", nlocal, &ucell);
    reader.read();

    int nw = ucell.atoms[0].nw;
    auto* ap_r = hc_read->find_pair(0, 0);
    ASSERT_NE(ap_r, nullptr);

    for (int k = 0; k < nw; k++)
    {
        double expected = 1.23456789 * (k + 1) * 0.1;
        // 4-digit precision => ~1e-4 tolerance
        EXPECT_NEAR(ap_r->get_pointer()[k * nw + k], expected, 5e-4)
            << "Low-precision round-trip failed at orbital " << k;
    }

    delete hc;
    delete hc_read;
    std::remove("./test_ohc_dir/hrs1_nao.csr");
    rmdir("./test_ohc_dir");
}

// ---- Test 4: nspin=2 two-file round-trip ----
TEST_F(OutputHContainerTest, Nspin2TwoFileConsistency)
{
    mkdir("./test_ohc_dir", 0755);

    auto* hc_up = create_hcontainer(1.0);
    auto* hc_down = create_hcontainer(3.0);
    write_full_csr("./test_ohc_dir/hrs1_nao.csr", hc_up, 0, 2);
    write_full_csr("./test_ohc_dir/hrs2_nao.csr", hc_down, 1, 2);

    // Read back spin-up
    auto* hc_read_up = create_hcontainer(0.0);
    hc_read_up->set_zero();
    hamilt::Read_HContainer<double> reader_up(hc_read_up, "./test_ohc_dir/hrs1_nao.csr", nlocal, &ucell);
    reader_up.read();

    // Read back spin-down
    auto* hc_read_down = create_hcontainer(0.0);
    hc_read_down->set_zero();
    hamilt::Read_HContainer<double> reader_down(hc_read_down, "./test_ohc_dir/hrs2_nao.csr", nlocal, &ucell);
    reader_down.read();

    int nw = ucell.atoms[0].nw;
    for (int k = 0; k < nw; k++)
    {
        double exp_up = 1.0 * (k + 1) * 0.1;
        double exp_down = 3.0 * (k + 1) * 0.1;
        EXPECT_NEAR(hc_read_up->find_pair(0, 0)->get_pointer()[k * nw + k], exp_up, 1e-8);
        EXPECT_NEAR(hc_read_down->find_pair(0, 0)->get_pointer()[k * nw + k], exp_down, 1e-8);
    }

    // Verify the two are independent
    EXPECT_GT(std::abs(hc_read_up->find_pair(0, 0)->get_pointer()[0]
                     - hc_read_down->find_pair(0, 0)->get_pointer()[0]), 1e-6);

    delete hc_up;
    delete hc_down;
    delete hc_read_up;
    delete hc_read_down;
    std::remove("./test_ohc_dir/hrs1_nao.csr");
    std::remove("./test_ohc_dir/hrs2_nao.csr");
    rmdir("./test_ohc_dir");
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
