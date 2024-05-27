#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "module_io/nscf_band.h"
#include "module_cell/parallel_kpoints.h"
#include "module_cell/klist.h"

Parallel_Kpoints::Parallel_Kpoints()
{
    nks_pool = nullptr;
    startk_pool = nullptr;
    whichpool = nullptr;
}

Parallel_Kpoints::~Parallel_Kpoints()
{
    delete[] nks_pool;
    delete[] startk_pool;
    delete[] whichpool;
}

K_Vectors::K_Vectors()
{
}

K_Vectors::~K_Vectors()
{
}

int K_Vectors::get_nks() const
{
    return this->nks;
}

int K_Vectors::get_nkstot() const
{
    return this->nkstot;
}

int K_Vectors::get_nkstot_ibz() const
{
    return this->nkstot_ibz;
}

int K_Vectors::get_nkstot_full() const
{
    return this->nkstot_full;
}

void K_Vectors::set_nks(int value) {
    this->nks = value;
}

void K_Vectors::set_nkstot(int value) {
    this->nkstot = value;
}

void K_Vectors::set_nkstot_ibz(int value) {
    this->nkstot_ibz = value;
}

void K_Vectors::set_nkstot_full(int value) {
    this->nkstot_full = value;
}

/************************************************
 *  unit test of nscf_band
 ***********************************************/

/**
 * - Tested Functions:
 *   - nscf_band()
 *     - output band structure in nscf calculation
 */

class BandTest : public ::testing::Test
{
protected:
    void SetUp() override {
        // Set up test data
        is = 0;
        out_band_dir = "test_band.dat";
        nks = 2;
        nband = 3;
        fermie = 0.0;
        ekb.create(nks, nband);
	    ekb(0,0) = -2.0;
	    ekb(0,1) = -1.0;
	    ekb(0,2) =  0.0;
	    ekb(1,0) =  1.0;
	    ekb(1,1) =  2.0;
	    ekb(1,2) =  3.0;
        kv = new K_Vectors;
        // specify the kpoints
        kv->kvec_c.resize(nks);
        kv->kvec_c[0] = ModuleBase::Vector3<double>(0.0, 0.0, 0.0);
        kv->kvec_c[1] = ModuleBase::Vector3<double>(1.0, 0.0, 0.0);
        kv->isk.resize(nks);
        kv->isk[0] = 0;
        kv->isk[1] = 1;
        kv->kl_segids.resize(nks);
        kv->kl_segids[0] = 0;
        kv->kl_segids[1] = 0;
        Pkpoints = new Parallel_Kpoints;
    }

    void TearDown() override {
        // Clean up test data
        delete kv;
        delete Pkpoints;
        std::remove(out_band_dir.c_str());
    }

    // Test data
    int is;
    std::string out_band_dir;
    int nks;
    int nband;
    double fermie;
    ModuleBase::matrix ekb;
    K_Vectors* kv;
    Parallel_Kpoints* Pkpoints;
};

TEST_F(BandTest, nscf_band)
{
    Pkpoints->nks_pool = new int[1];
    Pkpoints->nks_pool[0] = nks;
    Pkpoints->nkstot_np = nks;
    Pkpoints->nks_np = nks;
    // Call the function to be tested
    ModuleIO::nscf_band(is, out_band_dir, nband, fermie, 8, ekb, *kv, Pkpoints);

    // Check the output file
    std::ifstream ifs(out_band_dir);
    std::string str((std::istreambuf_iterator<char>(ifs)),std::istreambuf_iterator<char>());
    ASSERT_TRUE(ifs.is_open());
    EXPECT_THAT(str, testing::HasSubstr("   1 0.00000000 -27.21139600 -13.60569800 0.00000000"));
    ifs.close();
}
