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
        kv->isk.resize(nks);
        kv->isk[0] = 0;
        kv->isk[1] = 1;
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
    // Call the function to be tested
    ModuleIO::nscf_band(is, out_band_dir, nks, nband, fermie, ekb, *kv, Pkpoints);

    // Check the output file
    std::ifstream ifs(out_band_dir);
    std::string str((std::istreambuf_iterator<char>(ifs)),std::istreambuf_iterator<char>());
    ASSERT_TRUE(ifs.is_open());
    EXPECT_THAT(str, testing::HasSubstr("1       -27.2114       -13.6057              0"));
    ifs.close();
}

#include "mpi.h"
int main(int argc, char **argv)
{
        MPI_Init(&argc, &argv);
        testing::InitGoogleTest(&argc, argv);

        int nproc;
        int my_rank;

        MPI_Comm_size(MPI_COMM_WORLD,&nproc);
        MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

        int result = RUN_ALL_TESTS();

        MPI_Finalize();

        return result;
}