//---------------------------------------------
// TEST for PW_Basis_big & PW_Basis_big_k
//---------------------------------------------
#include "../pw_basis_k.h"
#ifdef __MPI
#include "test_tool.h"
#include "module_base/parallel_global.h"
#include "mpi.h"
#endif
#include "module_base/constants.h"
#include "module_base/global_function.h"
#include "pw_test.h"

using namespace std;
TEST_F(PWTEST,test_big)
{
    cout<<"Temporary: test for pw_basis_big and pw_basis_k_big. (They should be removed in the future)"<<endl;
    ModulePW::PW_Basis_Big pwtest(device_flag, precision_flag);
    ModulePW::PW_Basis_K_Big pwktest(device_flag, precision_flag);
    ModuleBase::Matrix3 latvec;
    int nx,ny,nz;  //f*G
    double wfcecut;
    double lat0;
    bool gamma_only;
    int nks;
    ModuleBase::Vector3<double> *kvec_d;
    //--------------------------------------------------
    lat0 = 3;
    ModuleBase::Matrix3 la(1, 1, 0, 0, 1, 0, 0, 0, 1);
    latvec = la;
    wfcecut = 20;
    nks = 2;
    kvec_d = new ModuleBase::Vector3<double>[nks];
    kvec_d[0].set(0,0,0.5);
    kvec_d[1].set(0.5,0.5,0.5);
    gamma_only = false;
    int distribution_type = 1;
    bool xprime = true;
    //--------------------------------------------------
    //init //Real parameters
#ifdef __MPI
    pwtest.initmpi(nproc_in_pool, rank_in_pool, POOL_WORLD);
    pwktest.initmpi(nproc_in_pool, rank_in_pool, POOL_WORLD);
#endif
    pwtest.setbxyz(2,2,2);
    pwktest.setbxyz(2,2,2);
    pwtest.initgrids(lat0,latvec, 11, 11, 11);
    EXPECT_EQ(pwtest.nx%2, 0);
    EXPECT_EQ(pwtest.ny%2, 0);
    EXPECT_EQ(pwtest.nz%2, 0);
    pwtest.initgrids(lat0,latvec, 2*wfcecut);
    pwtest.initgrids(lat0,latvec, 3*wfcecut);
    pwktest.initgrids(lat0,latvec, pwtest.nx, pwtest.ny, pwtest.nz);
    pwtest.initparameters(gamma_only,wfcecut,distribution_type,xprime);
    pwktest.initparameters(gamma_only,wfcecut,nks,kvec_d,distribution_type, xprime);
    static_cast<ModulePW::PW_Basis>(pwtest).setuptransform();
    pwktest.setuptransform();
    EXPECT_EQ(pwtest.nx%2, 0);
    EXPECT_EQ(pwtest.ny%2, 0);
    EXPECT_EQ(pwtest.nz%2, 0);
    EXPECT_EQ(pwktest.nx%2, 0);
    EXPECT_EQ(pwktest.ny%2, 0);
    EXPECT_EQ(pwktest.nz%2, 0);

    int bsize = 0;
    pwtest.autoset_big_cell_size(bsize, 12);
    EXPECT_EQ(bsize, 4);
    pwtest.autoset_big_cell_size(bsize, 12, 4);
    EXPECT_EQ(bsize, 3);
    pwtest.autoset_big_cell_size(bsize, 14, 4);
    EXPECT_EQ(bsize, 2);


    delete[] kvec_d;
    ModulePW::PW_Basis_Big *p_pw = new ModulePW::PW_Basis_Big(device_flag, precision_flag);
    ModulePW::PW_Basis_K_Big *p_pwk = new ModulePW::PW_Basis_K_Big(device_flag, precision_flag);
    delete p_pw;
    delete p_pwk;
    fftw_cleanup();
#ifdef __ENABLE_FLOAT_FFTW
    fftwf_cleanup();
#endif
}

class TestPW_Basis_Big : public ::testing::Test
{
    public:
    ModulePW::PW_Basis_Big pwtest = ModulePW::PW_Basis_Big();
};

// Test the function with nproc = 0 (bx and by)
TEST_F(TestPW_Basis_Big, BxByTest) {
    int b_size = 0;
    int nc_size = 12;
    pwtest.autoset_big_cell_size(b_size, nc_size);
    EXPECT_EQ(b_size, 4);
}

// Test the function with nproc > 0 (bz)
TEST_F(TestPW_Basis_Big, BzTest) {
    int b_size = 0;
    int nc_size = 12;
    int nproc = 2;
    pwtest.autoset_big_cell_size(b_size, nc_size, nproc);
    EXPECT_EQ(b_size, 3);
}

// Test the function with nproc > 0 (bz) and nc_size not factored by any candidate
TEST_F(TestPW_Basis_Big, BzNoFactorTest) {
    int b_size = 0;
    int nc_size = 11;
    int nproc = 2;
    pwtest.autoset_big_cell_size(b_size, nc_size, nproc);
    EXPECT_EQ(b_size, 4);
}

// Test the function with nproc > 0 (bz) and nc_size not factored by any candidate
TEST_F(TestPW_Basis_Big, BzNoFactorNoResultTest) {
    int b_size = 0;
    int nc_size = 11;
    int nproc = 3;
    pwtest.autoset_big_cell_size(b_size, nc_size, nproc);
    EXPECT_EQ(b_size, 4);
}

// Test the function with nproc > 0 (bz) and nc_size smaller than candidates
TEST_F(TestPW_Basis_Big, BzSmallTest) {
    int b_size = 0;
    int nc_size = 2;
    int nproc = 2;
    pwtest.autoset_big_cell_size(b_size, nc_size, nproc);
    EXPECT_EQ(b_size, 2);
}

// Test the function with nproc > 0 (bz) and nc_size smaller than candidates
TEST_F(TestPW_Basis_Big, BzSmallNoResultTest) {
    int b_size = 0;
    int nc_size = 2;
    int nproc = 3;
    pwtest.autoset_big_cell_size(b_size, nc_size, nproc);
    EXPECT_EQ(b_size, 2);
}

// Test the function with nproc > 0 (bz) and nc_size not divisible by nproc
TEST_F(TestPW_Basis_Big, BzNprocTest) {
    int b_size = 0;
    int nc_size = 12;
    int nproc = 3;
    pwtest.autoset_big_cell_size(b_size, nc_size, nproc);
    EXPECT_EQ(b_size, 4);
}

// Test the function with nproc > 0 (bz) and nc_size not divisible by nproc
TEST_F(TestPW_Basis_Big, BzNprocNoResultTest) {
    int b_size = 0;
    int nc_size = 12;
    int nproc = 5;
    pwtest.autoset_big_cell_size(b_size, nc_size, nproc);
    EXPECT_EQ(b_size, 3);
}