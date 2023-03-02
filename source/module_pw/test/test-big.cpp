//---------------------------------------------
// TEST for PW_Basis_big & PW_Basis_big_k
//---------------------------------------------
#include "../pw_basis_k.h"
#ifdef __MPI
#include "test_tool.h"
#include "../../module_base/parallel_global.h"
#include "mpi.h"
#endif
#include "../../module_base/constants.h"
#include "../../module_base/global_function.h"
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
    pwtest.setuptransform();
    pwktest.setuptransform();
    EXPECT_EQ(pwtest.nx%2, 0);
    EXPECT_EQ(pwtest.ny%2, 0);
    EXPECT_EQ(pwtest.nz%2, 0);
    EXPECT_EQ(pwktest.nx%2, 0);
    EXPECT_EQ(pwktest.ny%2, 0);
    EXPECT_EQ(pwktest.nz%2, 0);

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