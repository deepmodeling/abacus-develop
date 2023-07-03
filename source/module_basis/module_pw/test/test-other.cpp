//---------------------------------------------
// TEST for FFT
//---------------------------------------------
#include "../pw_basis_k.h"
#ifdef __MPI
#include "test_tool.h"
#include "module_base/parallel_global.h"
#include "mpi.h"
#endif
#include "module_base/constants.h"
#include "module_base/global_function.h"
#include "module_psi/kernels/types.h"
#include "pw_test.h"
#include "gmock/gmock.h"

using namespace std;
TEST_F(PWTEST,test_other)
{
    cout<<"Test other codes"<<endl;
    ModulePW::PW_Basis pwtest(device_flag, precision_flag);
    ModulePW::PW_Basis_K pwktest(device_flag, precision_flag);
    ModuleBase::Matrix3 latvec(0.2, 0, 0, 0, 1, 0, 0, 0, 1);
#ifdef __MPI
    pwtest.initmpi(nproc_in_pool, rank_in_pool, POOL_WORLD);
    pwktest.initmpi(nproc_in_pool, rank_in_pool, POOL_WORLD);
#endif
    pwtest.initgrids(3, latvec, 2, 10, 10);
    pwtest.initparameters(false, 20, 3); //distribute_type = 3
    testing::internal::CaptureStdout();
    EXPECT_EXIT(pwtest.setuptransform(), ::testing::ExitedWithCode(0), "");
    string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output,testing::HasSubstr("NOTICE"));

    int nks = 2;
    ModuleBase::Vector3<double> *kvec_d = new ModuleBase::Vector3<double>[nks];
    kvec_d[0].set(0,0,0.5);
    kvec_d[1].set(0.5,0.5,0.5);
    pwktest.set_precision("double");
    pwktest.initgrids(2, latvec, 4,4,4);
    pwktest.initparameters(true, 20, nks, kvec_d);
    pwktest.setuptransform();
    pwktest.collect_local_pw();
#ifdef __ENABLE_FLOAT_FFTW
    pwktest.set_precision("single");
#endif
    pwktest.initparameters(true, 8, nks, kvec_d);
    pwktest.setuptransform();
    pwktest.collect_local_pw();
    const int nrxx = pwktest.nrxx;
    complex<double> * rhor1 = new complex<double> [nrxx];
    complex<double> * rhor2 = new complex<double> [nrxx];
#ifdef __ENABLE_FLOAT_FFTW
    complex<float> * rhofr1 = new complex<float> [nrxx];
    complex<float> * rhofr2 = new complex<float> [nrxx];
#endif
    const psi::DEVICE_CPU * ctx;
    for(int ik = 0; ik < nks; ++ik)
    {
        const int npwk = pwktest.npwk[ik];
        complex<double> * rhog1 = new complex<double> [npwk];
        complex<double> * rhog2 = new complex<double> [npwk];
#ifdef __ENABLE_FLOAT_FFTW
        complex<float> * rhofg1 = new complex<float> [npwk];
        complex<float> * rhofg2 = new complex<float> [npwk];
#endif
        for(int ig = 0 ; ig < npwk ; ++ig)
        {
            rhog1[ig] = 1.0/(pwktest.getgk2(ik,ig)+1) + ModuleBase::IMAG_UNIT / (std::abs(pwktest.getgdirect(ik,ig).x+1) + 1);
            rhog2[ig] = 1.0/(pwktest.getgk2(ik,ig)+1) + ModuleBase::IMAG_UNIT / (std::abs(pwktest.getgdirect(ik,ig).x+1) + 1);
        }    
#ifdef __ENABLE_FLOAT_FFTW
        for(int ig = 0 ; ig < npwk ; ++ig)
        {
            rhofg1[ig] = 1.0/(pwktest.getgk2(ik,ig)+1) + ModuleBase::IMAG_UNIT / (std::abs(pwktest.getgdirect(ik,ig).x+1) + 1);
            rhofg2[ig] = 1.0/(pwktest.getgk2(ik,ig)+1) + ModuleBase::IMAG_UNIT / (std::abs(pwktest.getgdirect(ik,ig).x+1) + 1);
        }  
#endif

        pwktest.recip_to_real(ctx, rhog1, rhor1, ik);
        pwktest.recip2real(rhog2, rhor2, ik);
        for(int ir = 0 ; ir < nrxx; ++ir)
        {
            EXPECT_NEAR(std::abs(rhor1[ir]),std::abs(rhor2[ir]),1e-8);
        }
        pwktest.real_to_recip(ctx, rhor1, rhog1, ik);
        pwktest.real2recip(rhor2, rhog2, ik);
        for(int ig = 0 ; ig < npwk; ++ig)
        {
            EXPECT_NEAR(std::abs(rhog1[ig]),std::abs(rhog2[ig]),1e-8);
        }
#ifdef __ENABLE_FLOAT_FFTW
        pwktest.recip_to_real(ctx, rhofg1, rhofr1, ik);
        pwktest.recip2real(rhofg2, rhofr2, ik);
        for(int ir = 0 ; ir < nrxx; ++ir)
        {
            EXPECT_NEAR(std::abs(rhofr1[ir]),std::abs(rhofr2[ir]),1e-6);
        }
        pwktest.real_to_recip(ctx, rhofr1, rhofg1, ik);
        pwktest.real2recip(rhofr2, rhofg2, ik);
        for(int ig = 0 ; ig < npwk; ++ig)
        {
            EXPECT_NEAR(std::abs(rhofg1[ig]),std::abs(rhofg2[ig]),1e-6);
        }
#endif



        delete [] rhog1;
        delete [] rhog2;
#ifdef __ENABLE_FLOAT_FFTW
        delete [] rhofg1;
        delete [] rhofg2;
#endif
    }
    delete [] rhor1;
    delete [] rhor2;
#ifdef __ENABLE_FLOAT_FFTW
    delete [] rhofr1;
    delete [] rhofr2;
#endif


    double* d_kvec_c = pwktest.get_kvec_c_data<double>();
    double* d_gcar = pwktest.get_gcar_data<double>();
    double* d_gk2 = pwktest.get_gk2_data<double>();
#ifdef __ENABLE_FLOAT_FFTW
    float* s_kvec_c = pwktest.get_kvec_c_data<float>();
    float* s_gcar = pwktest.get_gcar_data<float>();
    float* s_gk2 = pwktest.get_gk2_data<float>();
#endif

    
    delete[] kvec_d;
    ModulePW::PW_Basis *p_pw = new ModulePW::PW_Basis(device_flag, precision_flag);
    ModulePW::PW_Basis_K *p_pwk = new ModulePW::PW_Basis_K(device_flag, precision_flag);
    delete p_pw;
    delete p_pwk;
    fftw_cleanup();
#ifdef __ENABLE_FLOAT_FFTW
    fftwf_cleanup();
#endif
}