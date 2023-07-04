//---------------------------------------------
// TEST for FFT
//---------------------------------------------
#include "../pw_basis.h"
#ifdef __MPI
#include "test_tool.h"
#include "module_base/parallel_global.h"
#include "mpi.h"
#endif
#include "module_base/constants.h"
#include "module_base/global_function.h"
#include "pw_test.h"

using namespace std;
TEST_F(PWTEST,test5_4_2)
{
    cout<<"dividemthd 1, gamma_only: off, xprime: false, full_pw: true, full_pw_dim: 2, check fft between complex and complex"<<endl;
    ModulePW::PW_Basis pwtest(device_flag, precision_flag);
    ModuleBase::Matrix3 latvec;
    int nx,ny,nz;  //f*G
    double wfcecut;
    double lat0;
    bool gamma_only;
    //--------------------------------------------------
    lat0 = 2;
    ModuleBase::Matrix3 la(1, 1, 0, 0, 4, 0, 0, 0, 1);
    latvec = la;
    wfcecut = 40;
    gamma_only = false;
    int distribution_type = 1;
    bool xprime = false;
    bool full_pw = true;
    int full_pw_dim = 2;
    //--------------------------------------------------
#ifdef __MPI
    pwtest.initmpi(nproc_in_pool, rank_in_pool, POOL_WORLD);
#endif
    pwtest.setfullpw(full_pw, full_pw_dim);
    pwtest.initgrids(lat0,latvec,wfcecut);
    pwtest.initparameters(gamma_only,wfcecut,distribution_type,xprime);
    pwtest.setuptransform();
    pwtest.collect_local_pw();

    const int npw = pwtest.npw;
    const int nrxx = pwtest.nrxx;
    const int nmaxgr = pwtest.nmaxgr;
    nx = pwtest.nx;
    ny = pwtest.ny;
    nz = pwtest.nz;
    int nplane = pwtest.nplane;

    double tpiba2 = ModuleBase::TWO_PI * ModuleBase::TWO_PI / lat0 / lat0;
    double ggecut = wfcecut / tpiba2;
    ModuleBase::Matrix3 GT,G,GGT;
    GT = latvec.Inverse();
	G  = GT.Transpose();
	GGT = G * GT;
    complex<double> *tmp = new complex<double> [nx*ny*nz];
    if(rank_in_pool == 0)
    {
        for(int ix = 0 ; ix < nx ; ++ix)
        {
            for(int iy = 0 ; iy < ny ; ++iy)
            {
                for(int iz = 0 ; iz < nz ; ++iz)
                {
                    tmp[ix*ny*nz + iy*nz + iz]=0.0;
                    double vx = ix + int(nx/2) - nx + 1;
                    double vy = iy + int(ny/2) - ny + 1;
                    double vz = iz + int(nz/2) - nz + 1;
                    ModuleBase::Vector3<double> v(vx,vy,vz);
                    double modulus = v * (GGT * v);
                    if (modulus <= ggecut || full_pw)
                    {
                        tmp[ix*ny*nz + iy*nz + iz]=1.0/(modulus+1) + ModuleBase::IMAG_UNIT / (std::abs(v.x+1) + 1);
                    }
                }
            }   
        }
        fftw_plan pp = fftw_plan_dft_3d(nx,ny,nz,(fftw_complex *) tmp, (fftw_complex *) tmp, FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_execute(pp);    
        fftw_destroy_plan(pp); 
        
        ModuleBase::Vector3<double> delta_g((double(int(nx/2))-1.)/nx, (double(int(ny/2))-1.)/ny, (double(int(nz/2))-1.)/nz); 
        for(int ixy = 0 ; ixy < nx * ny ; ++ixy)
        {
            for(int iz = 0 ; iz < nz ; ++iz)
            {
                int ix = ixy / ny;
                int iy = ixy % ny;
                ModuleBase::Vector3<double> real_r(ix, iy, iz);
                double phase_im = -delta_g * real_r;
                complex<double> phase(0,ModuleBase::TWO_PI * phase_im);
                tmp[ixy * nz + iz] *= exp(phase);
            }
        }
    }
#ifdef __MPI
    MPI_Bcast(tmp,2*nx*ny*nz,MPI_DOUBLE,0,POOL_WORLD);
#endif
    
    complex<double> * rhog = new complex<double> [npw];
    complex<double> * rhogr = new complex<double> [nmaxgr];
    complex<double> * rhogout = new complex<double> [npw];
    for(int ig = 0 ; ig < npw ; ++ig)
    {
        rhog[ig] = 1.0/(pwtest.gg[ig]+1) + ModuleBase::IMAG_UNIT / (std::abs(pwtest.gdirect[ig].x+1) + 1);
        rhogr[ig] = 1.0/(pwtest.gg[ig]+1) + ModuleBase::IMAG_UNIT / (std::abs(pwtest.gdirect[ig].x+1) + 1);
    }    
    complex<double> * rhor = new complex<double> [nrxx];
#ifdef __ENABLE_FLOAT_FFTW
    complex<float> * rhofg = new complex<float> [npw];
    complex<float> * rhofgr = new complex<float> [nmaxgr];
    complex<float> * rhofgout = new complex<float> [npw];
    for(int ig = 0 ; ig < npw ; ++ig)
    {
        rhofg[ig] = 1.0/(pwtest.gg[ig]+1) + ModuleBase::IMAG_UNIT / (std::abs(pwtest.gdirect[ig].x+1) + 1);
        rhofgr[ig] = 1.0/(pwtest.gg[ig]+1) + ModuleBase::IMAG_UNIT / (std::abs(pwtest.gdirect[ig].x+1) + 1);
    }    
    complex<float> * rhofr = new complex<float> [nrxx];
#endif
    
    pwtest.recip2real(rhog,rhor);//check out-of-place transform

    pwtest.recip2real(rhogr,rhogr);//check in-place transform

#ifdef __ENABLE_FLOAT_FFTW
    pwtest.recip2real(rhofg,rhofr);//check out-of-place transform

    pwtest.recip2real(rhofgr,rhofgr);//check in-place transform
#endif

    int startiz = pwtest.startz_current;
    for(int ixy = 0 ; ixy < nx * ny ; ++ixy)
    {
        for(int iz = 0 ; iz < nplane ; ++iz)
        {
            EXPECT_NEAR(tmp[ixy * nz + startiz + iz].real(),rhor[ixy*nplane+iz].real(),1e-6);
            EXPECT_NEAR(tmp[ixy * nz + startiz + iz].imag(),rhor[ixy*nplane+iz].imag(),1e-6);
            EXPECT_NEAR(tmp[ixy * nz + startiz + iz].real(),rhogr[ixy*nplane+iz].real(),1e-6);
            EXPECT_NEAR(tmp[ixy * nz + startiz + iz].imag(),rhogr[ixy*nplane+iz].imag(),1e-6);
#ifdef __ENABLE_FLOAT_FFTW
            EXPECT_NEAR(tmp[ixy * nz + startiz + iz].real(),rhofr[ixy*nplane+iz].real(),1e-4);
            EXPECT_NEAR(tmp[ixy * nz + startiz + iz].imag(),rhofr[ixy*nplane+iz].imag(),1e-4);
            EXPECT_NEAR(tmp[ixy * nz + startiz + iz].real(),rhofgr[ixy*nplane+iz].real(),1e-4);
            EXPECT_NEAR(tmp[ixy * nz + startiz + iz].imag(),rhofgr[ixy*nplane+iz].imag(),1e-4);
#endif
        }
    }

    
    
    pwtest.real2recip(rhor,rhogout);//check out-of-place transform

    pwtest.real2recip(rhogr,rhogr);//check in-place transform

#ifdef __ENABLE_FLOAT_FFTW
    pwtest.real2recip(rhofr,rhofgout);//check out-of-place transform

    pwtest.real2recip(rhofgr,rhofgr);//check in-place transform
#endif

    for(int ig = 0 ; ig < npw ; ++ig)
    {
        EXPECT_NEAR(rhog[ig].real(),rhogout[ig].real(),1e-6);
        EXPECT_NEAR(rhog[ig].imag(),rhogout[ig].imag(),1e-6);
        EXPECT_NEAR(rhogr[ig].real(),rhogout[ig].real(),1e-6);
        EXPECT_NEAR(rhogr[ig].imag(),rhogout[ig].imag(),1e-6);
#ifdef __ENABLE_FLOAT_FFTW
        EXPECT_NEAR(rhofg[ig].real(),rhofgout[ig].real(),1e-4);
        EXPECT_NEAR(rhofg[ig].imag(),rhofgout[ig].imag(),1e-4);
        EXPECT_NEAR(rhofgr[ig].real(),rhofgout[ig].real(),1e-4);
        EXPECT_NEAR(rhofgr[ig].imag(),rhofgout[ig].imag(),1e-4);
#endif
    }
  
    delete [] rhog;
    delete [] rhogout;
    delete [] rhor;
    delete [] tmp; 
    delete [] rhogr;

    fftw_cleanup();
#ifdef __ENABLE_FLOAT_FFTW
    delete [] rhofg;
    delete [] rhofgout;
    delete [] rhofr;
    delete [] rhofgr;
    fftwf_cleanup();
#endif
}