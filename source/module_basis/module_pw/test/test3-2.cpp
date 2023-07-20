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
TEST_F(PWTEST,test3_2)
{
    cout<<"dividemthd 1, gamma_only: off, xprime: true, check fft between complex and complex, reset ggecut to latecut"<<endl;
    ModulePW::PW_Basis pwtest(device_flag, precision_flag);
    ModuleBase::Matrix3 latvec;
    int nx,ny,nz;  //f*G
    double wfcecut;
    double lat0;
    bool gamma_only;
    //--------------------------------------------------
    lat0 = 3;
    ModuleBase::Matrix3 la(1, 1, 0, 0, 1, 0, 0, 0, 1);
    latvec = la;
    wfcecut = 20;
    gamma_only = false;
    int distribution_type = 1;
    bool xprime = true;
    //--------------------------------------------------
    //init //Real parameters
#ifdef __MPI
    pwtest.initmpi(nproc_in_pool, rank_in_pool, POOL_WORLD);
#endif
    pwtest.initgrids(lat0,latvec,wfcecut);
    pwtest.initparameters(gamma_only, 2 * wfcecut, distribution_type, xprime); // pwtest will reset 2*wfcecut to wfcecut
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
                    double vx = ix -  int(nx/2);
                    double vy = iy -  int(ny/2);
                    double vz = iz -  int(nz/2);
                    ModuleBase::Vector3<double> v(vx,vy,vz);
                    double modulus = v * (GGT * v);
                    if (modulus <= ggecut)
                    {
                        tmp[ix*ny*nz + iy*nz + iz]=1.0/(modulus+1) + ModuleBase::IMAG_UNIT / (std::abs(v.x+1) + 1);
                    }
                }
            }   
        }
        fftw_plan pp = fftw_plan_dft_3d(nx,ny,nz,(fftw_complex *) tmp, (fftw_complex *) tmp, FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_execute(pp);    
        fftw_destroy_plan(pp); 
        
        ModuleBase::Vector3<double> delta_g(double(int(nx/2))/nx, double(int(ny/2))/ny, double(int(nz/2))/nz); 
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
    ModuleBase::GlobalFunc::ZEROS(rhor, nrxx);
    
    pwtest.recip2real(rhog,rhor);//check out-of-place transform

    pwtest.recip2real(rhogr,rhogr);//check in-place transform

    int startiz = pwtest.startz_current;
    for(int ixy = 0 ; ixy < nx * ny ; ++ixy)
    {
        for(int iz = 0 ; iz < nplane ; ++iz)
        {
            EXPECT_NEAR(tmp[ixy * nz + startiz + iz].real(),rhor[ixy*nplane+iz].real(),1e-6);
            EXPECT_NEAR(tmp[ixy * nz + startiz + iz].imag(),rhor[ixy*nplane+iz].imag(),1e-6);
            EXPECT_NEAR(tmp[ixy * nz + startiz + iz].real(),rhogr[ixy*nplane+iz].real(),1e-6);
            EXPECT_NEAR(tmp[ixy * nz + startiz + iz].imag(),rhogr[ixy*nplane+iz].imag(),1e-6);
        }
    }

    
    ModuleBase::GlobalFunc::ZEROS(rhogout, npw);
    pwtest.real2recip<double>(rhor,rhogout,true, 1);//check out-of-place transform

    pwtest.real2recip(rhogr,rhogr);//check in-place transform

    for(int ig = 0 ; ig < npw ; ++ig)
    {
        EXPECT_NEAR(rhog[ig].real(),rhogout[ig].real(),1e-6);
        EXPECT_NEAR(rhog[ig].imag(),rhogout[ig].imag(),1e-6);
        EXPECT_NEAR(rhogr[ig].real(),rhogout[ig].real(),1e-6);
        EXPECT_NEAR(rhogr[ig].imag(),rhogout[ig].imag(),1e-6);
    }
  
    delete [] rhog;
    delete [] rhogout;
    delete [] rhor;
    delete [] tmp; 
    delete [] rhogr;

    fftw_cleanup();
#ifdef __ENABLE_FLOAT_FFTW
    fftwf_cleanup();
#endif
}