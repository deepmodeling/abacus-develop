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
#include "pw_test.h"

using namespace std;
TEST_F(PWTEST,test4_4)
{
    cout<<"dividemthd 2, gamma_only: off, xprime: true, 2 kpoints, check fft"<<endl;
    ModulePW::PW_Basis_K pwtest(device_flag, precision_flag);
    ModuleBase::Matrix3 latvec;
    int nx,ny,nz;  //f*G
    double wfcecut;
    double lat0;
    bool gamma_only;
    ModuleBase::Vector3<double> *kvec_d;
    int nks;
    //--------------------------------------------------
    lat0 = 4;
    ModuleBase::Matrix3 la(1, 5, 0, 0, 1, 0, 0, 0, 1);
    nks = 2;
    kvec_d = new ModuleBase::Vector3<double>[nks];
    kvec_d[0].set(0,0,0);
    kvec_d[1].set(0,0.5,0.5);
    latvec = la;
    wfcecut = 40;
    gamma_only = false;
    int distribution_type = 2;
    bool xprime = true;
    //--------------------------------------------------

    //init //real parameter
#ifdef __MPI
    pwtest.initmpi(nproc_in_pool, rank_in_pool, POOL_WORLD);
#endif
    pwtest.initgrids(lat0,latvec,4*wfcecut);
    pwtest.initparameters(gamma_only,wfcecut,nks,kvec_d,distribution_type,xprime);
    pwtest.setuptransform();
    pwtest.collect_local_pw();

    const int nrxx = pwtest.nrxx;
    const int nmaxgr = pwtest.nmaxgr;
    nx = pwtest.nx;
    ny = pwtest.ny;
    nz = pwtest.nz;
    const int nplane = pwtest.nplane;

    double tpiba2 = ModuleBase::TWO_PI * ModuleBase::TWO_PI / lat0 / lat0;
    double ggecut = wfcecut / tpiba2;
    ModuleBase::Matrix3 GT,G,GGT;
    GT = latvec.Inverse();
	G  = GT.Transpose();
	GGT = G * GT;
    complex<double> *tmp = new complex<double> [nx*ny*nz];
    complex<double> * rhor = new complex<double> [nrxx];
    complex<double> * rhogr = new complex<double> [nmaxgr];
#ifdef __ENABLE_FLOAT_FFTW
    complex<float> * rhofr = new complex<float> [nrxx];
    complex<float> * rhofgr = new complex<float> [nmaxgr];
#endif
    for(int ik  = 0; ik < nks; ++ik)
    {
        int npwk = pwtest.npwk[ik];
        if(rank_in_pool == 0)
        {
            ModuleBase::Vector3<double> kk = kvec_d[ik];
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
                        // double modulus = v * (GGT * v);
                        double modulusgk = (v+kk) * (GGT * (v+kk));
                        if (modulusgk <= ggecut)
                        {
                            tmp[ix*ny*nz + iy*nz + iz]=1.0/(modulusgk+1) + ModuleBase::IMAG_UNIT / (std::abs(v.x+1) + 1);
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
        complex<double> * rhog = new complex<double> [npwk];
        complex<double> * rhogout = new complex<double> [npwk];
#ifdef __ENABLE_FLOAT_FFTW
        complex<float> * rhofg = new complex<float> [npwk];
        complex<float> * rhofgout = new complex<float> [npwk];
#endif
        for(int ig = 0 ; ig < npwk ; ++ig)
        {
            rhog[ig] = 1.0/(pwtest.getgk2(ik,ig)+1) + ModuleBase::IMAG_UNIT / (std::abs(pwtest.getgdirect(ik,ig).x+1) + 1);
            rhogr[ig] = 1.0/(pwtest.getgk2(ik,ig)+1) + ModuleBase::IMAG_UNIT / (std::abs(pwtest.getgdirect(ik,ig).x+1) + 1);
        }    
#ifdef __ENABLE_FLOAT_FFTW
        for(int ig = 0 ; ig < npwk ; ++ig)
        {
            rhofg[ig] = 1.0/(pwtest.getgk2(ik,ig)+1) + ModuleBase::IMAG_UNIT / (std::abs(pwtest.getgdirect(ik,ig).x+1) + 1);
            rhofgr[ig] = 1.0/(pwtest.getgk2(ik,ig)+1) + ModuleBase::IMAG_UNIT / (std::abs(pwtest.getgdirect(ik,ig).x+1) + 1);
        }  
#endif

        ModuleBase::GlobalFunc::ZEROS(rhor, nrxx);
        pwtest.recip2real<double>(rhog,rhor,ik, true, 1); //check out-of-place transform

        pwtest.recip2real(rhogr,rhogr,ik); //check in-place transform

#ifdef __ENABLE_FLOAT_FFTW
        ModuleBase::GlobalFunc::ZEROS(rhofr, nrxx);
        pwtest.recip2real(rhofg,rhofr,ik, true, float(1)); //check out-of-place transform

        pwtest.recip2real(rhofgr,rhofgr,ik); //check in-place transform
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

        ModuleBase::GlobalFunc::ZEROS(rhogout, npwk);
        pwtest.real2recip<double>(rhor,rhogout,ik, true, 1);

        pwtest.real2recip(rhogr,rhogr,ik);

#ifdef __ENABLE_FLOAT_FFTW
        ModuleBase::GlobalFunc::ZEROS(rhofgout, npwk);
        pwtest.real2recip(rhofr,rhofgout,ik, true, float(1));

        pwtest.real2recip(rhofgr,rhofgr,ik);
#endif

        for(int ig = 0 ; ig < npwk ; ++ig)
        {
            EXPECT_NEAR(rhog[ig].real(),rhogout[ig].real(),1e-6);
            EXPECT_NEAR(rhog[ig].imag(),rhogout[ig].imag(),1e-6);
            EXPECT_NEAR(rhog[ig].real(),rhogr[ig].real(),1e-6);
            EXPECT_NEAR(rhog[ig].imag(),rhogr[ig].imag(),1e-6);
#ifdef __ENABLE_FLOAT_FFTW
            EXPECT_NEAR(rhofg[ig].real(),rhofgout[ig].real(),1e-4);
            EXPECT_NEAR(rhofg[ig].imag(),rhofgout[ig].imag(),1e-4);
            EXPECT_NEAR(rhofg[ig].real(),rhofgr[ig].real(),1e-4);
            EXPECT_NEAR(rhofg[ig].imag(),rhofgr[ig].imag(),1e-4);
#endif
        }


        delete [] rhog;
        delete [] rhogout;
#ifdef __ENABLE_FLOAT_FFTW
        delete [] rhofg;
        delete [] rhofgout;
#endif
        //check getgcar(ik,ig)
        for(int igl = 0 ; igl < npwk; ++igl)
        {
            EXPECT_NEAR(pwtest.getgcar(ik,igl).norm2(), (pwtest.getgdirect(ik,igl) * G).norm2(), 1e-8);
            EXPECT_NEAR(pwtest.getgpluskcar(ik,igl).norm2(), ((pwtest.getgdirect(ik,igl) + kvec_d[ik]) * G).norm2(), 1e-8);
        }

        //check igl2ig
        for(int igl = 0; igl < npwk ; ++igl)
        {        
            const int isz = pwtest.getigl2isz(ik,igl);
            for(int ig = 0 ; ig < pwtest.npw; ++ig)
            {
                if(isz == pwtest.ig2isz[ig]){
                    EXPECT_EQ(ig,pwtest.getigl2ig(ik,igl));}
            }
        }

        //check getig2ixyz_k
        pwtest.get_ig2ixyz_k();
        for(int igl = 0; igl < npwk ; ++igl)
        { 
            EXPECT_GE(pwtest.ig2ixyz_k_[igl + ik * pwtest.npwk_max], 0);
        }

    }
    delete []tmp; 
    delete [] rhor;
    delete[] kvec_d;
    delete[] rhogr;
    fftw_cleanup();
#ifdef __ENABLE_FLOAT_FFTW
    delete[] rhofr;
    delete[] rhofgr;
    fftwf_cleanup();
#endif
}