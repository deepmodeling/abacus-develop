//---------------------------------------------
// TEST for FFT
//---------------------------------------------
#include "../pw_basis.h"
#ifdef __MPI
#include "test_tool.h"
#include "mpi.h"
#endif
#include "../../module_base/constants.h"
#include "../../module_base/global_function.h"
#include <iostream>
#include <iomanip>
#include "../../module_base/timer.h"

using namespace std;
int main(int argc,char **argv)
{
    int nproc, myrank;
    int nproc_in_pool, kpar, mypool, rank_in_pool;
    ModulePW::PW_Basis pwtest;
    ModuleBase::Matrix3 latvec;
    int nx,ny,nz;  //f*G
    double wfcecut;
    double lat0;
    bool gamma_only;
    //--------------------------------------------------
    lat0 = 2;
    ModuleBase::Matrix3 la(1, 1, 0, 0, 2, 0, 0, 0, 2);
    latvec = la;
    wfcecut = 10;
    kpar = 1;
    gamma_only = false;
    //--------------------------------------------------
    
    //setup mpi
#ifdef __MPI
    setupmpi(argc,argv,nproc, myrank);
    divide_pools(nproc, myrank, nproc_in_pool, kpar, mypool, rank_in_pool);
#else
    nproc = nproc_in_pool = kpar = 1;
    myrank = mypool = rank_in_pool = 0;
#endif
    //cout<<nproc<<" d "<<myrank<<" d "<<nproc_in_pool<<" "<<kpar<<" "<<mypool<<" "<<rank_in_pool<<endl;
    ModuleBase::timer::start();
    
    //init
    pwtest.initgrids(lat0,latvec,wfcecut);
    //pwtest.initgrids(lat0,latvec,5,7,7);
    pwtest.initparameters(gamma_only,wfcecut,nproc_in_pool,rank_in_pool,1);
    pwtest.setuptransform();
    pwtest.collect_local_pw();

    int npw = pwtest.npw;
    int nrxx = pwtest.nrxx;
    nx = pwtest.nx;
    ny = pwtest.ny;
    nz = pwtest.nz;
    int nplane = pwtest.nplane;
    int nxyz = nx * ny * nz;
    if(myrank == 0) cout<<"FFT: "<<nx<<" "<<ny<<" "<<nz<<endl;
    double tpiba2 = ModuleBase::TWO_PI * ModuleBase::TWO_PI / lat0 / lat0;
    double ggecut = wfcecut / tpiba2;
    ModuleBase::Matrix3 GT,G,GGT;
    GT = latvec.Inverse();
	G  = GT.Transpose();
	GGT = G * GT;
    complex<double> *tmp = NULL;
    if(myrank == 0)
    {
        tmp = new complex<double> [nx*ny*nz];
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
                        tmp[ix*ny*nz + iy*nz + iz]=1.0/(modulus+1) + ModuleBase::IMAG_UNIT / (abs(v.x+1) + 1);
                        //tmp[ix*ny*nz + iy*nz + iz] = 1.0;
                    }
                }
            }   
        }
        fftw_plan pp = fftw_plan_dft_3d(nx,ny,nz,(fftw_complex *) tmp, (fftw_complex *) tmp, FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_execute(pp);    
        fftw_destroy_plan(pp); 
        
        //output
        cout << "reference\n";
        ModuleBase::Vector3<double> delta_g(double(int(nx/2))/nx, double(int(ny/2))/ny, double(int(ny/2))/nz); 
        for(int ixy = 0 ; ixy < nx * ny ; ixy+=5)
        {
            for(int iz = 0 ; iz < nz ; ++iz)
            {
                int ix = ixy / ny;
                int iy = ixy % ny;
                ModuleBase::Vector3<double> real_r(ix, iy, iz);
                double phase_im = -delta_g * real_r;
                complex<double> phase(0,ModuleBase::TWO_PI * phase_im);
                tmp[ixy * nz + iz] /= nxyz;
                tmp[ixy * nz + iz] *= exp(phase);
                cout<<setprecision(5)<<setiosflags(ios::left)<<setw(30)<<tmp[ixy * nz + iz];
            }
        }
        cout<<endl;
    }
    
    complex<float> * rhog = new complex<float> [npw];
    for(int ig = 0 ; ig < npw ; ++ig)
    {
        rhog[ig] = complex<float>(1.0/(pwtest.gg[ig]+1) + ModuleBase::IMAG_UNIT / (abs(pwtest.gdirect[ig].x+1) + 1));
        //rhog[ig] = 1.0/(pwtest.gg[ig]+1);
        //rhog[ig] = 1.0;
    }    
    complex<float> * rhor = new complex<float> [nrxx];
    pwtest.recip2real(rhog,rhor);
    if(myrank == 0)     cout << "new pw module\n";
#ifdef __MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    for(int ixy = 0 ; ixy < nx * ny ; ixy+=5)
    {
        for(int ip = 0 ; ip < nproc ; ++ip)
        {
        if (myrank == ip)
        {
            for(int iz = 0 ; iz < nplane ; ++iz)
            {
                cout<<setprecision(5)<<setiosflags(ios::left)<<setw(30)<<rhor[ixy*nplane+iz];
            }
        }
#ifdef __MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif          
        } 
    }
    
    if(myrank == 0)             cout<<endl<<endl;

    if(myrank == nproc - 1)
    {
        cout<<"before transform: "<<endl;
        for(int ig = 0 ; ig < npw ; ++ig)
        {
            cout<<rhog[ig]<<" ";
        }
        cout<<endl;
    }
    
    pwtest.real2recip(rhor,rhog);
    
    if(myrank == nproc - 1)
    {
        cout<<"after transform:"<<endl;
        for(int ig = 0 ; ig < npw ; ++ig)
        {
            cout<<rhog[ig]<<" ";
        }
        cout<<endl;
    }
    
    if(rank_in_pool==0) ModuleBase::timer::finish(GlobalV::ofs_running, true);

#ifdef __MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif     
    delete [] rhog;
    delete [] rhor;
    if(tmp!=NULL) delete []tmp; 
    fftw_cleanup();
#ifdef __MIX_PRECISION
    fftwf_cleanup();
#endif
#ifdef __MPI
    finishmpi();
#endif
    return 0;
}