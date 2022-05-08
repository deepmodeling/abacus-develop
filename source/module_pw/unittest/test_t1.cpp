//---------------------------------------------
// TEST for FFT
//---------------------------------------------
#include "../pw_basis.h"
#include "../../src_parallel/parallel_global.h"
#ifdef __MPI
#include "test_tool.h"
#include "mpi.h"
#endif
#include "../../module_base/constants.h"
#include "../../module_base/global_function.h"
#include <iostream>
#include <iomanip>
#include "time.h"
#include <gperftools/profiler.h>
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
    lat0 = 20;
    ModuleBase::Matrix3 la(1, 0, 0, 0, 1, 0, 0, 0, 1);
    latvec = la;
    wfcecut = 160;
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
    // MPI_Barrier(MPI_COMM_WORLD);
    // if(myrank==0) ProfilerStart("test0.prof");
    // if(myrank==1) ProfilerStart("test1.prof");
    
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
    if(myrank == 0) cout<<"FFT: "<<nx<<" "<<ny<<" "<<nz<<endl;
    ModuleBase::Matrix3 GT,G,GGT;
    GT = latvec.Inverse();
	G  = GT.Transpose();
	GGT = G * GT;
    complex<double> *tmp = NULL;
   
    complex<double> * rhog = new complex<double> [npw];
    for(int ig = 0 ; ig < npw ; ++ig)
    {
        // rhog[ig] = 1.0/(pwtest.gg[ig]+1) + ModuleBase::IMAG_UNIT / (abs(pwtest.gdirect[ig].x+1) + 1);
        rhog[ig] = 1.0/(pwtest.gg[ig]+1);
        //rhog[ig] = 1.0;
    }    

    complex<double> * rhor = new complex<double> [nrxx];

    clock_t fftstart, fftend;
#ifdef __MPI
        MPI_Barrier(POOL_WORLD);
#endif
    fftstart = clock();
    
    for (int i = 0; i < 1000; ++i)
    {
        pwtest.recip2real(rhog,rhor);
    }

#ifdef __MPI
        MPI_Barrier(POOL_WORLD);
#endif
    fftend = clock();
    double fftduration = (double)(fftend - fftstart)/CLOCKS_PER_SEC;

    cout<<"\n";
    cout<<"spend "<<fftduration<<"s\n";

    if(rank_in_pool==0) ModuleBase::timer::finish(GlobalV::ofs_running, true);
    
    if(tmp!=NULL) delete []tmp; 
    // MPI_Barrier(MPI_COMM_WORLD);
    // ProfilerStop();
    return 0;
}