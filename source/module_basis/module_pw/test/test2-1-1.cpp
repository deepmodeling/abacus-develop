#include "../pw_basis.h"
#ifdef __MPI
#include "test_tool.h"
#include "module_base/parallel_global.h"
#include "mpi.h"
#endif
#include "module_base/global_function.h"
#include "module_base/constants.h"
#include "pw_test.h"
extern int nproc_in_pool,rank_in_pool;
using namespace std;

TEST_F(PWTEST,test2_1_1)
{
    cout<<"dividemthd 2, gamma_only: off, check gcar,gdirect,gg,istot2ixy,ig2isz"<<endl;
    //--------------------------------------------------
    ModuleBase::Matrix3 latvec(1,20,0,1,1,0,0,0,5);
    bool gamma_only = false;
    double wfcecut = 70;
    double lat0 = 1;
    int distribution_type = 2;
    bool xprime = false;
    //--------------------------------------------------

    ModulePW::PW_Basis pwtest(device_flag, precision_flag);
#ifdef __MPI
    pwtest.initmpi(nproc_in_pool, rank_in_pool, POOL_WORLD);
#endif
    pwtest.initgrids(lat0, latvec, wfcecut);
    pwtest.initparameters(gamma_only, wfcecut, distribution_type,xprime);
    pwtest.setuptransform();
    pwtest.collect_local_pw();
    pwtest.collect_uniqgg();
    ModuleBase::Matrix3 GT,G,GGT;
    GT = latvec.Inverse();
	G  = GT.Transpose();
	GGT = G * GT;
    double tpiba2 = ModuleBase::TWO_PI * ModuleBase::TWO_PI / lat0 / lat0;
    double ggecut = wfcecut / tpiba2;

    //ref
    const int totnpw_ref = 911;
    const int totnst_ref = 95;
    const int nx_ref = 54;
    const int ny_ref = 3;
    const int nz_ref = 15;
    //some results for different number of processors
    int npw_per_ref[12][12]={
        {911},
        {468,443},
        {332,256,323},
        {248,220,204,239},
        {191,195,143,187,195},
        {156,176,136,120,176,147},
        {140,152,130,98,119,149,123},
        {116,132,122,98,82,122,140,99},
        {111,115,117,97,71,88,106,110,96},
        {102,102,108,90,72,71,83,99,99,85},
        {89,91,101,89,79,53,77,82,90,88,72},
        {82,74,92,84,70,66,52,68,84,92,80,67}
    };
    int nst_per_ref[12][12]={
        {95},
        {48,47},
        {32,32,31},
        {24,24,24,23},
        {19,19,19,19,19},
        {16,16,16,16,16,15},
        {14,14,14,14,13,13,13},
        {12,12,12,12,12,12,12,11},
        {11,11,11,11,11,10,10,10,10},
        {10,10,10,10,10,9,9,9,9,9},
        {9,9,9,9,9,9,9,8,8,8,8},
        {8,8,8,8,8,8,8,8,8,8,8,7}
    };
    int *npw_per = nullptr;
    if(rank_in_pool == 0)
    {
        npw_per = new int [nproc_in_pool];
    }
#ifdef __MPI
    MPI_Gather(&pwtest.npw,1,MPI_INT,npw_per,1,MPI_INT,0,POOL_WORLD);
#else
    if(rank_in_pool == 0) npw_per[0] = pwtest.npw;
#endif
    if(rank_in_pool == 0)
    {
        if(nproc_in_pool <= 12)
        {
            for(int ip = 0 ; ip < nproc_in_pool ; ++ip)
            {
                EXPECT_EQ(npw_per_ref[nproc_in_pool-1][ip], npw_per[ip]);
                EXPECT_EQ(nst_per_ref[nproc_in_pool-1][ip], pwtest.nst_per[ip]);
            }
        }
        else
        {
            cout<<"Please use mpi processors no more than 12."<<endl;
        }
        delete []npw_per;
    }
    

    //results
    int tot_npw = 0;
#ifdef __MPI
    MPI_Allreduce(&pwtest.npw, &tot_npw, 1, MPI_INT, MPI_SUM, POOL_WORLD);
#else
    tot_npw = pwtest.npw;
#endif
    EXPECT_EQ(pwtest.nx, nx_ref);
    EXPECT_EQ(pwtest.ny, ny_ref);
    EXPECT_EQ(pwtest.nz, nz_ref);
    EXPECT_EQ(pwtest.fftnx, nx_ref);
    EXPECT_EQ(pwtest.fftny, ny_ref);
    EXPECT_EQ(tot_npw, totnpw_ref);
    EXPECT_EQ(pwtest.npwtot, totnpw_ref);
    EXPECT_EQ(pwtest.nstot,totnst_ref);
    EXPECT_EQ(pwtest.nxyz, nx_ref*ny_ref*nz_ref);


    int *tmpx = new int[pwtest.nx*pwtest.ny*pwtest.nz];
    int *tmpy = new int[pwtest.nx*pwtest.ny*pwtest.nz];
    int *tmpz = new int[pwtest.nx*pwtest.ny*pwtest.nz];
    ModuleBase::GlobalFunc::ZEROS(tmpx,pwtest.nx*pwtest.ny*pwtest.nz);
    ModuleBase::GlobalFunc::ZEROS(tmpy,pwtest.nx*pwtest.ny*pwtest.nz);
    ModuleBase::GlobalFunc::ZEROS(tmpz,pwtest.nx*pwtest.ny*pwtest.nz);
    
    int * startnst = new int [nproc_in_pool];
    startnst[0] = 0;
    for(int ip = 1 ; ip < nproc_in_pool; ++ip)
    {
        startnst[ip] = startnst[ip-1] + pwtest.nst_per[ip-1];
    }

    for(int ig = 0 ; ig < pwtest.npw; ++ig)
    {
        int istot = pwtest.ig2isz[ig] / pwtest.nz + startnst[rank_in_pool];
        // int is = pwtest.ig2isz[ig] / pwtest.nz;
        int iz = pwtest.ig2isz[ig] % pwtest.nz;
        int iy = pwtest.istot2ixy[istot] % pwtest.ny;
        int ix = pwtest.istot2ixy[istot] / pwtest.ny;
        // int iy = pwtest.is2fftixy[is] % pwtest.ny;
        // int ix = pwtest.is2fftixy[is] / pwtest.ny;

        tmpx[iz+(iy+ix*pwtest.ny)*pwtest.nz] = int(pwtest.gdirect[ig].x);
        tmpy[iz+(iy+ix*pwtest.ny)*pwtest.nz] = int(pwtest.gdirect[ig].y);
        tmpz[iz+(iy+ix*pwtest.ny)*pwtest.nz] = int(pwtest.gdirect[ig].z);
    }
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE,tmpx,pwtest.nxyz,MPI_INT,MPI_SUM,POOL_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,tmpy,pwtest.nxyz,MPI_INT,MPI_SUM,POOL_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,tmpz,pwtest.nxyz,MPI_INT,MPI_SUM,POOL_WORLD);
#endif
    if(rank_in_pool==0)
    {
        for(int iz = 0 ; iz < pwtest.nz; ++iz)
        {
            for(int iy = 0 ; iy < pwtest.ny ; ++iy)
            {
                for(int ix = 0 ; ix < pwtest.nx ; ++ix)
                {
                    ModuleBase::Vector3<double> f;
                    f.x = ix;
                    f.y = iy;
                    f.z = iz;
                    if(iz >= int(pwtest.nz/2) +1) f.z -= pwtest.nz;
                    if(iy >= int(pwtest.ny/2) +1) f.y -= pwtest.ny;
                    if(ix >= int(pwtest.nx/2) +1) f.x -= pwtest.nx;
                    double modulus = f * (GGT * f);
                    if (modulus <= ggecut)
                    {
                        EXPECT_EQ(tmpx[iz + iy*pwtest.nz + ix*pwtest.ny*pwtest.nz], int(f.x));
                        EXPECT_EQ(tmpy[iz + iy*pwtest.nz + ix*pwtest.ny*pwtest.nz], int(f.y));
                        EXPECT_EQ(tmpz[iz + iy*pwtest.nz + ix*pwtest.ny*pwtest.nz], int(f.z));
                    }
                    
                }
            }
        }
    }
    for(int ig = 0 ;ig < pwtest.npw ; ++ig)
    {
        ModuleBase::Vector3<double> f;
        f.x = pwtest.gdirect[ig].x;
        f.y = pwtest.gdirect[ig].y;
        f.z = pwtest.gdirect[ig].z;
        ModuleBase::Vector3<double> gcar;
        gcar = f * G;
        double modulus = f*GGT*f;
        EXPECT_NEAR(gcar.x,pwtest.gcar[ig].x,1e-6);
        EXPECT_NEAR(gcar.y,pwtest.gcar[ig].y,1e-6);
        EXPECT_NEAR(gcar.z,pwtest.gcar[ig].z,1e-6);
        EXPECT_NEAR(modulus,pwtest.gg[ig],1e-6);
        EXPECT_NEAR(pwtest.gg[ig], pwtest.gg_uniq[pwtest.ig2igg[ig]],1e-8);
    }
    for(int igg = 1 ; igg < pwtest.ngg ; ++igg)
    {
        EXPECT_GT(pwtest.gg_uniq[igg], pwtest.gg_uniq[igg-1]);
    }
    if(pwtest.ig_gge0 >= 0) {EXPECT_NEAR(0.0, pwtest.gg[pwtest.ig_gge0], 1e-8);}
    delete [] startnst;
    delete [] tmpx;
    delete [] tmpy;
    delete [] tmpz;
}