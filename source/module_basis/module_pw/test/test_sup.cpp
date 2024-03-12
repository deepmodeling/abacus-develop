#include "../pw_basis.h"
#ifdef __MPI
#include "module_base/parallel_global.h"
#include "mpi.h"
#include "test_tool.h"
#endif
#include "module_base/constants.h"
#include "module_base/global_function.h"
#include "pw_test.h"
extern int nproc_in_pool, rank_in_pool;
using namespace std;

TEST_F(PWTEST, test_sup)
{
    cout << "dividemthd 3, gamma_only: off, check gcar,gdirect,gg,istot2ixy,ig2isz" << endl;
    //--------------------------------------------------
    ModuleBase::Matrix3 latvec(20, 1, 0, 0, 1, 0, 0, 0, 5);
    bool gamma_only = false;
    double wfcecut = 70;
    double wfcecutdense = 100;
    double lat0 = 1;
    int distribution_type = 1;
    bool xprime = false;
    //--------------------------------------------------

    // smooth grids
    ModulePW::PW_Basis pwsmooth(device_flag, precision_flag);
#ifdef __MPI
    pwsmooth.initmpi(nproc_in_pool, rank_in_pool, POOL_WORLD);
#endif
    pwsmooth.initgrids(lat0, latvec, wfcecut);
    pwsmooth.initparameters(gamma_only, wfcecut, distribution_type, xprime);
    pwsmooth.setuptransform();
    pwsmooth.collect_local_pw();
    pwsmooth.collect_uniqgg();

    // dense grids
    ModulePW::PW_Basis_Sup pwdense(device_flag, precision_flag);
#ifdef __MPI
    pwdense.initmpi(nproc_in_pool, rank_in_pool, POOL_WORLD);
#endif
    pwdense.initgrids(lat0, latvec, wfcecutdense);
    pwdense.initparameters(gamma_only, wfcecutdense, distribution_type, xprime);
    pwdense.setuptransform(&pwsmooth);
    pwdense.collect_local_pw();
    pwdense.collect_uniqgg();

    ModuleBase::Matrix3 GT, G, GGT;
    GT = latvec.Inverse();
    G = GT.Transpose();
    GGT = G * GT;
    double tpiba2 = ModuleBase::TWO_PI * ModuleBase::TWO_PI / lat0 / lat0;
    double ggecut = wfcecutdense / tpiba2;

    // ref
    const int totnpw_ref = 1751;
    const int totnst_ref = 161;
    const int nx_ref = 64;
    const int ny_ref = 3;
    const int nz_ref = 15;

    // some results for different number of processors
    int npw_per_ref[12][12] = {
        {1751},
        {874, 877},
        {583, 584, 584},
        {439, 436, 438, 438},
        {351, 350, 350, 350, 350},
        {293, 293, 293, 291, 291, 290},
        {249, 251, 251, 251, 251, 249, 249},
        {221, 218, 220, 218, 218, 218, 218, 220},
        {196, 196, 196, 194, 194, 194, 194, 193, 194},
        {177, 176, 176, 176, 174, 174, 174, 176, 174, 174},
        {161, 161, 161, 159, 159, 158, 158, 158, 159, 158, 159},
        {147, 148, 148, 145, 146, 146, 145, 145, 145, 146, 145, 145}
    };
    int nst_per_ref[12][12] = {
        {161},
        {80, 81},
        {53, 54, 54},
        {41, 40, 40, 40},
        {33, 32, 32, 32, 32},
        {27, 27, 27, 27, 27, 26},
        {23, 23, 23, 23, 23, 23, 23},
        {21, 20, 20, 20, 20, 20, 20, 20},
        {18, 18, 18, 18, 18, 18, 18, 17, 18},
        {17, 16, 16, 16, 16, 16, 16, 16, 16, 16},
        {15, 15, 15, 15, 15, 14, 14, 14, 15, 14, 15},
        {13, 14, 14, 13, 14, 14, 13, 13, 13, 14, 13, 13}
    };
    int* npw_per = nullptr;
    if (rank_in_pool == 0)
    {
        npw_per = new int[nproc_in_pool];
    }
#ifdef __MPI
    MPI_Gather(&pwdense.npw, 1, MPI_INT, npw_per, 1, MPI_INT, 0, POOL_WORLD);
#else
    if (rank_in_pool == 0)
        npw_per[0] = pwdense.npw;
#endif
    if (rank_in_pool == 0)
    {
        if (nproc_in_pool <= 12)
        {
            for (int ip = 0; ip < nproc_in_pool; ++ip)
            {
                EXPECT_EQ(npw_per_ref[nproc_in_pool - 1][ip], npw_per[ip]);
                EXPECT_EQ(nst_per_ref[nproc_in_pool - 1][ip], pwdense.nst_per[ip]);
            }
        }
        else
        {
            cout << "Please use mpi processors no more than 12." << endl;
        }
        delete[] npw_per;
    }

    // results
    int tot_npw = 0;
#ifdef __MPI
    MPI_Allreduce(&pwdense.npw, &tot_npw, 1, MPI_INT, MPI_SUM, POOL_WORLD);
#else
    tot_npw = pwdense.npw;
#endif
    EXPECT_EQ(pwdense.nx, nx_ref);
    EXPECT_EQ(pwdense.ny, ny_ref);
    EXPECT_EQ(pwdense.nz, nz_ref);
    EXPECT_EQ(pwdense.fftnx, nx_ref);
    EXPECT_EQ(pwdense.fftny, ny_ref);
    EXPECT_EQ(tot_npw, totnpw_ref);
    EXPECT_EQ(pwdense.npwtot, totnpw_ref);
    EXPECT_EQ(pwdense.nstot, totnst_ref);
    EXPECT_EQ(pwdense.nxyz, nx_ref * ny_ref * nz_ref);

    int* tmpx = new int[pwdense.nx * pwdense.ny * pwdense.nz];
    int* tmpy = new int[pwdense.nx * pwdense.ny * pwdense.nz];
    int* tmpz = new int[pwdense.nx * pwdense.ny * pwdense.nz];
    ModuleBase::GlobalFunc::ZEROS(tmpx, pwdense.nx * pwdense.ny * pwdense.nz);
    ModuleBase::GlobalFunc::ZEROS(tmpy, pwdense.nx * pwdense.ny * pwdense.nz);
    ModuleBase::GlobalFunc::ZEROS(tmpz, pwdense.nx * pwdense.ny * pwdense.nz);

    int* startnst = new int[nproc_in_pool];
    startnst[0] = 0;
    for (int ip = 1; ip < nproc_in_pool; ++ip)
    {
        startnst[ip] = startnst[ip - 1] + pwdense.nst_per[ip - 1];
    }

    for (int ig = 0; ig < pwdense.npw; ++ig)
    {
        int istot = pwdense.ig2isz[ig] / pwdense.nz + startnst[rank_in_pool];
        // int is = pwdense.ig2isz[ig] / pwdense.nz;
        int iz = pwdense.ig2isz[ig] % pwdense.nz;
        int iy = pwdense.istot2ixy[istot] % pwdense.ny;
        int ix = pwdense.istot2ixy[istot] / pwdense.ny;
        // int iy = pwdense.is2fftixy[is] % pwdense.ny;
        // int ix = pwdense.is2fftixy[is] / pwdense.ny;

        tmpx[iz + (iy + ix * pwdense.ny) * pwdense.nz] = int(pwdense.gdirect[ig].x);
        tmpy[iz + (iy + ix * pwdense.ny) * pwdense.nz] = int(pwdense.gdirect[ig].y);
        tmpz[iz + (iy + ix * pwdense.ny) * pwdense.nz] = int(pwdense.gdirect[ig].z);
    }
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, tmpx, pwdense.nxyz, MPI_INT, MPI_SUM, POOL_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, tmpy, pwdense.nxyz, MPI_INT, MPI_SUM, POOL_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, tmpz, pwdense.nxyz, MPI_INT, MPI_SUM, POOL_WORLD);
#endif
    if (rank_in_pool == 0)
    {
        for (int iz = 0; iz < pwdense.nz; ++iz)
        {
            for (int iy = 0; iy < pwdense.ny; ++iy)
            {
                for (int ix = 0; ix < pwdense.nx; ++ix)
                {
                    ModuleBase::Vector3<double> f;
                    f.x = ix;
                    f.y = iy;
                    f.z = iz;
                    if (iz >= int(pwdense.nz / 2) + 1)
                        f.z -= pwdense.nz;
                    if (iy >= int(pwdense.ny / 2) + 1)
                        f.y -= pwdense.ny;
                    if (ix >= int(pwdense.nx / 2) + 1)
                        f.x -= pwdense.nx;
                    double modulus = f * (GGT * f);
                    if (modulus <= ggecut)
                    {
                        EXPECT_EQ(tmpx[iz + iy * pwdense.nz + ix * pwdense.ny * pwdense.nz], int(f.x));
                        EXPECT_EQ(tmpy[iz + iy * pwdense.nz + ix * pwdense.ny * pwdense.nz], int(f.y));
                        EXPECT_EQ(tmpz[iz + iy * pwdense.nz + ix * pwdense.ny * pwdense.nz], int(f.z));
                    }
                }
            }
        }
    }
    for (int ig = 0; ig < pwdense.npw; ++ig)
    {
        ModuleBase::Vector3<double> f;
        f.x = pwdense.gdirect[ig].x;
        f.y = pwdense.gdirect[ig].y;
        f.z = pwdense.gdirect[ig].z;
        ModuleBase::Vector3<double> gcar;
        gcar = f * G;
        double modulus = f * GGT * f;
        EXPECT_NEAR(gcar.x, pwdense.gcar[ig].x, 1e-6);
        EXPECT_NEAR(gcar.y, pwdense.gcar[ig].y, 1e-6);
        EXPECT_NEAR(gcar.z, pwdense.gcar[ig].z, 1e-6);
        EXPECT_NEAR(modulus, pwdense.gg[ig], 1e-6);
        EXPECT_NEAR(pwdense.gg[ig], pwdense.gg_uniq[pwdense.ig2igg[ig]], 1e-8);
    }
    for (int igg = 1; igg < pwdense.ngg; ++igg)
    {
        EXPECT_GT(pwdense.gg_uniq[igg], pwdense.gg_uniq[igg - 1]);
    }
    if (pwdense.ig_gge0 >= 0)
    {
        EXPECT_NEAR(0.0, pwdense.gg[pwdense.ig_gge0], 1e-8);
    }
    delete[] startnst;
    delete[] tmpx;
    delete[] tmpy;
    delete[] tmpz;

    // the planewaves of dense grids must be consistent with the smooth grids
    for (int ig = 0; ig < pwsmooth.npw; ++ig)
    {
        for (int ipol = 0; ipol < 3; ipol++)
        {
            EXPECT_DOUBLE_EQ(pwsmooth.gcar[ig][ipol], pwdense.gcar[ig][ipol]);
        }
    }
}