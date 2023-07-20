#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "module_base/matrix3.h"
#include "module_base/parallel_global.h"
#include "module_elecstate/module_charge/charge.h"
Charge::Charge()
{
}
Charge::~Charge()
{
    delete[] rec;
    delete[] dis;
}
namespace elecstate
{
int tmp_xc_func_type = 3;
int get_xc_func_type()
{
    return tmp_xc_func_type;
}
} // namespace elecstate

auto sum_array = [](const double* v, const int& nv) {
    double sum = 0;
    for (int i = 0; i < nv; ++i)
        sum += v[i];
    return sum;
};
/************************************************
 *  unit test of module_charge/charge_mpi.cpp
 ***********************************************/

/**
 * - Tested Functions:
 *   - rho_mpi: Charge::rho_mpi():
 *     - test rho_mpi
 *   - reduce_diff_pools: Charge::reduce_diff_pools()
 *     - test reduce_diff_pools
 *     - using rhopw and GlobalV
 */

class ChargeMpiTest : public ::testing::Test
{
  protected:
    Charge* charge;
    std::string output;
    double lat0 = 4;
    ModuleBase::Matrix3 latvec;
    void SetUp() override
    {
        charge = new Charge;
    }
    void TearDown() override
    {
        delete charge;
    }
};

TEST_F(ChargeMpiTest, reduce_diff_pools1)
{
    if (GlobalV::NPROC >= 2 && GlobalV::NPROC % 2 == 0)
    {
        GlobalV::KPAR = 2;
        Parallel_Global::divide_pools();
        ModulePW::PW_Basis* rhopw = new ModulePW::PW_Basis();
        rhopw->initmpi(GlobalV::NPROC_IN_POOL, GlobalV::RANK_IN_POOL, POOL_WORLD);
        rhopw->initgrids(lat0, latvec, 40);
        rhopw->initparameters(false, 10);
        rhopw->setuptransform();

        int nz = rhopw->nz;
        const int nrxx = rhopw->nrxx;
        const int nxy = rhopw->nxy;
        const int nplane = rhopw->nplane;
        charge->nrxx = nrxx;
        double* array_rho = new double[nrxx];
        for (int ir = 0; ir < nxy; ++ir)
        {
            for (int iz = 0; iz < nplane; ++iz)
            {
                array_rho[nplane * ir + iz] = (rhopw->startz_current + iz + ir * nz) / double(nxy * nz);
            }
        }
        double refsum = sum_array(array_rho, nrxx);

        charge->init_chgmpi(nz, 1);
        charge->reduce_diff_pools(array_rho);
        double sum = sum_array(array_rho, nrxx);
        EXPECT_EQ(sum, refsum * GlobalV::KPAR);

        delete[] array_rho;
        delete rhopw;
    }
}

TEST_F(ChargeMpiTest, reduce_diff_pools2)
{
    if (GlobalV::NPROC >= 3)
    {
        GlobalV::KPAR = 3;
        Parallel_Global::divide_pools();
        ModulePW::PW_Basis* rhopw = new ModulePW::PW_Basis();
        rhopw->initmpi(GlobalV::NPROC_IN_POOL, GlobalV::RANK_IN_POOL, POOL_WORLD);
        rhopw->initgrids(lat0, latvec, 40);
        rhopw->initparameters(false, 10);
        rhopw->setuptransform();
        charge->rhopw = rhopw;

        int nz = rhopw->nz;
        const int nrxx = rhopw->nrxx;
        const int nxy = rhopw->nxy;
        const int nplane = rhopw->nplane;
        charge->nrxx = nrxx;
        double* array_ref = new double[nxy * nz];
        for (int ixyz = 0; ixyz < nxy * nz; ++ixyz)
        {
            array_ref[ixyz] = ixyz / double(nxy * nz);
        }
        double refsum = sum_array(array_ref, nxy * nz);

        double* array_rho = new double[nrxx];
        for (int ir = 0; ir < nxy; ++ir)
        {
            for (int iz = 0; iz < nplane; ++iz)
            {
                array_rho[nplane * ir + iz] = (rhopw->startz_current + iz + ir * nz) / double(nxy * nz);
            }
        }

        charge->init_chgmpi(nz, 1);
        charge->reduce_diff_pools(array_rho);
        double sum = sum_array(array_rho, nrxx);
        MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE, MPI_SUM, POOL_WORLD);
        EXPECT_EQ(sum, refsum * GlobalV::KPAR);

        delete[] array_rho;
        delete rhopw;
        delete[] array_ref;
    }
}

TEST_F(ChargeMpiTest, rho_mpi)
{
    if (GlobalV::NPROC >= 2)
    {
        GlobalV::KPAR = 2;
        Parallel_Global::divide_pools();
        ModulePW::PW_Basis* rhopw = new ModulePW::PW_Basis();
        rhopw->initmpi(GlobalV::NPROC_IN_POOL, GlobalV::RANK_IN_POOL, POOL_WORLD);
        rhopw->initgrids(lat0, latvec, 40);
        rhopw->initparameters(false, 10);
        rhopw->setuptransform();
        charge->rhopw = rhopw;
        GlobalV::NSPIN = 1;
        charge->rho = new double*[1];
        charge->kin_r = new double*[1];

        int nz = rhopw->nz;
        const int nrxx = rhopw->nrxx;
        const int nxy = rhopw->nxy;
        const int nplane = rhopw->nplane;
        charge->nrxx = nrxx;
        charge->rho[0] = new double[nrxx];
        charge->kin_r[0] = new double[nrxx];
        charge->init_chgmpi(nz, 1);
        charge->rho_mpi();

        delete[] charge->rho[0];
        delete[] charge->rho;
        delete[] charge->kin_r[0];
        delete[] charge->kin_r;
        delete rhopw;
    }

    GlobalV::KPAR = 1;
    charge->rho_mpi();
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &GlobalV::NPROC);
    MPI_Comm_rank(MPI_COMM_WORLD, &GlobalV::MY_RANK);
    int result = 0;
    testing::InitGoogleTest(&argc, argv);
    result = RUN_ALL_TESTS();
    MPI_Finalize();
    return result;
}
