#ifdef __MPI
#include "source_base/parallel_global.h"
#include <mpi.h>
#endif

#include "../sccs/sccs_pw_coulomb.h"

#include "source_base/constants.h"
#include "source_base/matrix3.h"
#include "source_basis/module_pw/pw_basis.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <vector>

namespace
{

class SccsPwCoulombTest : public testing::Test
{
  protected:
    SccsPwCoulombTest() : basis_("cpu", "double")
    {
    }

    void SetUp() override
    {
#ifdef __MPI
        basis_.initmpi(1, 0, POOL_WORLD);
#endif
        const ModuleBase::Matrix3 lattice(1.0, 0.0, 0.0,
                                          0.0, 1.0, 0.0,
                                          0.0, 0.0, 1.0);
        basis_.initgrids(10.0, lattice, 20.0);
        basis_.initparameters(false, 20.0, 1, false);
        basis_.setuptransform();
        basis_.collect_local_pw();
    }

    ModulePW::PW_Basis basis_;
};

TEST_F(SccsPwCoulombTest, RemovesConstantPeriodicMode)
{
    const ModuleSccs::PeriodicCoulombOperator coulomb(basis_, ModuleBase::TWO_PI / 10.0);
    const std::vector<double> charge(basis_.nrxx, 1.0);
    ModuleSccs::ElectrostaticField field;
    coulomb.apply(charge, field);

    ASSERT_EQ(field.potential.size(), charge.size());
    ASSERT_EQ(field.gradient.size(), charge.size());
    for (int ir = 0; ir < basis_.nrxx; ++ir)
    {
        EXPECT_NEAR(field.potential[ir], 0.0, 1.0e-12);
        EXPECT_NEAR(field.gradient[ir].x, 0.0, 1.0e-12);
        EXPECT_NEAR(field.gradient[ir].y, 0.0, 1.0e-12);
        EXPECT_NEAR(field.gradient[ir].z, 0.0, 1.0e-12);
    }
}

TEST_F(SccsPwCoulombTest, SolvesSingleFourierShell)
{
    std::vector<std::complex<double>> charge_g(basis_.npw);
    double selected_gg = 0.0;
    for (int ig = 0; ig < basis_.npw; ++ig)
    {
        const double gx = basis_.gdirect[ig].x;
        const double gy = basis_.gdirect[ig].y;
        const double gz = basis_.gdirect[ig].z;
        if (std::abs(std::abs(gx) - 1.0) < 1.0e-12 && std::abs(gy) < 1.0e-12
            && std::abs(gz) < 1.0e-12)
        {
            charge_g[ig] = 0.5;
            selected_gg = basis_.gg[ig];
        }
    }
    ASSERT_GT(selected_gg, 0.0);

    std::vector<double> charge(basis_.nrxx);
    basis_.recip2real(charge_g.data(), charge.data());
    const double tpiba = ModuleBase::TWO_PI / 10.0;
    const double factor = ModuleBase::FOUR_PI / (tpiba * tpiba * selected_gg);
    const ModuleSccs::PeriodicCoulombOperator coulomb(basis_, tpiba);
    ModuleSccs::ElectrostaticField field;
    coulomb.apply(charge, field);

    double maximum_error = 0.0;
    for (int ir = 0; ir < basis_.nrxx; ++ir)
    {
        maximum_error = std::max(maximum_error, std::abs(field.potential[ir] - factor * charge[ir]));
    }
    EXPECT_LT(maximum_error, 1.0e-11);
}

} // namespace

int main(int argc, char** argv)
{
#ifdef __MPI
    int process_count = 1;
    int thread_count = 1;
    int rank = 0;
    Parallel_Global::read_pal_param(argc, argv, process_count, thread_count, rank);
    POOL_WORLD = MPI_COMM_WORLD;
    KP_WORLD = MPI_COMM_NULL;
    INT_BGROUP = MPI_COMM_NULL;
    BP_WORLD = MPI_COMM_NULL;
    GRID_WORLD = MPI_COMM_NULL;
    DIAG_WORLD = MPI_COMM_NULL;
#endif
    testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
#ifdef __MPI
    Parallel_Global::finalize_mpi();
#endif
    return result;
}
