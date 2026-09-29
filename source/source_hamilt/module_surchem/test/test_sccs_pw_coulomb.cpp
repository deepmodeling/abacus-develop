#ifdef __MPI
#include "source_base/parallel_global.h"
#include <mpi.h>
#endif

#include "../sccs/sccs_pw_coulomb.h"
#include "../common/pw_grid.h"
#include "../sccs/sccs_charge.h"
#include "../sccs/sccs_pcc_coulomb.h"
#include "../sccs/sccs_pcc_2d_coulomb.h"

#include "source_base/constants.h"
#include "source_base/matrix3.h"
#include "source_basis/module_pw/pw_basis.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <memory>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

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
    std::vector<double> potential;
    coulomb.apply_potential(charge, potential);

    ASSERT_EQ(potential.size(), charge.size());
    for (int ir = 0; ir < basis_.nrxx; ++ir)
    {
        EXPECT_NEAR(potential[ir], 0.0, 1.0e-12);
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
    std::vector<double> potential;
    coulomb.apply_potential(charge, potential);

    double maximum_error = 0.0;
    for (int ir = 0; ir < basis_.nrxx; ++ir)
    {
        maximum_error = std::max(maximum_error, std::abs(potential[ir] - factor * charge[ir]));
    }
    EXPECT_LT(maximum_error, 1.0e-11);
}

TEST_F(SccsPwCoulombTest, RepeatedCallsDoNotRetainPreviousPotentials)
{
    const double tpiba = ModuleBase::TWO_PI / 10.0;
    const ModuleSccs::PeriodicCoulombOperator coulomb(basis_, tpiba);
    std::vector<double> charge(basis_.nrxx);
    for (int i = 0; i < basis_.nrxx; ++i)
    {
        charge[i] = std::sin(0.13 * i);
    }
    std::vector<double> expected;
    coulomb.apply_potential(charge, expected);
    const std::vector<double> zero_charge(basis_.nrxx, 0.0);
    std::vector<double> potential;
    for (int repeat = 0; repeat < 3; ++repeat)
    {
        coulomb.apply_potential(zero_charge, potential);
        for (int i = 0; i < basis_.nrxx; ++i)
        {
            EXPECT_DOUBLE_EQ(potential[i], 0.0);
        }
        coulomb.apply_potential(charge, potential);
        for (int i = 0; i < basis_.nrxx; ++i)
        {
            EXPECT_DOUBLE_EQ(potential[i], expected[i]);
        }
    }
}

TEST_F(SccsPwCoulombTest, CountsOneTransformPairPerPotential)
{
    const double tpiba = ModuleBase::TWO_PI / 10.0;
    const ModuleSccs::PeriodicCoulombOperator coulomb(basis_, tpiba);
    const std::vector<double> charge(basis_.nrxx, 1.0);
    std::vector<double> potential;
    coulomb.apply_potential(charge, potential);
    ModuleSccs::CoulombTransformCounts counts = coulomb.transform_counts();
    EXPECT_EQ(counts.forward_calls, 1);
    EXPECT_EQ(counts.inverse_calls, 1);
    coulomb.apply_potential(charge, potential);
    counts = coulomb.transform_counts();
    EXPECT_EQ(counts.forward_calls, 2);
    EXPECT_EQ(counts.inverse_calls, 2);
}

#ifdef _OPENMP
TEST_F(SccsPwCoulombTest, ParallelGridLoopsMatchSerialPotential)
{
    const int previous_threads = omp_get_max_threads();
    const double tpiba = ModuleBase::TWO_PI / 10.0;
    const ModuleSccs::PeriodicCoulombOperator coulomb(basis_, tpiba);
    std::vector<double> charge(basis_.nrxx);
    for (int ir = 0; ir < basis_.nrxx; ++ir)
    {
        charge[ir] = std::sin(0.13 * ir);
    }
    std::vector<double> serial;
    omp_set_num_threads(1);
    coulomb.apply_potential(charge, serial);
    std::vector<double> parallel;
    omp_set_num_threads(2);
    coulomb.apply_potential(charge, parallel);
    omp_set_num_threads(previous_threads);
    for (int ir = 0; ir < basis_.nrxx; ++ir)
    {
        EXPECT_DOUBLE_EQ(serial[ir], parallel[ir]);
    }
}
#endif

class SccsCoulombOperatorsTest : public testing::Test
{
  protected:
    void SetUp() override
    {
        basis.reset(new ModulePW::PW_Basis("cpu", "double"));
#ifdef __MPI
        basis->initmpi(1, 0, POOL_WORLD);
#endif
        const ModuleBase::Matrix3 lattice(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
        const double length = 12.0;
        basis->initgrids(length, lattice, 30.0);
        basis->initparameters(false, 30.0, 1, false);
        basis->setuptransform();
        basis->collect_local_pw();
        tpiba = ModuleBase::TWO_PI / length;
        dv = length * length * length / basis->nxyz;
        positions = ModuleSurchem::pw_grid_positions(*basis, lattice, length);
        geometry = ModulePcc::pcc_geometry(lattice, length, 1.0e-10);
        geometry_2d = ModulePcc::pcc_2d_geometry(lattice, length, 1.0e-10);
        charge.resize(basis->nrxx);
        other_charge.resize(basis->nrxx);
        for (int i = 0; i < basis->nrxx; ++i)
        {
            const double x = positions[i].x - 6.0;
            const double y = positions[i].y - 6.0;
            const double z = positions[i].z - 6.0;
            charge[i] = -0.02 * std::exp(-((x - 0.7) * (x - 0.7) + y * y + z * z) / 2.0);
            const double shifted_square = x * x + (y + 1.1) * (y + 1.1) + (z - 0.4) * (z - 0.4);
            other_charge[i] = (0.01 + 0.004 * y) * std::exp(-shifted_square / 3.0);
        }
    }

    std::unique_ptr<ModuleSccs::CoulombOperator> make_operator(const int boundary)
    {
        if (boundary == 0)
        {
            return std::unique_ptr<ModuleSccs::CoulombOperator>(
                new ModuleSccs::PeriodicCoulombOperator(*basis, tpiba));
        }
        if (boundary == 1)
        {
            return std::unique_ptr<ModuleSccs::CoulombOperator>(
                new ModuleSccs::PccCoulombOperator(*basis, tpiba, positions, dv, geometry, charge_reduction));
        }
        return std::unique_ptr<ModuleSccs::CoulombOperator>(
            new ModuleSccs::Pcc2dCoulombOperator(*basis, tpiba, positions, dv, geometry_2d, charge_reduction));
    }

    std::unique_ptr<ModulePW::PW_Basis> basis;
    std::vector<ModuleBase::Vector3<double>> positions;
    std::vector<double> charge;
    std::vector<double> other_charge;
    ModulePcc::PccGeometry geometry;
    ModulePcc::Pcc2dGeometry geometry_2d;
    ModuleSurchem::SerialChargeReduction charge_reduction;
    double tpiba = 0.0;
    double dv = 0.0;
};

// The sqrt-CG energy q^T A^-1 q / 2 needs A = sqrt(eps) G^-1 sqrt(eps) + F to
// be symmetric, so every Coulomb operator G, PCC terms included, must be.
TEST_F(SccsCoulombOperatorsTest, CoulombOperatorsAreSymmetric)
{
    for (int boundary = 0; boundary < 3; ++boundary)
    {
        SCOPED_TRACE(boundary);
        const std::unique_ptr<ModuleSccs::CoulombOperator> op = make_operator(boundary);
        std::vector<double> potential;
        std::vector<double> other_potential;
        op->apply_potential(charge, potential);
        op->apply_potential(other_charge, other_potential);
        double left = 0.0;
        double right = 0.0;
        for (int i = 0; i < basis->nrxx; ++i)
        {
            left += other_charge[i] * potential[i] * dv;
            right += charge[i] * other_potential[i] * dv;
        }
        const double scale = std::abs(left);
        const double tolerance = 1.0e-12 * scale;
        EXPECT_GT(scale, 1.0e-6);
        EXPECT_NEAR(left, right, tolerance);
    }
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
