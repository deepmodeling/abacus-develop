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

TEST_F(SccsPwCoulombTest, ForwardAndAdjointCallsDoNotRetainPreviousFields)
{
    const double tpiba = ModuleBase::TWO_PI / 10.0;
    const ModuleSccs::PeriodicCoulombOperator coulomb(basis_, tpiba);
    std::vector<double> charge(basis_.nrxx);
    std::vector<ModuleBase::Vector3<double>> probe(basis_.nrxx);
    for (int i = 0; i < basis_.nrxx; ++i)
    {
        charge[i] = std::sin(0.13 * i);
        probe[i].x = std::cos(0.17 * i);
        probe[i].y = std::sin(0.19 * i);
        probe[i].z = std::cos(0.23 * i);
    }
    ModuleSccs::ElectrostaticField expected;
    coulomb.apply(charge, expected);
    std::vector<double> expected_adjoint;
    coulomb.apply_gradient_adjoint(probe, expected_adjoint);
    const std::vector<double> zero_charge(basis_.nrxx, 0.0);
    const std::vector<ModuleBase::Vector3<double>> zero_probe(basis_.nrxx);
    ModuleSccs::ElectrostaticField field;
    std::vector<double> adjoint;
    for (int repeat = 0; repeat < 3; ++repeat)
    {
        coulomb.apply(zero_charge, field);
        coulomb.apply_gradient_adjoint(zero_probe, adjoint);
        for (int i = 0; i < basis_.nrxx; ++i)
        {
            EXPECT_DOUBLE_EQ(field.potential[i], 0.0);
            EXPECT_DOUBLE_EQ(adjoint[i], 0.0);
            for (int d = 0; d < 3; ++d)
            {
                EXPECT_DOUBLE_EQ(field.gradient[i][d], 0.0);
            }
        }
        coulomb.apply_gradient_adjoint(probe, adjoint);
        coulomb.apply(charge, field);
        for (int i = 0; i < basis_.nrxx; ++i)
        {
            EXPECT_DOUBLE_EQ(field.potential[i], expected.potential[i]);
            EXPECT_DOUBLE_EQ(adjoint[i], expected_adjoint[i]);
            for (int d = 0; d < 3; ++d)
            {
                EXPECT_DOUBLE_EQ(field.gradient[i][d], expected.gradient[i][d]);
            }
        }
    }
}

TEST_F(SccsPwCoulombTest, CountsForwardAndAdjointTransforms)
{
    const double tpiba = ModuleBase::TWO_PI / 10.0;
    const ModuleSccs::PeriodicCoulombOperator coulomb(basis_, tpiba);
    const std::vector<double> charge(basis_.nrxx, 1.0);
    ModuleSccs::ElectrostaticField field;
    coulomb.apply_gradient(charge, field.gradient);
    ModuleSccs::CoulombTransformCounts counts = coulomb.transform_counts();
    EXPECT_EQ(counts.forward_calls, 1);
    EXPECT_EQ(counts.inverse_calls, 3);
    coulomb.apply(charge, field);
    counts = coulomb.transform_counts();
    EXPECT_EQ(counts.forward_calls, 2);
    EXPECT_EQ(counts.inverse_calls, 7);
    std::vector<double> adjoint;
    coulomb.apply_gradient_adjoint(field.gradient, adjoint);
    counts = coulomb.transform_counts();
    EXPECT_EQ(counts.forward_calls, 5);
    EXPECT_EQ(counts.inverse_calls, 8);
}

TEST_F(SccsPwCoulombTest, ScalarPotentialMatchesFullFieldWithoutGradientTransforms)
{
    const double tpiba = ModuleBase::TWO_PI / 10.0;
    const ModuleSccs::PeriodicCoulombOperator coulomb(basis_, tpiba);
    std::vector<double> charge(basis_.nrxx);
    for (int ir = 0; ir < basis_.nrxx; ++ir)
    {
        charge[ir] = std::sin(0.13 * ir);
    }
    ModuleSccs::ElectrostaticField field;
    coulomb.apply(charge, field);
    std::vector<double> potential;
    coulomb.apply_potential(charge, potential);
    const ModuleSccs::CoulombTransformCounts counts = coulomb.transform_counts();
    EXPECT_EQ(counts.forward_calls, 2);
    EXPECT_EQ(counts.inverse_calls, 5);
    EXPECT_EQ(potential, field.potential);
    const std::vector<double> zero(basis_.nrxx, 0.0);
    coulomb.apply_potential(zero, potential);
    for (const double value : potential)
    {
        EXPECT_DOUBLE_EQ(value, 0.0);
    }
}

#ifdef _OPENMP
TEST_F(SccsPwCoulombTest, ParallelGridLoopsMatchSerialFields)
{
    const int previous_threads = omp_get_max_threads();
    const double tpiba = ModuleBase::TWO_PI / 10.0;
    const ModuleSccs::PeriodicCoulombOperator coulomb(basis_, tpiba);
    std::vector<double> charge(basis_.nrxx);
    for (int ir = 0; ir < basis_.nrxx; ++ir)
    {
        charge[ir] = std::sin(0.13 * ir);
    }
    ModuleSccs::ElectrostaticField serial;
    std::vector<double> serial_adjoint;
    omp_set_num_threads(1);
    coulomb.apply(charge, serial);
    coulomb.apply_gradient_adjoint(serial.gradient, serial_adjoint);
    ModuleSccs::ElectrostaticField parallel;
    std::vector<double> parallel_adjoint;
    omp_set_num_threads(2);
    coulomb.apply(charge, parallel);
    coulomb.apply_gradient_adjoint(parallel.gradient, parallel_adjoint);
    omp_set_num_threads(previous_threads);
    for (int ir = 0; ir < basis_.nrxx; ++ir)
    {
        EXPECT_DOUBLE_EQ(serial.potential[ir], parallel.potential[ir]);
        EXPECT_DOUBLE_EQ(serial_adjoint[ir], parallel_adjoint[ir]);
        for (int direction = 0; direction < 3; ++direction)
        {
            EXPECT_DOUBLE_EQ(serial.gradient[ir][direction], parallel.gradient[ir][direction]);
        }
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
        positions = ModuleSccs::pw_grid_positions(*basis, lattice, length);
        geometry = ModuleSccs::pcc_geometry(lattice, length, 1.0e-10);
        geometry_2d = ModuleSccs::pcc_2d_geometry(lattice, length, 1.0e-10);
        charge.resize(basis->nrxx);
        for (int i = 0; i < basis->nrxx; ++i)
        {
            const double x = positions[i].x - 6.0;
            const double y = positions[i].y - 6.0;
            const double z = positions[i].z - 6.0;
            charge[i] = -0.02 * std::exp(-((x - 0.7) * (x - 0.7) + y * y + z * z) / 2.0);
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
    ModuleSccs::PccGeometry geometry;
    ModuleSccs::Pcc2dGeometry geometry_2d;
    ModuleSccs::SerialChargeReduction charge_reduction;
    double tpiba = 0.0;
    double dv = 0.0;
};

TEST_F(SccsCoulombOperatorsTest, CoulombGradientAdjointsPreserveInnerProducts)
{
    for (int boundary = 0; boundary < 3; ++boundary)
    {
        SCOPED_TRACE(boundary);
        const auto op = make_operator(boundary);
        ModuleSccs::ElectrostaticField field;
        op->apply(charge, field);
        std::vector<ModuleBase::Vector3<double>> probe(basis->nrxx);
        for (int i = 0; i < basis->nrxx; ++i)
        {
            probe[i].x = 0.3 + std::sin(positions[i].x);
            probe[i].y = -0.5 + std::cos(positions[i].y);
            probe[i].z = positions[i].z - 6.0;
        }
        std::vector<double> transpose;
        op->apply_gradient_adjoint(probe, transpose);
        double left = 0.0;
        double right = 0.0;
        for (int i = 0; i < basis->nrxx; ++i)
        {
            left += (field.gradient[i] * probe[i]) * dv;
            right += charge[i] * transpose[i] * dv;
        }
        EXPECT_NEAR(left, right, 2.0e-12);
    }
}

TEST_F(SccsCoulombOperatorsTest, GradientOnlyMatchesFullFieldForAllBoundaries)
{
    for (int boundary = 0; boundary < 3; ++boundary)
    {
        SCOPED_TRACE(boundary);
        const auto op = make_operator(boundary);
        ModuleSccs::ElectrostaticField field;
        std::vector<ModuleBase::Vector3<double>> gradient;
        for (const double scale : {1.0, 0.0, -0.7})
        {
            auto source = charge;
            for (double& value : source)
            {
                value *= scale;
            }
            op->apply_gradient(source, gradient);
            op->apply(source, field);
            ASSERT_EQ(gradient.size(), field.gradient.size());
            for (std::size_t i = 0; i < gradient.size(); ++i)
            {
                for (int d = 0; d < 3; ++d)
                {
                    EXPECT_DOUBLE_EQ(gradient[i][d], field.gradient[i][d]);
                }
            }
        }
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
