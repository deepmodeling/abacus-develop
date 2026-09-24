#include "../sccs_poisson.h"

#include "source_base/constants.h"

#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace
{

class LocalResponseOperator : public ModuleSccs::CoulombOperator
{
  public:
    explicit LocalResponseOperator(const double response) : response_(response)
    {
    }

    void apply(const std::vector<double>& charge, ModuleSccs::ElectrostaticField& field) const override
    {
        field.potential = charge;
        field.gradient.assign(charge.size(), ModuleBase::Vector3<double>());
        for (std::size_t index = 0; index < charge.size(); ++index)
        {
            field.gradient[index].x = ModuleBase::FOUR_PI * response_ * charge[index];
        }
    }

  private:
    double response_ = 0.0;
};

class NonFiniteOperator : public ModuleSccs::CoulombOperator
{
  public:
    void apply(const std::vector<double>& charge, ModuleSccs::ElectrostaticField& field) const override
    {
        field.potential.assign(charge.size(), std::numeric_limits<double>::quiet_NaN());
        field.gradient.assign(charge.size(), ModuleBase::Vector3<double>());
    }
};

class TwoDomainReduction : public ModuleSccs::PolarizationReduction
{
  public:
    void reduce_residual(double& square_sum,
                         double& maximum,
                         double& point_count) const override
    {
        square_sum += 9.0;
        maximum = std::max(maximum, 3.0);
        point_count += 1.0;
    }
};

ModuleSccs::PolarizationSolverParameters converged_parameters()
{
    ModuleSccs::PolarizationSolverParameters parameters;
    parameters.max_iterations = 200;
    parameters.mixing = 0.7;
    parameters.tolerance_rms = 1.0e-12;
    parameters.tolerance_max = 1.0e-12;
    return parameters;
}

TEST(SccsPoisson, UniformDielectricHasAnalyticPolarizationCharge)
{
    const std::vector<double> solute_charge = {0.5, -0.25, 0.125, -0.375};
    const std::vector<double> epsilon(solute_charge.size(), 5.0);
    const std::vector<ModuleBase::Vector3<double>> grad_log_epsilon(solute_charge.size());
    const LocalResponseOperator coulomb(0.0);

    const ModuleSccs::PolarizationResult result
        = ModuleSccs::solve_polarization(solute_charge,
                                         epsilon,
                                         grad_log_epsilon,
                                         std::vector<double>(),
                                         converged_parameters(),
                                         coulomb);

    EXPECT_EQ(result.status, ModuleSccs::PolarizationStatus::Converged);
    EXPECT_GT(result.iterations, 1);
    for (std::size_t index = 0; index < solute_charge.size(); ++index)
    {
        EXPECT_NEAR(result.polarization_charge[index], -0.8 * solute_charge[index], 2.0e-12);
        EXPECT_NEAR(result.field.potential[index], 0.2 * solute_charge[index], 2.0e-12);
    }
}

TEST(SccsPoisson, VariableDielectricFixedPointMatchesAnalyticLocalModel)
{
    const std::vector<double> solute_charge = {0.4, -0.2, 0.1};
    const std::vector<double> epsilon(solute_charge.size(), 2.0);
    std::vector<ModuleBase::Vector3<double>> grad_log_epsilon(solute_charge.size());
    for (std::size_t index = 0; index < grad_log_epsilon.size(); ++index)
    {
        grad_log_epsilon[index].x = 1.0;
    }
    const double response = 0.2;
    const LocalResponseOperator coulomb(response);

    const ModuleSccs::PolarizationResult result
        = ModuleSccs::solve_polarization(solute_charge,
                                         epsilon,
                                         grad_log_epsilon,
                                         std::vector<double>(),
                                         converged_parameters(),
                                         coulomb);

    EXPECT_EQ(result.status, ModuleSccs::PolarizationStatus::Converged);
    const double expected_factor = (response - 0.5) / (1.0 - response);
    for (std::size_t index = 0; index < solute_charge.size(); ++index)
    {
        EXPECT_NEAR(result.polarization_charge[index], expected_factor * solute_charge[index], 2.0e-12);
    }
}

TEST(SccsPoisson, AcceleratedMixingReducesIterationCount)
{
    const std::vector<double> solute_charge = {0.4, -0.2, 0.1};
    const std::vector<double> epsilon(solute_charge.size(), 2.0);
    std::vector<ModuleBase::Vector3<double>> grad_log_epsilon(solute_charge.size());
    for (std::size_t index = 0; index < grad_log_epsilon.size(); ++index)
    {
        grad_log_epsilon[index].x = 1.0;
    }
    const double response = 0.95;
    const LocalResponseOperator coulomb(response);
    ModuleSccs::PolarizationSolverParameters parameters = converged_parameters();
    parameters.max_iterations = 1000;
    parameters.mixing = 0.5;
    parameters.mixing_history = 6;
    parameters.tolerance_rms = 1.0e-8;
    parameters.tolerance_max = 1.0e-8;

    parameters.mixing_method = "linear";
    const ModuleSccs::PolarizationResult linear
        = ModuleSccs::solve_polarization(solute_charge,
                                         epsilon,
                                         grad_log_epsilon,
                                         std::vector<double>(),
                                         parameters,
                                         coulomb);
    ASSERT_EQ(linear.status, ModuleSccs::PolarizationStatus::Converged);
    EXPECT_GT(linear.iterations, 500);
    EXPECT_DOUBLE_EQ(linear.final_mixing, parameters.mixing);
    EXPECT_EQ(linear.mixing_restarts, 0);

    const std::vector<std::string> accelerated_methods = {"pulay", "anderson"};
    for (std::size_t method = 0; method < accelerated_methods.size(); ++method)
    {
        parameters.mixing_method = accelerated_methods[method];
        const ModuleSccs::PolarizationResult accelerated
            = ModuleSccs::solve_polarization(solute_charge,
                                             epsilon,
                                             grad_log_epsilon,
                                             std::vector<double>(),
                                             parameters,
                                             coulomb);
        EXPECT_EQ(accelerated.status, ModuleSccs::PolarizationStatus::Converged)
            << accelerated_methods[method];
        EXPECT_LT(accelerated.iterations, linear.iterations)
            << accelerated_methods[method];
        EXPECT_LT(accelerated.iterations, 20)
            << accelerated_methods[method];
        for (std::size_t index = 0; index < solute_charge.size(); ++index)
        {
            EXPECT_NEAR(accelerated.polarization_charge[index],
                        linear.polarization_charge[index],
                        3.0e-7)
                << accelerated_methods[method];
        }
    }
}

TEST(SccsPoisson, AdaptiveMixingRaisesAndLowersTheDampingFactor)
{
    const std::vector<double> solute_charge(1, 1.0);
    const std::vector<double> epsilon(1, 2.0);
    std::vector<ModuleBase::Vector3<double>> gradient(1);
    ModuleSccs::PolarizationSolverParameters parameters = converged_parameters();
    parameters.mixing_method = "linear";
    parameters.adaptive_mixing = true;
    parameters.mixing = 0.5;
    parameters.mixing_min = 0.1;
    parameters.mixing_max = 0.8;
    parameters.max_iterations = 8;
    parameters.tolerance_rms = 1.0e-30;
    parameters.tolerance_max = 1.0e-30;

    const LocalResponseOperator contractive_coulomb(0.0);
    const ModuleSccs::PolarizationResult contractive
        = ModuleSccs::solve_polarization(solute_charge,
                                         epsilon,
                                         gradient,
                                         std::vector<double>(),
                                         parameters,
                                         contractive_coulomb);
    EXPECT_EQ(contractive.status, ModuleSccs::PolarizationStatus::MaxIterations);
    EXPECT_GT(contractive.final_mixing, parameters.mixing);
    EXPECT_LE(contractive.final_mixing, parameters.mixing_max);
    EXPECT_EQ(contractive.mixing_restarts, 0);

    gradient[0].x = 1.0;
    const LocalResponseOperator divergent_coulomb(4.0);
    parameters.max_iterations = 4;
    const ModuleSccs::PolarizationResult divergent
        = ModuleSccs::solve_polarization(solute_charge,
                                         epsilon,
                                         gradient,
                                         std::vector<double>(),
                                         parameters,
                                         divergent_coulomb);
    EXPECT_EQ(divergent.status, ModuleSccs::PolarizationStatus::MaxIterations);
    EXPECT_LT(divergent.final_mixing, parameters.mixing);
    EXPECT_GE(divergent.final_mixing, parameters.mixing_min);
    EXPECT_EQ(divergent.mixing_restarts, 2);
}

TEST(SccsPoisson, ResidualIsMeasuredBeforeMixing)
{
    const std::vector<double> solute_charge(1, 1.0);
    const std::vector<double> epsilon(1, 2.0);
    const std::vector<ModuleBase::Vector3<double>> grad_log_epsilon(1);
    const LocalResponseOperator coulomb(0.0);
    ModuleSccs::PolarizationSolverParameters parameters = converged_parameters();
    parameters.max_iterations = 1;
    parameters.mixing = 1.0e-6;

    const ModuleSccs::PolarizationResult result
        = ModuleSccs::solve_polarization(solute_charge,
                                         epsilon,
                                         grad_log_epsilon,
                                         std::vector<double>(),
                                         parameters,
                                         coulomb);

    EXPECT_EQ(result.status, ModuleSccs::PolarizationStatus::MaxIterations);
    EXPECT_DOUBLE_EQ(result.residual_rms, 0.5);
    EXPECT_DOUBLE_EQ(result.residual_max, 0.5);
}

TEST(SccsPoisson, ReportsNonFiniteCoulombOutput)
{
    const std::vector<double> solute_charge(2, 0.0);
    const std::vector<double> epsilon(2, 1.0);
    const std::vector<ModuleBase::Vector3<double>> grad_log_epsilon(2);
    const NonFiniteOperator coulomb;

    const ModuleSccs::PolarizationResult result
        = ModuleSccs::solve_polarization(solute_charge,
                                         epsilon,
                                         grad_log_epsilon,
                                         std::vector<double>(),
                                         converged_parameters(),
                                         coulomb);

    EXPECT_EQ(result.status, ModuleSccs::PolarizationStatus::NonFinite);
    EXPECT_EQ(result.iterations, 1);
}

TEST(SccsPoisson, UsesGlobalResidualReduction)
{
    const std::vector<double> solute_charge(1, 1.0);
    const std::vector<double> epsilon(1, 2.0);
    const std::vector<ModuleBase::Vector3<double>> grad_log_epsilon(1);
    const LocalResponseOperator coulomb(0.0);
    const TwoDomainReduction reduction;
    ModuleSccs::PolarizationSolverParameters parameters = converged_parameters();
    parameters.max_iterations = 1;

    const ModuleSccs::PolarizationResult result
        = ModuleSccs::solve_polarization(solute_charge,
                                         epsilon,
                                         grad_log_epsilon,
                                         std::vector<double>(),
                                         parameters,
                                         coulomb,
                                         reduction);

    EXPECT_EQ(result.status, ModuleSccs::PolarizationStatus::MaxIterations);
    EXPECT_NEAR(result.residual_rms, std::sqrt((0.25 + 9.0) / 2.0), 1.0e-14);
    EXPECT_DOUBLE_EQ(result.residual_max, 3.0);
}

TEST(SccsPoisson, RejectsInvalidInputs)
{
    const LocalResponseOperator coulomb(0.0);
    EXPECT_THROW(ModuleSccs::solve_polarization(std::vector<double>(),
                                                std::vector<double>(),
                                                std::vector<ModuleBase::Vector3<double>>(),
                                                std::vector<double>(),
                                                converged_parameters(),
                                                coulomb),
                 std::invalid_argument);

    std::vector<double> solute_charge(1, 0.0);
    std::vector<double> epsilon(1, 0.5);
    std::vector<ModuleBase::Vector3<double>> gradient(1);
    EXPECT_THROW(ModuleSccs::solve_polarization(solute_charge,
                                                epsilon,
                                                gradient,
                                                std::vector<double>(),
                                                converged_parameters(),
                                                coulomb),
                 std::domain_error);

    ModuleSccs::PolarizationSolverParameters invalid_parameters = converged_parameters();
    invalid_parameters.mixing_method = "unknown";
    EXPECT_THROW(ModuleSccs::solve_polarization(solute_charge,
                                                std::vector<double>(1, 1.0),
                                                gradient,
                                                std::vector<double>(),
                                                invalid_parameters,
                                                coulomb),
                 std::invalid_argument);
    invalid_parameters = converged_parameters();
    invalid_parameters.mixing_history = 1;
    EXPECT_THROW(ModuleSccs::solve_polarization(solute_charge,
                                                std::vector<double>(1, 1.0),
                                                gradient,
                                                std::vector<double>(),
                                                invalid_parameters,
                                                coulomb),
                 std::invalid_argument);
    invalid_parameters = converged_parameters();
    invalid_parameters.adaptive_mixing = true;
    invalid_parameters.mixing_min = 0.8;
    invalid_parameters.mixing_max = 0.2;
    EXPECT_THROW(ModuleSccs::solve_polarization(solute_charge,
                                                std::vector<double>(1, 1.0),
                                                gradient,
                                                std::vector<double>(),
                                                invalid_parameters,
                                                coulomb),
                 std::invalid_argument);
}

} // namespace
