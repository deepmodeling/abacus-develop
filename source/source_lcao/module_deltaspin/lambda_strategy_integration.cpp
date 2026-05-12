#include "spin_constrain.h"

#include "lambda_update_strategies.h"

/**
 * @file lambda_strategy_integration.cpp
 * @brief Integration of alternative lambda strategies into SpinConstrain.
 *
 * @par Status: INCOMPLETE
 * This file references members (strategy_type_, strategy_) that are NOT
 * declared in spin_constrain.h. The code will not compile as-is.
 * To enable, add the following to spin_constrain.h private section:
 *
 *   enum class LambdaStrategyType { BFGS, LinearResponse, AugmentedLagrangian, HybridDelayed };
 *   LambdaStrategyType strategy_type_ = LambdaStrategyType::BFGS;
 *   std::unique_ptr<LambdaUpdateStrategy> strategy_;
 *
 * And add the file to CMakeLists.txt.
 *
 * @par Purpose
 * Bridges the alternative strategy implementations (lambda_update_strategies.h)
 * to the SpinConstrain class. Allows runtime selection of lambda update algorithm.
 */

namespace spinconstrain
{

/**
 * @brief Set the lambda update strategy type.
 *
 * @details Creates the appropriate strategy object based on the enum value.
 * For BFGS (default), sets strategy_ = nullptr (uses hard-coded lambda_loop.cpp).
 *
 * @param type Strategy type to use
 */
template <typename TK>
void SpinConstrain<TK>::set_strategy_type(LambdaStrategyType type)
{
    strategy_type_ = type;
    switch(type)
    {
        case LambdaStrategyType::BFGS:
            strategy_ = nullptr; // Use hard-coded BFGS in lambda_loop.cpp
            break;
        case LambdaStrategyType::LinearResponse:
            strategy_ = std::unique_ptr<LambdaUpdateStrategy>(
                new LinearResponseUpdate());
            break;
        case LambdaStrategyType::AugmentedLagrangian:
            strategy_ = std::unique_ptr<LambdaUpdateStrategy>(
                new AugmentedLagrangianUpdate());
            break;
        case LambdaStrategyType::HybridDelayed:
            strategy_ = std::unique_ptr<LambdaUpdateStrategy>(
                new HybridDelayedUpdate());
            break;
        default:
            strategy_ = nullptr;
            strategy_type_ = LambdaStrategyType::BFGS;
            break;
    }
}

/**
 * @brief Configure parameters for the active strategy.
 *
 * @param mu_init Initial penalty parameter (AugmentedLagrangian, HybridDelayed)
 * @param mu_max Maximum penalty parameter
 * @param mu_growth Penalty growth factor
 * @param mix_beta Mixing parameter (LinearResponse)
 * @param sc_scf_thr SCF charge convergence threshold (HybridDelayed)
 */
template <typename TK>
void SpinConstrain<TK>::set_strategy_params(double mu_init, double mu_max,
                                             double mu_growth, double mix_beta,
                                             double sc_scf_thr)
{
    if (!strategy_) return; // BFGS uses hard-coded parameters

    if (strategy_type_ == LambdaStrategyType::LinearResponse)
    {
        if (auto* lr = dynamic_cast<LinearResponseUpdate*>(strategy_.get()))
        {
            // mix_beta is the primary tunable parameter for LinearResponse
            // chi_min, chi_max, lambda_max keep defaults
            *lr = LinearResponseUpdate(0.01, 100.0, mix_beta, 10.0);
        }
    }
    else if (strategy_type_ == LambdaStrategyType::AugmentedLagrangian)
    {
        if (auto* al = dynamic_cast<AugmentedLagrangianUpdate*>(strategy_.get()))
        {
            *al = AugmentedLagrangianUpdate(mu_init, mu_max, mu_growth, 5, 10.0);
        }
    }
    else if (strategy_type_ == LambdaStrategyType::HybridDelayed)
    {
        if (auto* hd = dynamic_cast<HybridDelayedUpdate*>(strategy_.get()))
        {
            *hd = HybridDelayedUpdate(sc_scf_thr, mu_init, mu_max, mu_growth, 5, 10, 10.0);
        }
    }
}

// Explicit template instantiation
template class SpinConstrain<std::complex<double>>;
template class SpinConstrain<double>;

} // namespace spinconstrain
