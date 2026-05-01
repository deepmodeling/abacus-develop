#include "spin_constrain.h"

#include "lambda_update_strategies.h"

namespace spinconstrain
{

template <typename TK>
void SpinConstrain<TK>::set_strategy_type(LambdaStrategyType type)
{
    strategy_type_ = type;
    switch(type)
    {
        case LambdaStrategyType::BFGS:
            strategy_ = nullptr;
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

template <typename TK>
void SpinConstrain<TK>::set_strategy_params(double mu_init, double mu_max,
                                             double mu_growth, double mix_beta,
                                             double sc_scf_thr)
{
    if (!strategy_) return;

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
