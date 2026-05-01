#include "lambda_update_strategies.h"
#include <sstream>
#include <cstring>

namespace spinconstrain
{

// ===================================================================
// Helper functions
// ===================================================================

double compute_rms_error(const std::vector<ModuleBase::Vector3<double>>& Mi,
                         const std::vector<ModuleBase::Vector3<double>>& target_mag,
                         const std::vector<ModuleBase::Vector3<int>>& constrain,
                         int nat)
{
    double sum = 0.0;
    int n_count = 0;
    for (int ia = 0; ia < nat; ++ia)
    {
        for (int ic = 0; ic < 3; ++ic)
        {
            if (constrain[ia][ic] != 0)
            {
                double diff = Mi[ia][ic] - target_mag[ia][ic];
                sum += diff * diff;
                ++n_count;
            }
        }
    }
    if (n_count == 0) return 0.0;
    return std::sqrt(sum / n_count);
}

int count_converged(const std::vector<ModuleBase::Vector3<double>>& Mi,
                    const std::vector<ModuleBase::Vector3<double>>& target_mag,
                    const std::vector<ModuleBase::Vector3<int>>& constrain,
                    double sc_thr,
                    int nat)
{
    int count = 0;
    for (int ia = 0; ia < nat; ++ia)
    {
        for (int ic = 0; ic < 3; ++ic)
        {
            if (constrain[ia][ic] != 0)
            {
                double diff = Mi[ia][ic] - target_mag[ia][ic];
                if (std::abs(diff) < sc_thr)
                {
                    ++count;
                }
            }
        }
    }
    return count;
}

void cap_lambda(std::vector<ModuleBase::Vector3<double>>& lambda,
                const std::vector<ModuleBase::Vector3<int>>& constrain,
                double lambda_max,
                int nat)
{
    for (int ia = 0; ia < nat; ++ia)
    {
        for (int ic = 0; ic < 3; ++ic)
        {
            if (constrain[ia][ic] != 0)
            {
                if (lambda[ia][ic] > lambda_max) lambda[ia][ic] = lambda_max;
                if (lambda[ia][ic] < -lambda_max) lambda[ia][ic] = -lambda_max;
            }
        }
    }
}

// ===================================================================
// Scheme B: Linear Response (One-Step) Update
// ===================================================================

LinearResponseUpdate::LinearResponseUpdate(double chi_min,
                                           double chi_max,
                                           double mix_beta,
                                           double lambda_max)
    : chi_min_(chi_min), chi_max_(chi_max), mix_beta_(mix_beta),
      lambda_max_(lambda_max), converged_(false), last_rms_(1e30)
{
}

LambdaUpdateResult LinearResponseUpdate::update_lambda(
    std::vector<ModuleBase::Vector3<double>>& lambda,
    const std::vector<ModuleBase::Vector3<double>>& Mi,
    const std::vector<ModuleBase::Vector3<double>>& target_mag,
    const std::vector<ModuleBase::Vector3<int>>& constrain,
    double sc_thr,
    int iter,
    int nat)
{
    LambdaUpdateResult result;
    result.n_atoms = nat;

    // Ensure response matrix is properly sized
    if (static_cast<int>(chi_.size()) != nat)
    {
        chi_.assign(nat, ModuleBase::Vector3<double>(1.0, 1.0, 1.0));
    }

    // Estimate chi from history if we have enough iterations
    if (iter >= 2 && static_cast<int>(Mi_history_.size()) >= 2)
    {
        const std::vector<ModuleBase::Vector3<double>>& Mi_old = Mi_history_[Mi_history_.size() - 2];
        const std::vector<ModuleBase::Vector3<double>>& lambda_old = lambda_history_[lambda_history_.size() - 2];
        for (int ia = 0; ia < nat; ++ia)
        {
            for (int ic = 0; ic < 3; ++ic)
            {
                if (constrain[ia][ic] == 0) continue;
                double dlambda = lambda[ia][ic] - lambda_old[ia][ic];
                double dM = Mi[ia][ic] - Mi_old[ia][ic];
                if (std::abs(dlambda) > 1e-8)
                {
                    double chi_new = dM / dlambda;
                    if (chi_new > chi_min_ && chi_new < chi_max_)
                    {
                        chi_[ia][ic] = chi_new;
                    }
                }
            }
        }
    }

    // Update lambda: lambda += mix_beta * (M_target - M) / chi
    for (int ia = 0; ia < nat; ++ia)
    {
        for (int ic = 0; ic < 3; ++ic)
        {
            if (constrain[ia][ic] == 0) continue;
            double residual = target_mag[ia][ic] - Mi[ia][ic];
            double delta = residual / chi_[ia][ic];
            lambda[ia][ic] += mix_beta_ * delta;
        }
    }

    // Cap lambda
    cap_lambda(lambda, constrain, lambda_max_, nat);

    // Save history
    Mi_history_.push_back(Mi);
    lambda_history_.push_back(lambda);
    // Keep only last 5 entries
    if (static_cast<int>(Mi_history_.size()) > 5)
    {
        Mi_history_.erase(Mi_history_.begin());
        lambda_history_.erase(lambda_history_.begin());
    }

    // Compute result
    result.rms_error = compute_rms_error(Mi, target_mag, constrain, nat);
    result.n_converged = count_converged(Mi, target_mag, constrain, sc_thr, nat);

    double max_l = 0.0;
    for (int ia = 0; ia < nat; ++ia)
    {
        for (int ic = 0; ic < 3; ++ic)
        {
            if (constrain[ia][ic] != 0)
            {
                max_l = std::max(max_l, std::abs(lambda[ia][ic]));
            }
        }
    }
    result.max_lambda = max_l;

    converged_ = (result.rms_error < sc_thr);
    result.status = converged_ ? "converged" : "updating";

    return result;
}

// ===================================================================
// Scheme C: Augmented Lagrangian Update
// ===================================================================

AugmentedLagrangianUpdate::AugmentedLagrangianUpdate(double mu_init,
                                                      double mu_max,
                                                      double mu_growth,
                                                      int mu_update_interval,
                                                      double lambda_max)
    : mu_(mu_init), mu_init_(mu_init), mu_max_(mu_max),
      mu_growth_(mu_growth), mu_update_interval_(mu_update_interval),
      lambda_max_(lambda_max), converged_(false), last_iter_(0)
{
}

LambdaUpdateResult AugmentedLagrangianUpdate::update_lambda(
    std::vector<ModuleBase::Vector3<double>>& lambda,
    const std::vector<ModuleBase::Vector3<double>>& Mi,
    const std::vector<ModuleBase::Vector3<double>>& target_mag,
    const std::vector<ModuleBase::Vector3<int>>& constrain,
    double sc_thr,
    int iter,
    int nat)
{
    LambdaUpdateResult result;
    result.n_atoms = nat;
    last_iter_ = iter;

    // Dual variable update: lambda += mu * (M - M_target)
    for (int ia = 0; ia < nat; ++ia)
    {
        for (int ic = 0; ic < 3; ++ic)
        {
            if (constrain[ia][ic] == 0) continue;
            double violation = Mi[ia][ic] - target_mag[ia][ic];
            lambda[ia][ic] += mu_ * violation;
        }
    }

    // Cap lambda
    cap_lambda(lambda, constrain, lambda_max_, nat);

    // Grow mu periodically
    if (iter > 0 && iter % mu_update_interval_ == 0)
    {
        mu_ = std::min(mu_max_, mu_ * mu_growth_);
    }

    // Compute result
    result.rms_error = compute_rms_error(Mi, target_mag, constrain, nat);
    result.n_converged = count_converged(Mi, target_mag, constrain, sc_thr, nat);

    double max_l = 0.0;
    for (int ia = 0; ia < nat; ++ia)
    {
        for (int ic = 0; ic < 3; ++ic)
        {
            if (constrain[ia][ic] != 0)
            {
                max_l = std::max(max_l, std::abs(lambda[ia][ic]));
            }
        }
    }
    result.max_lambda = max_l;

    converged_ = (result.rms_error < sc_thr);
    result.status = converged_ ? "converged" : "updating";

    return result;
}

// ===================================================================
// Scheme D: Hybrid Delayed Update
// ===================================================================

HybridDelayedUpdate::HybridDelayedUpdate(double sc_scf_thr,
                                          double mu_init,
                                          double mu_max,
                                          double mu_growth,
                                          int mu_update_interval,
                                          int max_inner_steps,
                                          double lambda_max)
    : sc_scf_thr_(sc_scf_thr), drho_(1e30), mu_(mu_init), mu_init_(mu_init),
      mu_max_(mu_max), mu_growth_(mu_growth),
      mu_update_interval_(mu_update_interval),
      max_inner_steps_(max_inner_steps), lambda_max_(lambda_max),
      converged_(false), inner_steps_(0), phase_("early")
{
}

LambdaUpdateResult HybridDelayedUpdate::update_lambda(
    std::vector<ModuleBase::Vector3<double>>& lambda,
    const std::vector<ModuleBase::Vector3<double>>& Mi,
    const std::vector<ModuleBase::Vector3<double>>& target_mag,
    const std::vector<ModuleBase::Vector3<int>>& constrain,
    double sc_thr,
    int iter,
    int nat)
{
    LambdaUpdateResult result;
    result.n_atoms = nat;

    // Phase decision
    if (drho_ > sc_scf_thr_ * 100)
    {
        // Early phase: skip lambda update
        phase_ = "early";
        result.rms_error = compute_rms_error(Mi, target_mag, constrain, nat);
        result.n_converged = 0;
        result.max_lambda = 0.0;
        for (int ia = 0; ia < nat; ++ia)
        {
            for (int ic = 0; ic < 3; ++ic)
            {
                if (constrain[ia][ic] != 0)
                {
                    result.max_lambda = std::max(result.max_lambda, std::abs(lambda[ia][ic]));
                }
            }
        }
        converged_ = (result.rms_error < sc_thr);
        result.status = "skipped_early";
        return result;
    }
    else if (drho_ > sc_scf_thr_)
    {
        // Mid phase: Augmented Lagrangian lightweight update
        phase_ = "mid";
        for (int ia = 0; ia < nat; ++ia)
        {
            for (int ic = 0; ic < 3; ++ic)
            {
                if (constrain[ia][ic] == 0) continue;
                double violation = Mi[ia][ic] - target_mag[ia][ic];
                lambda[ia][ic] += mu_ * violation;
            }
        }
        cap_lambda(lambda, constrain, lambda_max_, nat);

        if (iter > 0 && iter % mu_update_interval_ == 0)
        {
            mu_ = std::min(mu_max_, mu_ * mu_growth_);
        }
    }
    else
    {
        // Late phase: Augmented Lagrangian + inner loop fallback
        phase_ = "late";
        for (int ia = 0; ia < nat; ++ia)
        {
            for (int ic = 0; ic < 3; ++ic)
            {
                if (constrain[ia][ic] == 0) continue;
                double violation = Mi[ia][ic] - target_mag[ia][ic];
                lambda[ia][ic] += mu_ * violation;
            }
        }
        cap_lambda(lambda, constrain, lambda_max_, nat);

        if (iter > 0 && iter % mu_update_interval_ == 0)
        {
            mu_ = std::min(mu_max_, mu_ * mu_growth_);
        }

        // Check if fallback to inner loop is needed
        double rms = compute_rms_error(Mi, target_mag, constrain, nat);
        if (rms > sc_thr * 10 && inner_steps_ < max_inner_steps_)
        {
            result.status = "fallback_triggered";
            inner_steps_++;
        }
    }

    // Compute result
    result.rms_error = compute_rms_error(Mi, target_mag, constrain, nat);
    result.n_converged = count_converged(Mi, target_mag, constrain, sc_thr, nat);

    double max_l = 0.0;
    for (int ia = 0; ia < nat; ++ia)
    {
        for (int ic = 0; ic < 3; ++ic)
        {
            if (constrain[ia][ic] != 0)
            {
                max_l = std::max(max_l, std::abs(lambda[ia][ic]));
            }
        }
    }
    result.max_lambda = max_l;

    converged_ = (result.rms_error < sc_thr);
    if (result.status != "fallback_triggered")
    {
        if (converged_)
        {
            result.status = "converged";
        }
        else
        {
            result.status = std::string("updating_") + phase_;
        }
    }

    return result;
}

} // namespace spinconstrain
