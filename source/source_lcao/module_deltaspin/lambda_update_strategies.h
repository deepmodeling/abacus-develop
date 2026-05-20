#ifndef LAMBDA_UPDATE_STRATEGIES_H
#define LAMBDA_UPDATE_STRATEGIES_H

#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <limits>

#include "source_base/vector3.h"

/**
 * @file lambda_update_strategies.h
 * @brief Alternative lambda update strategies (NOT currently used in production).
 *
 * @par Status
 * These strategies are experimental and NOT compiled into the library
 * (not listed in CMakeLists.txt). The production code uses the hard-coded
 * BFGS-like optimizer in lambda_loop.cpp. These strategies are provided
 * for future development and experimentation.
 *
 * @par Available strategies
 * - BFGS (default): Hard-coded in lambda_loop.cpp, uses conjugate gradient
 * - LinearResponse (Scheme B): Estimates susceptibility chi = dM/dlambda
 * - AugmentedLagrangian (Scheme C): Dual ascent with penalty method
 * - HybridDelayed (Scheme D): Three-phase approach based on charge convergence
 *
 * @par Integration
 * To use these strategies, the following members need to be added to
 * SpinConstrain (in spin_constrain.h):
 *   LambdaStrategyType strategy_type_;
 *   std::unique_ptr<LambdaUpdateStrategy> strategy_;
 * And set_strategy_type()/set_strategy_params() need to be called from
 * init_sc() or the ESolver layer.
 */
namespace spinconstrain
{

/**
 * @brief Result struct for lambda update operations.
 */
struct LambdaUpdateResult
{
    int n_atoms;
    double rms_error;            ///< RMS of |M - M_target| after update
    double max_lambda;           ///< max |lambda| across all atoms/components
    int n_converged;             ///< number of (atom, component) pairs converged
    std::string status;          ///< "converged", "updating", "fallback_triggered"
};

/**
 * @brief Pure abstract base class for lambda update strategies.
 *
 * @par Design pattern
 * Strategy pattern: different update algorithms can be swapped at runtime
 * by creating a concrete subclass and passing it to SpinConstrain.
 */
class LambdaUpdateStrategy
{
  public:
    virtual ~LambdaUpdateStrategy() = default;

    /**
     * @brief Update lambda values based on current magnetic moments.
     *
     * @param lambda Current Lagrange multipliers (modified in place)
     * @param Mi Current magnetic moments
     * @param target_mag Target magnetic moments
     * @param constrain Constraint flags (0=free, 1=constrained)
     * @param sc_thr Convergence threshold
     * @param iter Current iteration number
     * @param nat Number of atoms
     * @return Result struct with convergence info
     */
    virtual LambdaUpdateResult update_lambda(std::vector<ModuleBase::Vector3<double>>& lambda,
                                             const std::vector<ModuleBase::Vector3<double>>& Mi,
                                             const std::vector<ModuleBase::Vector3<double>>& target_mag,
                                             const std::vector<ModuleBase::Vector3<int>>& constrain,
                                             double sc_thr,
                                             int iter,
                                             int nat) = 0;

    virtual std::string name() const = 0;
    virtual bool is_converged() const = 0;
};

/**
 * @brief Compute RMS error of |M - M_target| (respecting constrain flags).
 */
double compute_rms_error(const std::vector<ModuleBase::Vector3<double>>& Mi,
                         const std::vector<ModuleBase::Vector3<double>>& target_mag,
                         const std::vector<ModuleBase::Vector3<int>>& constrain,
                         int nat);

/**
 * @brief Count converged components (|Mi - M_target| < sc_thr).
 */
int count_converged(const std::vector<ModuleBase::Vector3<double>>& Mi,
                    const std::vector<ModuleBase::Vector3<double>>& target_mag,
                    const std::vector<ModuleBase::Vector3<int>>& constrain,
                    double sc_thr,
                    int nat);

/**
 * @brief Apply absolute cap to lambda values to prevent divergence.
 */
void cap_lambda(std::vector<ModuleBase::Vector3<double>>& lambda,
                const std::vector<ModuleBase::Vector3<int>>& constrain,
                double lambda_max,
                int nat);

// ===================================================================
// Scheme B: Linear Response (One-Step) Update
// ===================================================================

/**
 * @brief Linear response lambda update: lambda += mix_beta * (M_target - M) / chi.
 *
 * @par Algorithm
 * Estimates the magnetic susceptibility chi = dM/dlambda from the history
 * of the last 2 iterations:
 *   chi = (Mi_current - Mi_previous) / (lambda_current - lambda_previous)
 * Then updates:
 *   lambda += mix_beta * (M_target - Mi) / chi
 *
 * @par Parameters
 * - chi_min: minimum susceptibility (prevents division by small numbers)
 * - chi_max: maximum susceptibility (prevents unstable updates)
 * - mix_beta: mixing parameter (0.3 = conservative, 1.0 = aggressive)
 * - lambda_max: absolute cap on lambda values
 */
class LinearResponseUpdate : public LambdaUpdateStrategy
{
  public:
    LinearResponseUpdate(double chi_min = 0.01,
                         double chi_max = 100.0,
                         double mix_beta = 0.3,
                         double lambda_max = 10.0);

    LambdaUpdateResult update_lambda(std::vector<ModuleBase::Vector3<double>>& lambda,
                                     const std::vector<ModuleBase::Vector3<double>>& Mi,
                                     const std::vector<ModuleBase::Vector3<double>>& target_mag,
                                     const std::vector<ModuleBase::Vector3<int>>& constrain,
                                     double sc_thr,
                                     int iter,
                                     int nat) override;

    std::string name() const override { return "LinearResponse"; }
    bool is_converged() const override { return converged_; }

    const std::vector<ModuleBase::Vector3<double>>& get_chi() const { return chi_; }

  private:
    double chi_min_;
    double chi_max_;
    double mix_beta_;
    double lambda_max_;
    bool converged_;
    double last_rms_;
    std::vector<ModuleBase::Vector3<double>> chi_; ///< Estimated susceptibility dM/dlambda
    std::vector<std::vector<ModuleBase::Vector3<double>>> Mi_history_; ///< Last 5 Mi values
    std::vector<std::vector<ModuleBase::Vector3<double>>> lambda_history_; ///< Last 5 lambda values
};

// ===================================================================
// Scheme C: Augmented Lagrangian Update
// ===================================================================

/**
 * @brief Augmented Lagrangian lambda update: lambda += mu * (M - M_target).
 *
 * @par Algorithm
 * Simple dual ascent method. The penalty parameter mu grows periodically
 * to enforce the constraint more strongly:
 *   lambda += mu * (Mi - M_target)
 *   mu *= mu_growth every mu_update_interval iterations
 *
 * @par Parameters
 * - mu_init: initial penalty parameter
 * - mu_max: maximum penalty (prevents numerical instability)
 * - mu_growth: growth factor (1.5 = moderate)
 * - mu_update_interval: iterations between mu updates
 * - lambda_max: absolute cap on lambda values
 */
class AugmentedLagrangianUpdate : public LambdaUpdateStrategy
{
  public:
    AugmentedLagrangianUpdate(double mu_init = 0.1,
                              double mu_max = 10.0,
                              double mu_growth = 1.5,
                              int mu_update_interval = 5,
                              double lambda_max = 10.0);

    LambdaUpdateResult update_lambda(std::vector<ModuleBase::Vector3<double>>& lambda,
                                     const std::vector<ModuleBase::Vector3<double>>& Mi,
                                     const std::vector<ModuleBase::Vector3<double>>& target_mag,
                                     const std::vector<ModuleBase::Vector3<int>>& constrain,
                                     double sc_thr,
                                     int iter,
                                     int nat) override;

    std::string name() const override { return "AugmentedLagrangian"; }
    bool is_converged() const override { return converged_; }

    double get_mu() const { return mu_; }
    void reset_mu() { mu_ = mu_init_; }

  private:
    double mu_; ///< Current penalty parameter
    double mu_init_;
    double mu_max_;
    double mu_growth_;
    int mu_update_interval_;
    double lambda_max_;
    bool converged_;
    int last_iter_;
};

// ===================================================================
// Scheme D: Hybrid Delayed Update
// ===================================================================

/**
 * @brief Hybrid delayed update: three-phase approach based on charge convergence.
 *
 * @par Algorithm
 * Phase decision based on drho (charge density change):
 * - Early phase (drho > sc_scf_thr * 100): Skip lambda update entirely.
 *   The charge density is changing too rapidly for lambda optimization.
 * - Mid phase (sc_scf_thr < drho < sc_scf_thr * 100): Lightweight augmented
 *   Lagrangian update with small mu.
 * - Late phase (drho < sc_scf_thr): Full augmented Lagrangian with fallback
 *   to inner loop if RMS error is still large.
 *
 * @par Parameters
 * - sc_scf_thr: SCF charge convergence threshold (phase decision boundary)
 * - mu_init, mu_max, mu_growth: Augmented Lagrangian parameters
 * - max_inner_steps: maximum fallback inner loop iterations
 * - lambda_max: absolute cap on lambda values
 */
class HybridDelayedUpdate : public LambdaUpdateStrategy
{
  public:
    HybridDelayedUpdate(double sc_scf_thr = 1e-3,
                        double mu_init = 0.1,
                        double mu_max = 10.0,
                        double mu_growth = 1.5,
                        int mu_update_interval = 5,
                        int max_inner_steps = 10,
                        double lambda_max = 10.0);

    void set_drho(double drho) { drho_ = drho; }

    LambdaUpdateResult update_lambda(std::vector<ModuleBase::Vector3<double>>& lambda,
                                     const std::vector<ModuleBase::Vector3<double>>& Mi,
                                     const std::vector<ModuleBase::Vector3<double>>& target_mag,
                                     const std::vector<ModuleBase::Vector3<int>>& constrain,
                                     double sc_thr,
                                     int iter,
                                     int nat) override;

    std::string name() const override { return "HybridDelayed"; }
    bool is_converged() const override { return converged_; }

    std::string get_phase() const { return phase_; }
    void reset() { mu_ = mu_init_; inner_steps_ = 0; phase_ = "early"; }

  private:
    double sc_scf_thr_;
    double drho_; ///< Current charge density change
    double mu_;
    double mu_init_;
    double mu_max_;
    double mu_growth_;
    int mu_update_interval_;
    int max_inner_steps_;
    double lambda_max_;
    bool converged_;
    int inner_steps_; ///< Count of fallback inner loop iterations
    std::string phase_; ///< Current phase: "early", "mid", "late"
};

} // namespace spinconstrain

#endif // LAMBDA_UPDATE_STRATEGIES_H
