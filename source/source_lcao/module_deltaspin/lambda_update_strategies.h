#ifndef LAMBDA_UPDATE_STRATEGIES_H
#define LAMBDA_UPDATE_STRATEGIES_H

#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <limits>

#include "source_base/vector3.h"

namespace spinconstrain
{

/**
 * @brief Result struct for lambda update operations
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
 * @brief Pure abstract base class for lambda update strategies
 */
class LambdaUpdateStrategy
{
  public:
    virtual ~LambdaUpdateStrategy() = default;

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
 * @brief Compute RMS error of |M - M_target| (respecting constrain flags)
 */
double compute_rms_error(const std::vector<ModuleBase::Vector3<double>>& Mi,
                         const std::vector<ModuleBase::Vector3<double>>& target_mag,
                         const std::vector<ModuleBase::Vector3<int>>& constrain,
                         int nat);

/**
 * @brief Count converged components
 */
int count_converged(const std::vector<ModuleBase::Vector3<double>>& Mi,
                    const std::vector<ModuleBase::Vector3<double>>& target_mag,
                    const std::vector<ModuleBase::Vector3<int>>& constrain,
                    double sc_thr,
                    int nat);

/**
 * @brief Apply absolute cap to lambda values
 */
void cap_lambda(std::vector<ModuleBase::Vector3<double>>& lambda,
                const std::vector<ModuleBase::Vector3<int>>& constrain,
                double lambda_max,
                int nat);

// ===================================================================
// Scheme B: Linear Response (One-Step) Update
// ===================================================================

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
    std::vector<ModuleBase::Vector3<double>> chi_;
    std::vector<std::vector<ModuleBase::Vector3<double>>> Mi_history_;
    std::vector<std::vector<ModuleBase::Vector3<double>>> lambda_history_;
};

// ===================================================================
// Scheme C: Augmented Lagrangian Update
// ===================================================================

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
    double mu_;
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
    double drho_;
    double mu_;
    double mu_init_;
    double mu_max_;
    double mu_growth_;
    int mu_update_interval_;
    int max_inner_steps_;
    double lambda_max_;
    bool converged_;
    int inner_steps_;
    std::string phase_;
};

} // namespace spinconstrain

#endif // LAMBDA_UPDATE_STRATEGIES_H
