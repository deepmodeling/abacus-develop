#ifndef ABACUS_SOURCE_DFTB_CHARGE_MIXING_H
#define ABACUS_SOURCE_DFTB_CHARGE_MIXING_H

#include <string>
#include <vector>

namespace ModuleDFTB
{

struct DftbChargeMixerParameters
{
    std::string method = "linear";
    double mixing_parameter = 0.2;
    int history = 6;
    // Johnson modified-Broyden controls, matching DFTB+ defaults.
    double inverse_jacobi_weight = 0.01;
    double minimal_weight = 1.0;
    double maximal_weight = 1.0e5;
    double weight_factor = 1.0e-2;
};

/**
 * Stateful SCC charge mixer for linear, Pulay/DIIS, and Johnson modified-Broyden mixing.
 * Residuals use the fixed-point convention q_out - q_in; each returned vector is
 * projected back to the requested total charge.
 */
class DftbChargeMixer
{
  public:
    explicit DftbChargeMixer(const DftbChargeMixerParameters& parameters);

    std::vector<double> mix(const std::vector<double>& charges,
                            const std::vector<double>& residual,
                            double target_total_charge);
    const std::string& last_step() const { return last_step_; }

  private:
    DftbChargeMixerParameters parameters_;
    std::string last_step_;
    std::vector<std::vector<double>> charge_history_;
    std::vector<std::vector<double>> residual_history_;
    std::vector<double> previous_charges_;
    std::vector<double> previous_residual_;
    std::vector<std::vector<double>> normalized_residual_differences_;
    std::vector<std::vector<double>> update_differences_;
    std::vector<double> broyden_weights_;
    bool has_previous_ = false;
};

} // namespace ModuleDFTB

#endif // ABACUS_SOURCE_DFTB_CHARGE_MIXING_H
