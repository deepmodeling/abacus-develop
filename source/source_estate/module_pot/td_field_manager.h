#ifndef TD_FIELD_MANAGER_H
#define TD_FIELD_MANAGER_H

#include "source_base/vector3.h"
#include "td_field.h"

#include <memory>
#include <string>
#include <vector>

struct Input_para;

namespace elecstate
{

/**
 * @brief Own and advance all time-dependent electric fields in RT-TDDFT.
 *
 * The manager preserves per-occurrence values for output while also summing
 * fields that share a Cartesian direction. It is the single source of step,
 * electric-field, and vector-potential state for all spatial gauges.
 */
class TDFieldManager
{
  public:
    /** @brief Sample E at step*dt without integrating A (length gauge or initial state). */
    void prepare_sample(const int step);

    /** @brief Integrate [left_step*dt, (left_step+1)*dt] and sample E at sample_step*dt.
     * Repeating the same interval and sample is idempotent. Propagation A is the
     * endpoint average, not an exact midpoint sample.
     */
    void prepare_interval(const int left_step, const int sample_step);

    /** @brief Supply file-based propagation samples before preparing any steps. */
    void set_A_samples(const std::vector<ModuleBase::Vector3<double>>& samples_ha);

    /** @brief Return integrated endpoints in Hartree units; unavailable for file-based propagation samples. */
    const ModuleBase::Vector3<double>& A_left_ha() const;
    const ModuleBase::Vector3<double>& A_right_ha() const;
    /** @brief Return the endpoint average or supplied file sample used for propagation, in Hartree units. */
    const ModuleBase::Vector3<double>& A_prop_ha() const;

    /**
     * @brief Restore the electronic step and vector-potential state.
     *
     * @param file_dir Directory containing `Restart_td.txt`.
     */
    void read_restart(const std::string& file_dir);
    /** @brief Save three annotated data rows describing the interval state in Hartree units. */
    void write_restart(const std::string& file_dir) const;

    /** @brief Return the spatial-gauge selector supplied by `td_stype`. */
    int gauge() const;

    /** @brief Return the current zero-based electronic-step index. */
    int current_step() const;

    /** @brief Return the electronic time step in Hartree atomic time units. */
    double dt_ha() const;

    /** @brief Return the first reduced-coordinate cut of the length gauge. */
    double length_cut1() const;

    /** @brief Return the second reduced-coordinate cut of the length gauge. */
    double length_cut2() const;

    /** @brief Return whether the configured field is active at this step. */
    bool active() const;

    /** @brief Return all configured fields in input-occurrence order. */
    const std::vector<TDField>& fields() const;

    /** @brief Return per-occurrence field samples for the current step. */
    const std::vector<double>& field_vals_ha() const;

    /** @brief Return the direction-summed sampled electric field in Hartree atomic units. */
    const ModuleBase::Vector3<double>& efield_ha() const;

  private:
    TDFieldManager(bool enabled,
                   int gauge,
                   int start_step,
                   int end_step,
                   double dt,
                   double length_cut1,
                   double length_cut2,
                   std::vector<TDField> fields);

    bool enabled_;
    int gauge_;
    int start_step_;
    int end_step_;
    double dt_ha_;
    double length_cut1_;
    double length_cut2_;
    std::vector<TDField> fields_;
    int current_step_;
    bool active_;
    std::vector<double> field_vals_ha_;
    int interval_left_ = -1;
    std::vector<ModuleBase::Vector3<double>> A_samples_ha_;
    ModuleBase::Vector3<double> A_left_ha_;
    ModuleBase::Vector3<double> A_right_ha_;
    ModuleBase::Vector3<double> A_prop_ha_;
    ModuleBase::Vector3<double> efield_ha_;
    void sample_field(const int step);
    void select_A_prop();

    friend std::shared_ptr<TDFieldManager> create_td_field_manager(const Input_para& input);
};

/**
 * @brief Build all field profiles and convert user input to propagation units.
 *
 * Parameters specific to a waveform are paired with `td_ttype` by occurrence,
 * while `td_vext_dire` is paired by the overall field index.
 *
 * @param input Validated ABACUS input parameters.
 * @return Shared manager used by the ESolver and time-dependent potential.
 */
std::shared_ptr<TDFieldManager> create_td_field_manager(const Input_para& input);

} // namespace elecstate

#endif
