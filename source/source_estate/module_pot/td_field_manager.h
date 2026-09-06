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
    /**
     * @brief Advance one length-gauge step and sample every field at its start.
     */
    void advance_length_gauge();

    /**
     * @brief Advance one velocity/hybrid-gauge step using Simpson integration.
     */
    void advance_vector_gauge();

    /**
     * @brief Restore the electronic step and vector-potential state.
     *
     * @param file_dir Directory containing `Restart_td.txt`.
     */
    void read_restart(const std::string& file_dir);

    /** @brief Return the spatial-gauge selector supplied by `td_stype`. */
    int gauge() const;

    /** @brief Return the current zero-based electronic-step index. */
    int current_step() const;

    /** @brief Return the electronic time step in internal atomic time units. */
    double dt() const;

    /** @brief Return the first reduced-coordinate cut of the length gauge. */
    double length_cut1() const;

    /** @brief Return the second reduced-coordinate cut of the length gauge. */
    double length_cut2() const;

    /** @brief Return whether the configured field is active at this step. */
    bool active() const;

    /** @brief Return all configured fields in input-occurrence order. */
    const std::vector<TDField>& fields() const;

    /** @brief Return per-occurrence field samples for the current step. */
    const std::vector<double>& field_values() const;

    /** @brief Return the midpoint vector potential in propagation units. */
    const ModuleBase::Vector3<double>& vector_potential() const;

    /** @brief Return the integrated vector-potential change for this step. */
    const ModuleBase::Vector3<double>& vector_potential_laststep() const;

    /** @brief Return the direction-summed instantaneous hybrid-gauge field. */
    const ModuleBase::Vector3<double>& electric_field() const;

    /** @brief Return the direction-summed field used by force evaluation. */
    const ModuleBase::Vector3<double>& total_electric_field() const;

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
    double dt_;
    double length_cut1_;
    double length_cut2_;
    std::vector<TDField> fields_;
    int current_step_;
    bool active_;
    std::vector<double> field_values_;
    ModuleBase::Vector3<double> vector_potential_;
    ModuleBase::Vector3<double> vector_potential_laststep_;
    ModuleBase::Vector3<double> electric_field_;
    ModuleBase::Vector3<double> total_electric_field_;

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
