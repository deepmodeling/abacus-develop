#ifndef TD_FIELD_H
#define TD_FIELD_H

#include "source_base/vector3.h"

#include <memory>
#include <string>
#include <vector>

struct Input_para;

namespace elecstate
{

/**
 * @brief Sampling location used to evaluate one time-dependent field profile.
 *
 * A vector-gauge electronic step is divided into Simpson subintervals. The
 * sample represents node `simpson_node` in electronic step `electronic_step`;
 * length-gauge evaluation uses node zero directly.
 */
struct TDFieldSample
{
    /**
     * @brief Construct a field sampling point.
     *
     * @param electronic_step_in Zero-based electronic-step index.
     * @param simpson_node_in Node index in the current Simpson interval.
     * @param subdivisions_in Number of subintervals in one electronic step.
     * @param time_in Physical sampling time in internal atomic time units.
     */
    TDFieldSample(const int electronic_step_in, const int simpson_node_in, const int subdivisions_in, const double time_in)
        : electronic_step(electronic_step_in), simpson_node(simpson_node_in), subdivisions(subdivisions_in), time(time_in)
    {
    }

    /**
     * @brief Return the sampling location in continuous electronic-step units.
     */
    double step_position() const;

    int electronic_step; ///< Zero-based electronic-step index.
    int simpson_node;    ///< Node index within the current electronic step.
    int subdivisions;    ///< Number of subintervals in one electronic step.
    double time;         ///< Sampling time in internal atomic time units.
};

/**
 * @brief Scalar time profile of one configured electric field.
 *
 * Implementations return the field in ABACUS internal propagation units. A
 * profile contains no Cartesian direction; direction handling is owned by
 * TDField and TDFieldManager.
 */
class TDFieldProfile
{
  public:
    /** @brief Destroy the polymorphic field profile. */
    virtual ~TDFieldProfile() = default;

    /**
     * @brief Evaluate the electric field at one sampling point.
     *
     * @param sample Electronic-step and Simpson-node sampling information.
     * @return Scalar field value in internal propagation units.
     */
    virtual double electric_field(const TDFieldSample& sample) const = 0;
};

/** @brief Carrier cosine multiplied by a Gaussian envelope. */
class TDGaussianProfile : public TDFieldProfile
{
  public:
    /**
     * @brief Construct a Gaussian field profile in internal units.
     *
     * @param omega Angular frequency.
     * @param phase Carrier phase in radians.
     * @param sigma Gaussian width in atomic time units.
     * @param center_step Pulse center in electronic-step units.
     * @param amplitude Peak field amplitude.
     * @param dt Electronic time step in atomic time units.
     */
    TDGaussianProfile(double omega, double phase, double sigma, double center_step, double amplitude, double dt);

    /** @copydoc TDFieldProfile::electric_field */
    double electric_field(const TDFieldSample& sample) const override;

  private:
    double omega_;
    double phase_;
    double sigma_;
    double center_step_;
    double amplitude_;
    double dt_;
};

/** @brief Carrier cosine with linear-rise, plateau, and linear-fall envelope. */
class TDTrapezoidProfile : public TDFieldProfile
{
  public:
    /**
     * @brief Construct a trapezoidal-envelope field profile.
     *
     * @param omega Angular frequency in internal units.
     * @param phase Carrier phase in radians.
     * @param rise_end_step End of the linear rise in electronic-step units.
     * @param plateau_end_step End of the plateau in electronic-step units.
     * @param fall_end_step End of the linear fall in electronic-step units.
     * @param amplitude Peak field amplitude in internal units.
     */
    TDTrapezoidProfile(double omega, double phase, double rise_end_step, double plateau_end_step, double fall_end_step, double amplitude);

    /** @copydoc TDFieldProfile::electric_field */
    double electric_field(const TDFieldSample& sample) const override;

  private:
    double omega_;
    double phase_;
    double rise_end_step_;
    double plateau_end_step_;
    double fall_end_step_;
    double amplitude_;
};

/** @brief Cosine carrier multiplied by a squared-sine envelope. */
class TDTrigonometricProfile : public TDFieldProfile
{
  public:
    /**
     * @brief Construct a trigonometric field profile.
     *
     * @param omega1 Carrier angular frequency in internal units.
     * @param omega2 Envelope angular frequency in internal units.
     * @param phase1 Carrier phase in radians.
     * @param phase2 Envelope phase in radians.
     * @param amplitude Peak field amplitude in internal units.
     */
    TDTrigonometricProfile(double omega1, double omega2, double phase1, double phase2, double amplitude);

    /** @copydoc TDFieldProfile::electric_field */
    double electric_field(const TDFieldSample& sample) const override;

  private:
    double omega1_;
    double omega2_;
    double phase1_;
    double phase2_;
    double amplitude_;
};

/** @brief Constant field switched off at a specified electronic step. */
class TDHeavisideProfile : public TDFieldProfile
{
  public:
    /**
     * @brief Construct a Heaviside field profile.
     *
     * @param switch_step First electronic step for which the field is zero.
     * @param amplitude Field amplitude before the switch.
     */
    TDHeavisideProfile(double switch_step, double amplitude);

    /** @copydoc TDFieldProfile::electric_field */
    double electric_field(const TDFieldSample& sample) const override;

  private:
    double switch_step_;
    double amplitude_;
};

/**
 * @brief Compact supersine pulse derived from an analytic vector potential.
 *
 * The electric field includes both the carrier term and the envelope-derivative
 * term required by the relation between electric and vector potentials.
 */
class TDSupersineProfile : public TDFieldProfile
{
  public:
    /**
     * @brief Construct a supersine field profile.
     *
     * @param omega Carrier angular frequency in internal units.
     * @param phase Carrier phase in radians.
     * @param sigma Supersine envelope-shape parameter.
     * @param start_step First pulse-boundary electronic step.
     * @param end_step Last pulse-boundary electronic step.
     * @param amplitude Field scale in internal units.
     * @param dt Electronic time step in atomic time units.
     */
    TDSupersineProfile(double omega, double phase, double sigma, int start_step, int end_step, double amplitude, double dt);

    /** @copydoc TDFieldProfile::electric_field */
    double electric_field(const TDFieldSample& sample) const override;

  private:
    double omega_;
    double phase_;
    double sigma_;
    int start_step_;
    int end_step_;
    double amplitude_;
    double dt_;
};

/**
 * @brief One configured field occurrence and its Cartesian direction.
 *
 * Each occurrence remains independent even when multiple fields share a
 * direction. The class is move-only because it owns a polymorphic profile.
 */
class TDField
{
  public:
    /**
     * @brief Construct one directed time-dependent field.
     *
     * @param direction Zero-based absolute Cartesian direction.
     * @param profile Scalar field profile owned by this object.
     * @param subdivisions Simpson subinterval count for one electronic step.
     */
    TDField(int direction, std::unique_ptr<TDFieldProfile> profile, int subdivisions);

    /** @brief Move-construct a field while transferring profile ownership. */
    TDField(TDField&& other) = default;

    /** @brief Move-assign a field while transferring profile ownership. */
    TDField& operator=(TDField&& other) = default;
    TDField(const TDField&) = delete;
    TDField& operator=(const TDField&) = delete;

    /** @brief Return the zero-based absolute Cartesian direction. */
    int direction() const;

    /** @brief Return the Simpson subinterval count for one electronic step. */
    int subdivisions() const;

    /** @brief Evaluate this field's scalar profile at a sampling point. */
    double electric_field(const TDFieldSample& sample) const;

  private:
    int direction_;
    std::unique_ptr<TDFieldProfile> profile_;
    int subdivisions_;
};

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
