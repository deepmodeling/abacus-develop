#ifndef TD_FIELD_H
#define TD_FIELD_H

#include <memory>

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

} // namespace elecstate

#endif
