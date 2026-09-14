#ifndef TD_FIELD_PROFILES_H
#define TD_FIELD_PROFILES_H

#include "td_field.h"

namespace elecstate
{

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

} // namespace elecstate

#endif
