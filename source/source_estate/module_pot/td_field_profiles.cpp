#include "td_field_profiles.h"

#include "source_base/constants.h"

#include <cmath>

namespace elecstate
{

TDGaussianProfile::TDGaussianProfile(const double omega,
                                     const double phase,
                                     const double sigma,
                                     const double center_step,
                                     const double amplitude,
                                     const double dt)
    : omega_(omega), phase_(phase), sigma_(sigma), center_step_(center_step), amplitude_(amplitude), dt_(dt)
{
}

double TDGaussianProfile::electric_field(const TDFieldSample& sample) const
{
    const double relative_time = (sample.step_position() - center_step_) * dt_;
    return std::cos(omega_ * relative_time + phase_) * std::exp(-0.5 * relative_time * relative_time / (sigma_ * sigma_)) * amplitude_;
}

TDTrapezoidProfile::TDTrapezoidProfile(const double omega,
                                       const double phase,
                                       const double rise_end_step,
                                       const double plateau_end_step,
                                       const double fall_end_step,
                                       const double amplitude)
    : omega_(omega), phase_(phase), rise_end_step_(rise_end_step), plateau_end_step_(plateau_end_step), fall_end_step_(fall_end_step),
      amplitude_(amplitude)
{
}

double TDTrapezoidProfile::electric_field(const TDFieldSample& sample) const
{
    double envelope = 0.0;
    // Segment selection follows the electronic step, while step_position()
    // resolves the Simpson nodes used inside a selected linear segment.
    if (sample.electronic_step < rise_end_step_)
    {
        envelope = sample.step_position() / rise_end_step_;
    }
    else if (sample.electronic_step < plateau_end_step_)
    {
        envelope = 1.0;
    }
    else if (sample.electronic_step < fall_end_step_)
    {
        envelope = (fall_end_step_ - sample.step_position()) / (fall_end_step_ - plateau_end_step_);
    }
    return envelope * amplitude_ * std::cos(omega_ * sample.time + phase_);
}

TDTrigonometricProfile::TDTrigonometricProfile(const double omega1,
                                               const double omega2,
                                               const double phase1,
                                               const double phase2,
                                               const double amplitude)
    : omega1_(omega1), omega2_(omega2), phase1_(phase1), phase2_(phase2), amplitude_(amplitude)
{
}

double TDTrigonometricProfile::electric_field(const TDFieldSample& sample) const
{
    const double envelope = std::sin(omega2_ * sample.time + phase2_);
    return amplitude_ * std::cos(omega1_ * sample.time + phase1_) * envelope * envelope;
}

TDHeavisideProfile::TDHeavisideProfile(const double switch_step, const double amplitude) : switch_step_(switch_step), amplitude_(amplitude)
{
}

double TDHeavisideProfile::electric_field(const TDFieldSample& sample) const
{
    return sample.electronic_step < switch_step_ ? amplitude_ : 0.0;
}

TDSupersineProfile::TDSupersineProfile(const double omega,
                                       const double phase,
                                       const double sigma,
                                       const int start_step,
                                       const int end_step,
                                       const double amplitude,
                                       const double dt)
    : omega_(omega), phase_(phase), sigma_(sigma), start_step_(start_step), end_step_(end_step), amplitude_(amplitude), dt_(dt)
{
}

double TDSupersineProfile::electric_field(const TDFieldSample& sample) const
{
    // Integer substep coordinates make both pulse boundaries exactly zero and
    // avoid floating-point comparisons of independently constructed times.
    const long long start_substep = static_cast<long long>(start_step_) * sample.subdivisions;
    const long long end_substep = static_cast<long long>(end_step_) * sample.subdivisions;
    const long long sample_substep = static_cast<long long>(sample.electronic_step) * sample.subdivisions + sample.simpson_node;
    const long long local_substep = sample_substep - start_substep;
    const long long duration_substeps = end_substep - start_substep;
    const double duration = (end_step_ - start_step_) * dt_;

    double envelope = 0.0;
    double envelope_derivative = 0.0;
    if (local_substep > 0 && local_substep < duration_substeps)
    {
        const double x = static_cast<double>(local_substep) / duration_substeps;
        const double center_offset = x - 0.5;
        if (std::abs(center_offset) <= 1.0e-14)
        {
            // Symmetry fixes the envelope and its derivative at pulse center.
            envelope = 1.0;
        }
        else
        {
            const double sine = std::sin(ModuleBase::PI * x);
            const double log_sine = std::log(sine);
            const double absolute_offset = std::abs(center_offset);
            // Evaluate the variable-power envelope in logarithmic form for
            // better behavior near the compact-support boundaries.
            envelope = std::exp(ModuleBase::PI * absolute_offset * log_sine / sigma_);
            if (envelope > 0.0)
            {
                // Differentiate log(envelope) first, then use f' = f (log f)'.
                const double logarithmic_derivative
                    = ModuleBase::PI * std::copysign(1.0, center_offset) * log_sine / sigma_
                      + ModuleBase::PI * ModuleBase::PI * absolute_offset * std::cos(ModuleBase::PI * x) / (sigma_ * sine);
                envelope_derivative = envelope * logarithmic_derivative / duration;
            }
        }
    }

    const double time_from_center = (static_cast<double>(local_substep) / sample.subdivisions - 0.5 * (end_step_ - start_step_)) * dt_;
    const double carrier_phase = omega_ * time_from_center + phase_;
    // The derivative term is required so that this field is exactly -dA/dt
    // for the analytic compact-support supersine vector potential.
    return amplitude_ * (envelope * std::cos(carrier_phase) + envelope_derivative * std::sin(carrier_phase) / omega_);
}

} // namespace elecstate
