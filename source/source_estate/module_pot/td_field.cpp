#include "td_field.h"

#include "source_base/constants.h"
#include "source_base/math_integral.h"
#include "source_base/tool_quit.h"
#include "source_io/module_parameter/input_parameter.h"

#include <algorithm>
#include <cmath>
#include <fstream>

namespace
{

int integration_subdivisions(const double omega, const double dt, const int gauge)
{
    // Length gauge samples only the beginning of each electronic step and does
    // not integrate the field in time.
    if (gauge == 0)
    {
        return 1;
    }

    // Preserve the legacy frequency-dependent resolution while enforcing the
    // positive, even number of subintervals required by Simpson integration.
    int subdivisions = static_cast<int>(100.0 * std::abs(omega) * dt / ModuleBase::PI);
    subdivisions += subdivisions % 2 == 0 ? 2 : 1;
    return std::max(2, subdivisions);
}

double angular_frequency(const double frequency)
{
    // User frequencies are supplied in fs^-1; profiles use atomic time.
    return frequency * 2.0 * ModuleBase::PI * ModuleBase::AU_to_FS;
}

double field_amplitude(const double amplitude)
{
    // Convert the user-visible V/Angstrom scale to the propagation field unit.
    return amplitude * ModuleBase::BOHR_TO_A / ModuleBase::Ry_to_eV;
}

} // namespace

namespace elecstate
{

double TDFieldSample::step_position() const
{
    return electronic_step + static_cast<double>(simpson_node) / subdivisions;
}

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

TDField::TDField(const int direction, std::unique_ptr<TDFieldProfile> profile, const int subdivisions)
    : direction_(direction), profile_(std::move(profile)), subdivisions_(subdivisions)
{
}

int TDField::direction() const
{
    return direction_;
}

int TDField::subdivisions() const
{
    return subdivisions_;
}

double TDField::electric_field(const TDFieldSample& sample) const
{
    return profile_->electric_field(sample);
}

TDFieldManager::TDFieldManager(const bool enabled,
                               const int gauge,
                               const int start_step,
                               const int end_step,
                               const double dt,
                               const double length_cut1,
                               const double length_cut2,
                               std::vector<TDField> fields)
    : enabled_(enabled), gauge_(gauge), start_step_(start_step), end_step_(end_step), dt_(dt), length_cut1_(length_cut1),
      length_cut2_(length_cut2), fields_(std::move(fields)), current_step_(-1), active_(false), field_values_(fields_.size(), 0.0)
{
    vector_potential_.set(0.0, 0.0, 0.0);
    vector_potential_laststep_.set(0.0, 0.0, 0.0);
    electric_field_.set(0.0, 0.0, 0.0);
    total_electric_field_.set(0.0, 0.0, 0.0);
}

void TDFieldManager::advance_length_gauge()
{
    ++current_step_;
    active_ = enabled_ && current_step_ >= start_step_ && current_step_ <= end_step_;
    std::fill(field_values_.begin(), field_values_.end(), 0.0);
    total_electric_field_.set(0.0, 0.0, 0.0);
    if (!active_)
    {
        return;
    }

    for (std::size_t index = 0; index < fields_.size(); ++index)
    {
        const TDField& field = fields_[index];
        const TDFieldSample sample(current_step_, 0, field.subdivisions(), current_step_ * dt_);
        // Keep every occurrence for output, but sum repeated directions for
        // the physical length-gauge potential and ionic force.
        field_values_[index] = field.electric_field(sample);
        total_electric_field_[field.direction()] += field_values_[index];
    }
}

void TDFieldManager::advance_vector_gauge()
{
    ++current_step_;
    // Finish the second half of the previous interval before integrating the
    // current interval. vector_potential_ therefore remains a midpoint value.
    vector_potential_ = vector_potential_ + vector_potential_laststep_ / 2.0;
    vector_potential_laststep_.set(0.0, 0.0, 0.0);
    electric_field_.set(0.0, 0.0, 0.0);
    total_electric_field_.set(0.0, 0.0, 0.0);
    std::fill(field_values_.begin(), field_values_.end(), 0.0);
    active_ = enabled_ && current_step_ >= start_step_ && current_step_ <= end_step_;
    if (!active_)
    {
        return;
    }

    for (std::size_t index = 0; index < fields_.size(); ++index)
    {
        const TDField& field = fields_[index];
        const int subdivisions = field.subdivisions();
        const double integration_dt = dt_ / subdivisions;
        std::vector<double> samples(subdivisions + 1, 0.0);
        for (int node = 0; node <= subdivisions; ++node)
        {
            const double time = (current_step_ + static_cast<double>(node) / subdivisions) * dt_;
            samples[node] = field.electric_field(TDFieldSample(current_step_, node, subdivisions, time));
        }

        // Integrate E over [n*dt, (n+1)*dt]. The minus sign implements
        // A(t+dt)-A(t) = -integral E(t') dt'.
        double integral = 0.0;
        ModuleBase::Integral::Simpson_Integral(subdivisions + 1, samples.data(), integration_dt, integral);
        vector_potential_laststep_[field.direction()] -= integral;
        // Output and the hybrid-gauge scalar potential use E at the interval
        // start rather than an average over Simpson nodes.
        field_values_[index] = samples.front();
        if (gauge_ == 2)
        {
            electric_field_[field.direction()] += samples.front();
        }
    }

    // Advance from the interval endpoint to its midpoint representation.
    vector_potential_ = vector_potential_ + vector_potential_laststep_ / 2.0;
    if (gauge_ == 2)
    {
        total_electric_field_ = electric_field_;
    }
}

void TDFieldManager::read_restart(const std::string& file_dir)
{
    std::ifstream file((file_dir + "Restart_td.txt").c_str());
    if (!file)
    {
        ModuleBase::WARNING_QUIT("TDFieldManager::read_restart", "No Restart_td.txt!");
    }

    int restart_step = -1;
    file >> restart_step;
    file >> vector_potential_[0] >> vector_potential_[1] >> vector_potential_[2];
    file >> vector_potential_laststep_[0] >> vector_potential_laststep_[1] >> vector_potential_laststep_[2];
    // Retain the legacy restart-file sign convention expected by the first
    // half-step update in advance_vector_gauge().
    vector_potential_laststep_ = -vector_potential_laststep_;
    current_step_ = restart_step - 1;
}

int TDFieldManager::gauge() const
{
    return gauge_;
}

int TDFieldManager::current_step() const
{
    return current_step_;
}

double TDFieldManager::dt() const
{
    return dt_;
}

double TDFieldManager::length_cut1() const
{
    return length_cut1_;
}

double TDFieldManager::length_cut2() const
{
    return length_cut2_;
}

bool TDFieldManager::active() const
{
    return active_;
}

const std::vector<TDField>& TDFieldManager::fields() const
{
    return fields_;
}

const std::vector<double>& TDFieldManager::field_values() const
{
    return field_values_;
}

const ModuleBase::Vector3<double>& TDFieldManager::vector_potential() const
{
    return vector_potential_;
}

const ModuleBase::Vector3<double>& TDFieldManager::vector_potential_laststep() const
{
    return vector_potential_laststep_;
}

const ModuleBase::Vector3<double>& TDFieldManager::electric_field() const
{
    return electric_field_;
}

const ModuleBase::Vector3<double>& TDFieldManager::total_electric_field() const
{
    return total_electric_field_;
}

std::shared_ptr<TDFieldManager> create_td_field_manager(const Input_para& input)
{
    // An explicitly supplied electronic time step takes precedence; otherwise
    // derive it from the ionic time step and the electronic-step count.
    const double dt
        = input.td_dt != -1.0 ? input.td_dt / ModuleBase::AU_to_FS : input.mdp.md_dt / input.estep_per_md / ModuleBase::AU_to_FS;
    // Each waveform-specific parameter vector is indexed by occurrences of
    // that waveform, not by the overall position in td_ttype.
    std::vector<std::size_t> occurrences(5, 0);
    std::vector<TDField> fields;
    fields.reserve(input.td_ttype.size());

    for (std::size_t field_index = 0; field_index < input.td_ttype.size(); ++field_index)
    {
        const int field_type = input.td_ttype[field_index];
        const std::size_t occurrence = occurrences[field_type]++;
        std::unique_ptr<TDFieldProfile> profile;
        int subdivisions = 1;
        if (field_type == 0)
        {
            const double omega = angular_frequency(input.td_gauss_freq.at(occurrence));
            subdivisions = integration_subdivisions(omega, dt, input.td_stype);
            profile.reset(new TDGaussianProfile(omega,
                                                input.td_gauss_phase.at(occurrence),
                                                input.td_gauss_sigma.at(occurrence) / ModuleBase::AU_to_FS,
                                                input.td_gauss_t0.at(occurrence),
                                                field_amplitude(input.td_gauss_amp.at(occurrence)),
                                                dt));
        }
        else if (field_type == 1)
        {
            const double omega = angular_frequency(input.td_trape_freq.at(occurrence));
            subdivisions = integration_subdivisions(omega, dt, input.td_stype);
            profile.reset(new TDTrapezoidProfile(omega,
                                                 input.td_trape_phase.at(occurrence),
                                                 input.td_trape_t1.at(occurrence),
                                                 input.td_trape_t2.at(occurrence),
                                                 input.td_trape_t3.at(occurrence),
                                                 field_amplitude(input.td_trape_amp.at(occurrence))));
        }
        else if (field_type == 2)
        {
            const double omega1 = angular_frequency(input.td_trigo_freq1.at(occurrence));
            subdivisions = integration_subdivisions(omega1, dt, input.td_stype);
            profile.reset(new TDTrigonometricProfile(omega1,
                                                     angular_frequency(input.td_trigo_freq2.at(occurrence)),
                                                     input.td_trigo_phase1.at(occurrence),
                                                     input.td_trigo_phase2.at(occurrence),
                                                     field_amplitude(input.td_trigo_amp.at(occurrence))));
        }
        else if (field_type == 3)
        {
            subdivisions = input.td_stype == 0 ? 1 : 2;
            profile.reset(new TDHeavisideProfile(input.td_heavi_t0.at(occurrence), field_amplitude(input.td_heavi_amp.at(occurrence))));
        }
        else if (field_type == 4)
        {
            const double omega = angular_frequency(input.td_supersine_freq.at(occurrence));
            subdivisions = integration_subdivisions(omega, dt, input.td_stype);
            profile.reset(new TDSupersineProfile(omega,
                                                 input.td_supersine_phase.at(occurrence),
                                                 input.td_supersine_sigma.at(occurrence),
                                                 input.td_supersine_tstart.at(occurrence),
                                                 input.td_supersine_tend.at(occurrence),
                                                 field_amplitude(input.td_supersine_amp.at(occurrence)),
                                                 dt));
        }

        fields.push_back(TDField(input.td_vext_dire.at(field_index) - 1, std::move(profile), subdivisions));
    }

    return std::shared_ptr<TDFieldManager>(new TDFieldManager(input.td_vext,
                                                              input.td_stype,
                                                              input.td_tstart,
                                                              input.td_tend,
                                                              dt,
                                                              input.td_lcut1,
                                                              input.td_lcut2,
                                                              std::move(fields)));
}

} // namespace elecstate
