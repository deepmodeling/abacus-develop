#include "td_field_manager.h"

#include "source_base/constants.h"
#include "source_base/math_integral.h"
#include "source_base/tool_quit.h"
#include "source_io/module_parameter/input_parameter.h"
#include "td_field_profiles.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <utility>

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
            const double omega = angular_frequency(input.td_supsine_freq.at(occurrence));
            subdivisions = integration_subdivisions(omega, dt, input.td_stype);
            profile.reset(new TDSupersineProfile(omega,
                                                 input.td_supsine_phase.at(occurrence),
                                                 input.td_supsine_sigma.at(occurrence),
                                                 input.td_supsine_tstart.at(occurrence),
                                                 input.td_supsine_tend.at(occurrence),
                                                 field_amplitude(input.td_supsine_amp.at(occurrence)),
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
