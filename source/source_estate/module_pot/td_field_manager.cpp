#include "td_field_manager.h"

#include "source_base/constants.h"
#include "source_base/math_integral.h"
#include "source_base/tool_quit.h"
#include "source_io/module_parameter/input_parameter.h"
#include "td_field_profiles.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <utility>

namespace
{

std::vector<std::string> restart_rows(std::istream& input)
{
    std::vector<std::string> rows;
    std::string line;
    while (std::getline(input, line))
    {
        const std::size_t comment = line.find('#');
        if (comment != std::string::npos)
        {
            line.erase(comment);
        }
        if (line.find_first_not_of(" \t\r\n") == std::string::npos)
        {
            continue;
        }
        rows.push_back(line);
        if (rows.size() > 3)
        {
            break;
        }
    }
    if (rows.size() != 3)
    {
        ModuleBase::WARNING_QUIT("TDFieldManager::read_restart", "Expected three data rows with 5, 3 and 3 fields.");
    }
    return rows;
}

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
    return amplitude * ModuleBase::BOHR_TO_A / ModuleBase::Hartree_to_eV;
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
    : enabled_(enabled), gauge_(gauge), start_step_(start_step), end_step_(end_step), dt_ha_(dt), length_cut1_(length_cut1),
      length_cut2_(length_cut2), fields_(std::move(fields)), current_step_(-1), active_(false)
{
    field_vals_ha_.resize(fields_.size(), 0.0);
}

void TDFieldManager::sample_field(const int step)
{
    current_step_ = step;
    active_ = enabled_ && step >= start_step_ && step <= end_step_;
    std::fill(field_vals_ha_.begin(), field_vals_ha_.end(), 0.0);
    efield_ha_.set(0.0, 0.0, 0.0);
    if (!active_)
    {
        return;
    }
    for (std::size_t index = 0; index < fields_.size(); ++index)
    {
        const TDField& field = fields_[index];
        const double time_ha = step * dt_ha_;
        const TDFieldSample sample(step, 0, field.subdivisions(), time_ha);
        field_vals_ha_[index] = field.electric_field(sample);
        efield_ha_[field.direction()] += field_vals_ha_[index];
    }
}

void TDFieldManager::prepare_sample(const int step)
{
    if (step == current_step_)
    {
        return;
    }
    if (step != current_step_ + 1)
    {
        ModuleBase::WARNING_QUIT("TDFieldManager::prepare_sample", "Invalid electronic-step sequence.");
    }
    sample_field(step);
}

void TDFieldManager::prepare_interval(const int left_step, const int sample_step)
{
    if (left_step == interval_left_ && sample_step == current_step_)
    {
        return;
    }
    if (left_step != interval_left_ + 1 || (sample_step != left_step && sample_step != left_step + 1))
    {
        ModuleBase::WARNING_QUIT("TDFieldManager::prepare_interval", "Invalid field interval sequence.");
    }
    A_left_ha_ = A_right_ha_;
    if (A_samples_ha_.empty() && enabled_ && left_step >= start_step_ && left_step <= end_step_)
    {
        for (const TDField& field: fields_)
        {
            const int subdivisions = field.subdivisions();
            const double integration_dt = dt_ha_ / subdivisions;
            std::vector<double> samples(subdivisions + 1);
            for (int node = 0; node <= subdivisions; ++node)
            {
                const double time_ha = (left_step + static_cast<double>(node) / subdivisions) * dt_ha_;
                const TDFieldSample sample(left_step, node, subdivisions, time_ha);
                samples[node] = field.electric_field(sample);
            }
            double integral = 0.0;
            const int count = subdivisions + 1;
            ModuleBase::Integral::Simpson_Integral(count, samples.data(), integration_dt, integral);
            A_right_ha_[field.direction()] -= integral;
        }
    }
    interval_left_ = left_step;
    sample_field(sample_step);
    select_A_prop();
}

void TDFieldManager::read_restart(const std::string& file_dir)
{
    std::ifstream file((file_dir + "Restart_td.txt").c_str());
    if (!file)
    {
        ModuleBase::WARNING_QUIT("TDFieldManager::read_restart", "Cannot open Restart_td.txt.");
    }
    const std::vector<std::string> rows = restart_rows(file);
    std::istringstream state(rows[0]);
    std::istringstream left(rows[1]);
    std::istringstream right(rows[2]);
    int gauge = -1;
    int file_source = -1;
    double dt = 0.0;
    if (!(state >> current_step_ >> interval_left_ >> gauge >> dt >> file_source) || !(left >> A_left_ha_.x >> A_left_ha_.y >> A_left_ha_.z)
        || !(right >> A_right_ha_.x >> A_right_ha_.y >> A_right_ha_.z) || current_step_ < 0 || interval_left_ < -1
        || (gauge_ == 0 && interval_left_ != -1) || (gauge_ != 0 && interval_left_ != current_step_ && interval_left_ != current_step_ - 1)
        || file_source != !A_samples_ha_.empty() || gauge != gauge_ || !std::isfinite(dt) || std::abs(dt - dt_ha_) > 1.0e-12 * dt_ha_)
    {
        ModuleBase::WARNING_QUIT("TDFieldManager::read_restart", "Invalid or incompatible field restart state.");
    }
    state >> std::ws;
    left >> std::ws;
    right >> std::ws;
    if (!state.eof() || !left.eof() || !right.eof())
    {
        ModuleBase::WARNING_QUIT("TDFieldManager::read_restart", "Unexpected extra data in Restart_td.txt.");
    }
    for (int d = 0; d < 3; ++d)
    {
        if (!std::isfinite(A_left_ha_[d]) || !std::isfinite(A_right_ha_[d]))
        {
            ModuleBase::WARNING_QUIT("TDFieldManager::read_restart", "Non-finite restart vector potential.");
        }
    }
    sample_field(current_step_);
    select_A_prop();
}

void TDFieldManager::write_restart(const std::string& file_dir) const
{
    std::ofstream file((file_dir + "Restart_td.txt").c_str());
    file << std::setprecision(17) << "# Hartree atomic units; step indices start at 0.\n"
         << "# gauge: 0=length, 1=velocity, 2=hybrid; source: 0=field, 1=file\n"
         << "# step  left_step(-1=none)  gauge  dt  source\n"
         << current_step_ << " " << interval_left_ << " " << gauge_ << " " << dt_ha_ << " " << !A_samples_ha_.empty() << "\n"
         << A_left_ha_.x << " " << A_left_ha_.y << " " << A_left_ha_.z << "  # A_left: x y z\n"
         << A_right_ha_.x << " " << A_right_ha_.y << " " << A_right_ha_.z << "  # A_right: x y z\n";
    if (!file)
    {
        ModuleBase::WARNING_QUIT("TDFieldManager::write_restart", "Cannot write field restart.");
    }
}

void TDFieldManager::set_A_samples(const std::vector<ModuleBase::Vector3<double>>& samples_ha)
{
    if (current_step_ != -1 || samples_ha.empty())
    {
        ModuleBase::WARNING_QUIT("TDFieldManager::set_A_samples", "Supply nonempty propagation samples before initialization.");
    }
    for (const ModuleBase::Vector3<double>& A: samples_ha)
    {
        for (int d = 0; d < 3; ++d)
        {
            if (!std::isfinite(A[d]))
            {
                ModuleBase::WARNING_QUIT("TDFieldManager::set_A_samples", "Non-finite propagation sample.");
            }
        }
    }
    A_samples_ha_ = samples_ha;
}

void TDFieldManager::select_A_prop()
{
    if (A_samples_ha_.empty())
    {
        A_prop_ha_ = (A_left_ha_ + A_right_ha_) * 0.5;
    }
    else
    {
        const std::size_t index = std::min(static_cast<std::size_t>(current_step_), A_samples_ha_.size() - 1);
        A_prop_ha_ = A_samples_ha_[index];
    }
}

int TDFieldManager::gauge() const
{
    return gauge_;
}
int TDFieldManager::current_step() const
{
    return current_step_;
}
double TDFieldManager::dt_ha() const
{
    return dt_ha_;
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
const std::vector<double>& TDFieldManager::field_vals_ha() const
{
    return field_vals_ha_;
}
const ModuleBase::Vector3<double>& TDFieldManager::A_left_ha() const
{
    if (!A_samples_ha_.empty())
    {
        ModuleBase::WARNING_QUIT("TDFieldManager::A_left_ha", "File propagation samples do not define endpoints.");
    }
    return A_left_ha_;
}
const ModuleBase::Vector3<double>& TDFieldManager::A_right_ha() const
{
    if (!A_samples_ha_.empty())
    {
        ModuleBase::WARNING_QUIT("TDFieldManager::A_right_ha", "File propagation samples do not define endpoints.");
    }
    return A_right_ha_;
}
const ModuleBase::Vector3<double>& TDFieldManager::A_prop_ha() const
{
    return A_prop_ha_;
}
const ModuleBase::Vector3<double>& TDFieldManager::efield_ha() const
{
    return efield_ha_;
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
