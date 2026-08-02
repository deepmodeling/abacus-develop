#include "source_lcao/module_ri/sternheimer_supercell_perturbation.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>

namespace ModuleRI
{
namespace
{

int parse_integer(const std::string& text, const std::string& key)
{
    std::size_t parsed = 0;
    int value = 0;
    try
    {
        value = std::stoi(text, &parsed);
    }
    catch (const std::exception&)
    {
        throw std::invalid_argument("Invalid integer for supercell translation key " + key + ".");
    }
    if (parsed != text.size())
    {
        throw std::invalid_argument("Invalid integer for supercell translation key " + key + ".");
    }
    return value;
}

double parse_double(const std::string& text, const std::string& key)
{
    std::size_t parsed = 0;
    double value = 0.0;
    try
    {
        value = std::stod(text, &parsed);
    }
    catch (const std::exception&)
    {
        throw std::invalid_argument("Invalid floating-point value for supercell translation key " + key + ".");
    }
    if (parsed != text.size() || !std::isfinite(value))
    {
        throw std::invalid_argument("Invalid floating-point value for supercell translation key " + key + ".");
    }
    return value;
}

template <typename T, typename Parser>
std::array<T, 3> parse_triplet(const std::string& text,
                               const char separator,
                               const std::string& key,
                               Parser parser)
{
    std::array<T, 3> result{};
    std::istringstream input(text);
    std::string item;
    for (int i = 0; i != 3; ++i)
    {
        if (!std::getline(input, item, separator) || item.empty())
        {
            throw std::invalid_argument("Invalid triplet for supercell translation key " + key + ".");
        }
        result[static_cast<std::size_t>(i)] = parser(item, key);
    }
    if (std::getline(input, item, separator))
    {
        throw std::invalid_argument("Invalid triplet for supercell translation key " + key + ".");
    }
    return result;
}

void validate_config(const SternheimerSupercellTranslationSum& config)
{
    if (config.atoms_per_primitive <= 0 || config.basis_atom < 0
        || config.basis_atom >= config.atoms_per_primitive || config.channel_within_atom < 0)
    {
        throw std::invalid_argument("Invalid supercell translation atom or channel selector.");
    }
    for (int dimension = 0; dimension != 3; ++dimension)
    {
        const int repeat = config.repeats[static_cast<std::size_t>(dimension)];
        const double q = config.primitive_qpoint[static_cast<std::size_t>(dimension)];
        if (repeat <= 0 || std::abs(q * repeat - std::round(q * repeat)) > 1.0e-10)
        {
            throw std::invalid_argument("Supercell translation q point is not commensurate with the repeats.");
        }
    }
}

const SternheimerABFBlochGridChannel& channel_within_atom(
    const std::vector<SternheimerABFBlochGridChannel>& channels,
    const int atom_index,
    const int ordinal)
{
    int current = 0;
    const SternheimerABFBlochGridChannel* selected = nullptr;
    for (const auto& channel: channels)
    {
        if (channel.atom_index != atom_index)
        {
            continue;
        }
        if (current == ordinal)
        {
            if (selected != nullptr)
            {
                throw std::invalid_argument("Duplicated supercell translation channel ordinal.");
            }
            selected = &channel;
        }
        ++current;
    }
    if (selected == nullptr)
    {
        throw std::invalid_argument("Missing equivalent-atom channel in supercell translation sum.");
    }
    return *selected;
}

} // namespace

int sternheimer_supercell_primitive_cell_count(const SternheimerSupercellTranslationSum& config)
{
    validate_config(config);
    return config.repeats[0] * config.repeats[1] * config.repeats[2];
}

SternheimerSupercellTranslationSum parse_sternheimer_supercell_translation_sum(
    const std::string& specification)
{
    std::map<std::string, std::string> fields;
    std::istringstream input(specification);
    std::string field;
    while (std::getline(input, field, ','))
    {
        const std::size_t separator = field.find('=');
        if (separator == std::string::npos || separator == 0 || separator + 1 == field.size())
        {
            throw std::invalid_argument("Invalid supercell translation perturbation field.");
        }
        if (!fields.insert({field.substr(0, separator), field.substr(separator + 1)}).second)
        {
            throw std::invalid_argument("Duplicated supercell translation perturbation key.");
        }
    }
    const std::set<std::string> required{
        "repeats", "primitive_q", "atoms_per_primitive", "basis_atom", "channel_within_atom"};
    if (fields.size() != required.size())
    {
        throw std::invalid_argument("Incomplete supercell translation perturbation specification.");
    }
    for (const auto& entry: fields)
    {
        if (required.count(entry.first) == 0)
        {
            throw std::invalid_argument("Unknown supercell translation perturbation key " + entry.first + ".");
        }
    }

    SternheimerSupercellTranslationSum config;
    config.repeats = parse_triplet<int>(fields.at("repeats"), 'x', "repeats", parse_integer);
    config.primitive_qpoint
        = parse_triplet<double>(fields.at("primitive_q"), ':', "primitive_q", parse_double);
    config.atoms_per_primitive = parse_integer(fields.at("atoms_per_primitive"), "atoms_per_primitive");
    config.basis_atom = parse_integer(fields.at("basis_atom"), "basis_atom");
    config.channel_within_atom = parse_integer(fields.at("channel_within_atom"), "channel_within_atom");
    validate_config(config);
    return config;
}

SternheimerABFBlochGridChannel combine_sternheimer_supercell_translation_channel(
    const std::vector<SternheimerABFBlochGridChannel>& channels,
    const SternheimerSupercellTranslationSum& config)
{
    const int cell_count = sternheimer_supercell_primitive_cell_count(config);
    if (channels.empty())
    {
        throw std::invalid_argument("Supercell translation sum requires auxiliary channels.");
    }

    SternheimerABFBlochGridChannel result;
    bool initialized = false;
    const double normalization = 1.0 / std::sqrt(static_cast<double>(cell_count));
    const double two_pi = 2.0 * std::acos(-1.0);
    for (int ix = 0; ix != config.repeats[0]; ++ix)
    {
        for (int iy = 0; iy != config.repeats[1]; ++iy)
        {
            for (int iz = 0; iz != config.repeats[2]; ++iz)
            {
                const int cell_index = (ix * config.repeats[1] + iy) * config.repeats[2] + iz;
                const int atom_index = cell_index * config.atoms_per_primitive + config.basis_atom;
                const auto& channel = channel_within_atom(channels, atom_index, config.channel_within_atom);
                if (!initialized)
                {
                    result = channel;
                    result.channel_index = 0;
                    result.atom_index = config.basis_atom;
                    result.atom_local_index = config.channel_within_atom;
                    result.label = "supercell_translation_sum:" + channel.label;
                    result.potential_r.assign(channel.potential_r.size(), std::complex<double>(0.0, 0.0));
                    initialized = true;
                }
                if (channel.potential_r.size() != result.potential_r.size()
                    || channel.type_index != result.type_index
                    || channel.angular_momentum != result.angular_momentum
                    || channel.radial_index != result.radial_index
                    || channel.magnetic_index != result.magnetic_index)
                {
                    throw std::invalid_argument("Incompatible equivalent-atom channels in supercell translation sum.");
                }
                const double angle = two_pi * (config.primitive_qpoint[0] * ix
                                               + config.primitive_qpoint[1] * iy
                                               + config.primitive_qpoint[2] * iz);
                const std::complex<double> phase = normalization * std::exp(std::complex<double>(0.0, angle));
                for (std::size_t ir = 0; ir != result.potential_r.size(); ++ir)
                {
                    result.potential_r[ir] += phase * channel.potential_r[ir];
                }
            }
        }
    }
    result.max_abs = 0.0;
    for (const auto& value: result.potential_r)
    {
        result.max_abs = std::max(result.max_abs, std::abs(value));
    }
    return result;
}

std::vector<SternheimerABFBlochGridChannel>
combine_all_sternheimer_supercell_translation_channels(
    const std::vector<SternheimerABFBlochGridChannel>& channels,
    const SternheimerSupercellTranslationSum& config)
{
    sternheimer_supercell_primitive_cell_count(config);
    if (channels.empty())
    {
        throw std::invalid_argument("Supercell translation sum requires auxiliary channels.");
    }

    std::vector<SternheimerABFBlochGridChannel> combined;
    for (int basis_atom = 0; basis_atom != config.atoms_per_primitive; ++basis_atom)
    {
        const int channels_on_atom = static_cast<int>(std::count_if(
            channels.begin(), channels.end(), [basis_atom](const auto& channel) {
                return channel.atom_index == basis_atom;
            }));
        if (channels_on_atom <= 0)
        {
            throw std::invalid_argument(
                "Supercell translation sum found no auxiliary channels on a primitive atom.");
        }
        for (int ordinal = 0; ordinal != channels_on_atom; ++ordinal)
        {
            SternheimerSupercellTranslationSum channel_config = config;
            channel_config.basis_atom = basis_atom;
            channel_config.channel_within_atom = ordinal;
            SternheimerABFBlochGridChannel channel
                = combine_sternheimer_supercell_translation_channel(channels, channel_config);
            channel.channel_index = static_cast<int>(combined.size());
            channel.atom_index = basis_atom;
            channel.atom_local_index = ordinal;
            combined.push_back(std::move(channel));
        }
    }
    return combined;
}

} // namespace ModuleRI
