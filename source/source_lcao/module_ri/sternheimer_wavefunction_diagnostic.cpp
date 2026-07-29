#include "source_lcao/module_ri/sternheimer_wavefunction_diagnostic.h"

#include <cstdint>
#include <fstream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>

namespace ModuleRI
{
namespace
{

constexpr std::int32_t kMagic = -53291471;
constexpr std::int32_t kVersion = 1;

template <typename T>
void write_scalar(std::ofstream& out, const T& value, const std::string& filename)
{
    out.write(reinterpret_cast<const char*>(&value), sizeof(T));
    if (!out.good())
    {
        throw std::runtime_error("Failed to write Sternheimer wavefunction diagnostic " + filename + ".");
    }
}

template <typename T>
T read_scalar(std::ifstream& in, const std::string& filename)
{
    T value{};
    in.read(reinterpret_cast<char*>(&value), sizeof(T));
    if (!in.good())
    {
        throw std::runtime_error("Truncated Sternheimer wavefunction diagnostic " + filename + ".");
    }
    return value;
}

std::size_t checked_grid_size(const SternheimerWavefunctionDiagnostic::Metadata& metadata)
{
    if (metadata.nx <= 0 || metadata.ny <= 0 || metadata.nz <= 0)
    {
        throw std::invalid_argument("Sternheimer wavefunction diagnostic requires positive grid dimensions.");
    }
    const std::size_t nx = static_cast<std::size_t>(metadata.nx);
    const std::size_t ny = static_cast<std::size_t>(metadata.ny);
    const std::size_t nz = static_cast<std::size_t>(metadata.nz);
    if (nx > std::numeric_limits<std::size_t>::max() / ny
        || nx * ny > std::numeric_limits<std::size_t>::max() / nz)
    {
        throw std::overflow_error("Sternheimer wavefunction diagnostic grid size overflows size_t.");
    }
    return nx * ny * nz;
}

void validate_record(const SternheimerWavefunctionDiagnostic::Record& record)
{
    const std::size_t grid_size = checked_grid_size(record.metadata);
    std::set<std::string> names;
    for (const auto& named_vector: record.vectors)
    {
        if (named_vector.first.empty())
        {
            throw std::invalid_argument("Sternheimer wavefunction diagnostic vector name is empty.");
        }
        if (!names.insert(named_vector.first).second)
        {
            throw std::invalid_argument("Sternheimer wavefunction diagnostic vector name is duplicated.");
        }
        if (named_vector.second.size() != grid_size)
        {
            throw std::invalid_argument("Sternheimer wavefunction diagnostic vector size differs from the grid.");
        }
    }
}

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
        throw std::invalid_argument("Invalid integer for Sternheimer wavefunction diagnostic key " + key + ".");
    }
    if (parsed != text.size())
    {
        throw std::invalid_argument("Invalid integer for Sternheimer wavefunction diagnostic key " + key + ".");
    }
    return value;
}

} // namespace

bool SternheimerWavefunctionDiagnostic::Selector::matches(const int candidate_iq,
                                                           const int candidate_ik_full,
                                                           const int candidate_ib,
                                                           const int candidate_ifrequency,
                                                           const int candidate_channel) const
{
    return iq == candidate_iq && ik_full == candidate_ik_full && ib == candidate_ib
           && ifrequency == candidate_ifrequency && channel == candidate_channel;
}

SternheimerWavefunctionDiagnostic::Configuration
SternheimerWavefunctionDiagnostic::parse_configuration(const std::string& specification)
{
    std::map<std::string, std::string> fields;
    std::istringstream input(specification);
    std::string field;
    while (std::getline(input, field, ','))
    {
        const std::size_t separator = field.find('=');
        if (separator == std::string::npos || separator == 0 || separator + 1 == field.size())
        {
            throw std::invalid_argument("Invalid Sternheimer wavefunction diagnostic field.");
        }
        const std::string key = field.substr(0, separator);
        const std::string value = field.substr(separator + 1);
        if (!fields.insert({key, value}).second)
        {
            throw std::invalid_argument("Duplicated Sternheimer wavefunction diagnostic key " + key + ".");
        }
    }
    const std::set<std::string> required{"iq", "ik", "ib", "ifreq", "channel", "out"};
    if (fields.size() != required.size())
    {
        throw std::invalid_argument("Sternheimer wavefunction diagnostic requires iq,ik,ib,ifreq,channel,out.");
    }
    for (const auto& entry: fields)
    {
        if (required.count(entry.first) == 0)
        {
            throw std::invalid_argument("Unknown Sternheimer wavefunction diagnostic key " + entry.first + ".");
        }
    }

    Configuration config;
    config.selector.iq = parse_integer(fields.at("iq"), "iq");
    config.selector.ik_full = parse_integer(fields.at("ik"), "ik");
    config.selector.ib = parse_integer(fields.at("ib"), "ib");
    config.selector.ifrequency = parse_integer(fields.at("ifreq"), "ifreq");
    config.selector.channel = parse_integer(fields.at("channel"), "channel");
    config.output_filename = fields.at("out");
    if (config.selector.iq <= 0 || config.selector.ifrequency <= 0 || config.selector.ik_full < 0
        || config.selector.ib < 0 || config.selector.channel < 0 || config.output_filename.empty())
    {
        throw std::invalid_argument("Invalid Sternheimer wavefunction diagnostic selector range.");
    }
    return config;
}

void SternheimerWavefunctionDiagnostic::write(const std::string& filename, const Record& record)
{
    validate_record(record);
    std::ofstream out(filename.c_str(), std::ios::binary | std::ios::trunc);
    if (!out.good())
    {
        throw std::runtime_error("Failed to open Sternheimer wavefunction diagnostic " + filename + ".");
    }

    write_scalar(out, kMagic, filename);
    write_scalar(out, kVersion, filename);
    const std::int32_t indices[] = {record.metadata.nx,
                                    record.metadata.ny,
                                    record.metadata.nz,
                                    record.metadata.iq,
                                    record.metadata.ik_full,
                                    record.metadata.ib,
                                    record.metadata.ifrequency,
                                    record.metadata.channel};
    for (const std::int32_t value: indices)
    {
        write_scalar(out, value, filename);
    }
    for (const double value: record.metadata.lattice)
    {
        write_scalar(out, value, filename);
    }
    for (const double value: record.metadata.qpoint)
    {
        write_scalar(out, value, filename);
    }
    for (const double value: record.metadata.source_kpoint)
    {
        write_scalar(out, value, filename);
    }
    for (const double value: record.metadata.target_kpoint)
    {
        write_scalar(out, value, filename);
    }
    const double scalars[] = {record.metadata.omega_ha,
                              record.metadata.omega_ry,
                              record.metadata.volume_element,
                              record.metadata.reference_eigenvalue_ry,
                              record.metadata.weighted_occupation,
                              record.metadata.rhs_norm,
                              record.metadata.solver_relative_residual,
                              record.metadata.equation_relative_residual,
                              record.metadata.diagonal_branch_element.real(),
                              record.metadata.diagonal_branch_element.imag()};
    for (const double value: scalars)
    {
        write_scalar(out, value, filename);
    }
    if (record.vectors.size() > static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max()))
    {
        throw std::overflow_error("Too many Sternheimer wavefunction diagnostic vectors.");
    }
    const std::int32_t vector_count = static_cast<std::int32_t>(record.vectors.size());
    write_scalar(out, vector_count, filename);
    for (const auto& named_vector: record.vectors)
    {
        if (named_vector.first.size() > static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max())
            || named_vector.second.size() > static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max()))
        {
            throw std::overflow_error("Sternheimer wavefunction diagnostic vector metadata overflows.");
        }
        const std::int32_t name_size = static_cast<std::int32_t>(named_vector.first.size());
        const std::int64_t vector_size = static_cast<std::int64_t>(named_vector.second.size());
        write_scalar(out, name_size, filename);
        out.write(named_vector.first.data(), name_size);
        write_scalar(out, vector_size, filename);
        out.write(reinterpret_cast<const char*>(named_vector.second.data()),
                  static_cast<std::streamsize>(named_vector.second.size() * sizeof(Complex)));
        if (!out.good())
        {
            throw std::runtime_error("Failed to write Sternheimer wavefunction diagnostic " + filename + ".");
        }
    }
}

SternheimerWavefunctionDiagnostic::Record SternheimerWavefunctionDiagnostic::read(const std::string& filename)
{
    std::ifstream in(filename.c_str(), std::ios::binary);
    if (!in.good())
    {
        throw std::runtime_error("Failed to open Sternheimer wavefunction diagnostic " + filename + ".");
    }
    if (read_scalar<std::int32_t>(in, filename) != kMagic
        || read_scalar<std::int32_t>(in, filename) != kVersion)
    {
        throw std::runtime_error("Unsupported Sternheimer wavefunction diagnostic format in " + filename + ".");
    }

    Record record;
    int* indices[] = {&record.metadata.nx,
                      &record.metadata.ny,
                      &record.metadata.nz,
                      &record.metadata.iq,
                      &record.metadata.ik_full,
                      &record.metadata.ib,
                      &record.metadata.ifrequency,
                      &record.metadata.channel};
    for (int* value: indices)
    {
        *value = read_scalar<std::int32_t>(in, filename);
    }
    for (double& value: record.metadata.lattice)
    {
        value = read_scalar<double>(in, filename);
    }
    for (double& value: record.metadata.qpoint)
    {
        value = read_scalar<double>(in, filename);
    }
    for (double& value: record.metadata.source_kpoint)
    {
        value = read_scalar<double>(in, filename);
    }
    for (double& value: record.metadata.target_kpoint)
    {
        value = read_scalar<double>(in, filename);
    }
    double diagonal_branch_real = 0.0;
    double diagonal_branch_imag = 0.0;
    double* scalars[] = {&record.metadata.omega_ha,
                         &record.metadata.omega_ry,
                         &record.metadata.volume_element,
                         &record.metadata.reference_eigenvalue_ry,
                         &record.metadata.weighted_occupation,
                         &record.metadata.rhs_norm,
                         &record.metadata.solver_relative_residual,
                         &record.metadata.equation_relative_residual,
                         &diagonal_branch_real,
                         &diagonal_branch_imag};
    for (double* value: scalars)
    {
        *value = read_scalar<double>(in, filename);
    }
    record.metadata.diagonal_branch_element = Complex(diagonal_branch_real, diagonal_branch_imag);
    const std::size_t grid_size = checked_grid_size(record.metadata);

    const std::int32_t vector_count = read_scalar<std::int32_t>(in, filename);
    if (vector_count < 0)
    {
        throw std::runtime_error("Negative vector count in Sternheimer wavefunction diagnostic " + filename + ".");
    }
    record.vectors.reserve(static_cast<std::size_t>(vector_count));
    for (std::int32_t ivector = 0; ivector != vector_count; ++ivector)
    {
        const std::int32_t name_size = read_scalar<std::int32_t>(in, filename);
        if (name_size <= 0 || name_size > 4096)
        {
            throw std::runtime_error("Invalid vector name in Sternheimer wavefunction diagnostic " + filename + ".");
        }
        std::string name(static_cast<std::size_t>(name_size), '\0');
        in.read(&name[0], name_size);
        if (!in.good())
        {
            throw std::runtime_error("Truncated Sternheimer wavefunction diagnostic " + filename + ".");
        }
        const std::int64_t vector_size = read_scalar<std::int64_t>(in, filename);
        if (vector_size < 0
            || static_cast<std::uint64_t>(vector_size) != static_cast<std::uint64_t>(grid_size))
        {
            throw std::runtime_error("Vector size differs from the grid in Sternheimer wavefunction diagnostic "
                                     + filename + ".");
        }
        Vector values(static_cast<std::size_t>(vector_size));
        in.read(reinterpret_cast<char*>(values.data()),
                static_cast<std::streamsize>(values.size() * sizeof(Complex)));
        if (!in.good())
        {
            throw std::runtime_error("Truncated Sternheimer wavefunction diagnostic " + filename + ".");
        }
        record.vectors.push_back({std::move(name), std::move(values)});
    }
    validate_record(record);
    return record;
}

} // namespace ModuleRI
