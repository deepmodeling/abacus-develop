#include "source_dftb/skf_data.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace ModuleDFTB
{
namespace
{

std::string trim(const std::string& input)
{
    const auto first = input.find_first_not_of(" \t\r\n");
    if (first == std::string::npos)
    {
        return {};
    }
    const auto last = input.find_last_not_of(" \t\r\n");
    return input.substr(first, last - first + 1);
}

std::vector<std::string> read_lines(const std::string& filename)
{
    std::ifstream input(filename);
    if (!input)
    {
        throw std::runtime_error("Cannot open SKF file: " + filename);
    }
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(input, line))
    {
        lines.push_back(line);
    }
    if (!input.eof())
    {
        throw std::runtime_error("I/O error while reading SKF file: " + filename);
    }
    return lines;
}

std::vector<double> parse_numbers(const std::string& line,
                                  const std::string& filename,
                                  std::size_t line_number)
{
    std::string normalized = line;
    for (char& ch : normalized)
    {
        if (ch == ',' || ch == '\t')
        {
            ch = ' ';
        }
        else if (ch == 'd' || ch == 'D')
        {
            ch = 'E';
        }
    }

    std::istringstream stream(normalized);
    std::vector<double> values;
    std::string token;
    while (stream >> token)
    {
        std::size_t parsed = 0;
        double value = 0.0;
        try
        {
            value = std::stod(token, &parsed);
        }
        catch (const std::exception&)
        {
            throw std::runtime_error(filename + ":" + std::to_string(line_number)
                                     + ": expected numeric SKF data, found '" + token + "'");
        }
        if (parsed != token.size() || !std::isfinite(value))
        {
            throw std::runtime_error(filename + ":" + std::to_string(line_number)
                                     + ": invalid or non-finite SKF value '" + token + "'");
        }
        values.push_back(value);
    }
    return values;
}

std::vector<double> require_numbers(const std::vector<std::string>& lines,
                                    std::size_t index,
                                    std::size_t minimum,
                                    const std::string& filename,
                                    const std::string& description)
{
    if (index >= lines.size())
    {
        throw std::runtime_error(filename + ": missing " + description);
    }
    auto values = parse_numbers(lines[index], filename, index + 1);
    if (values.size() < minimum)
    {
        throw std::runtime_error(filename + ":" + std::to_string(index + 1) + ": " + description
                                 + " has " + std::to_string(values.size()) + " numeric fields; expected at least "
                                 + std::to_string(minimum));
    }
    return values;
}

std::size_t find_spline_tag(const std::vector<std::string>& lines,
                            std::size_t begin,
                            const std::string& filename)
{
    for (std::size_t i = begin; i < lines.size(); ++i)
    {
        if (trim(lines[i]) == "Spline")
        {
            return i;
        }
    }
    throw std::runtime_error(filename + ": legacy SKF repulsive Spline section was not found");
}

void parse_repulsive_spline(const std::vector<std::string>& lines,
                            std::size_t spline_tag,
                            const std::string& filename,
                            SkfRepulsiveSpline& spline)
{
    const auto header = require_numbers(lines, spline_tag + 1, 2, filename, "spline count and cutoff");
    const auto n_intervals_as_double = header[0];
    if (n_intervals_as_double < 2.0 || std::floor(n_intervals_as_double) != n_intervals_as_double)
    {
        throw std::runtime_error(filename + ": invalid repulsive spline interval count");
    }
    const auto n_intervals = static_cast<std::size_t>(n_intervals_as_double);
    spline.cutoff_bohr = header[1];

    const auto exponential = require_numbers(lines, spline_tag + 2, 3, filename, "spline exponential coefficients");
    std::copy_n(exponential.begin(), 3, spline.exponential_coefficients.begin());

    spline.interval_starts_bohr.resize(n_intervals);
    spline.interval_ends_bohr.resize(n_intervals);
    spline.cubic_coefficients.resize(n_intervals - 1);
    for (std::size_t interval = 0; interval + 1 < n_intervals; ++interval)
    {
        const auto fields = require_numbers(lines, spline_tag + 3 + interval, 6, filename, "cubic spline interval");
        spline.interval_starts_bohr[interval] = fields[0];
        spline.interval_ends_bohr[interval] = fields[1];
        std::copy_n(fields.begin() + 2, 4, spline.cubic_coefficients[interval].begin());
    }

    const auto tail_index = spline_tag + 3 + n_intervals - 1;
    const auto tail = require_numbers(lines, tail_index, 8, filename, "repulsive polynomial tail");
    spline.interval_starts_bohr[n_intervals - 1] = tail[0];
    spline.interval_ends_bohr[n_intervals - 1] = tail[1];
    std::copy_n(tail.begin() + 2, 6, spline.tail_coefficients.begin());
    spline.cutoff_bohr = spline.interval_ends_bohr.back();

    if (!(spline.cutoff_bohr > 0.0))
    {
        throw std::runtime_error(filename + ": repulsive spline cutoff must be positive");
    }
    for (std::size_t i = 0; i < n_intervals; ++i)
    {
        if (!(spline.interval_ends_bohr[i] > spline.interval_starts_bohr[i]))
        {
            throw std::runtime_error(filename + ": repulsive spline interval has non-positive width");
        }
        if (i + 1 < n_intervals
            && std::abs(spline.interval_ends_bohr[i] - spline.interval_starts_bohr[i + 1]) > 1.0e-8)
        {
            throw std::runtime_error(filename + ": repulsive spline intervals are discontinuous in radius");
        }
    }
}

} // namespace

SkfRepulsiveValue SkfRepulsiveSpline::evaluate(double distance_bohr) const
{
    if (distance_bohr < 0.0 || distance_bohr >= cutoff_bohr || interval_starts_bohr.empty())
    {
        return {};
    }

    if (distance_bohr < interval_starts_bohr.front())
    {
        const double exponential = std::exp(-exponential_coefficients[0] * distance_bohr
                                            + exponential_coefficients[1]);
        SkfRepulsiveValue value;
        value.energy_hartree = exponential + exponential_coefficients[2];
        value.derivative_hartree_per_bohr = -exponential_coefficients[0] * exponential;
        return value;
    }

    if (distance_bohr >= interval_starts_bohr.back())
    {
        const double dr = distance_bohr - interval_starts_bohr.back();
        double energy = tail_coefficients[5];
        for (int order = 4; order >= 0; --order)
        {
            energy = energy * dr + tail_coefficients[static_cast<std::size_t>(order)];
        }
        double derivative = 5.0 * tail_coefficients[5];
        for (int order = 4; order >= 1; --order)
        {
            derivative = derivative * dr + static_cast<double>(order) * tail_coefficients[static_cast<std::size_t>(order)];
        }
        SkfRepulsiveValue value;
        value.energy_hartree = energy;
        value.derivative_hartree_per_bohr = derivative;
        return value;
    }

    const auto upper = std::upper_bound(interval_starts_bohr.begin(), interval_starts_bohr.end(), distance_bohr);
    const auto interval = static_cast<std::size_t>(upper - interval_starts_bohr.begin() - 1);
    const double dr = distance_bohr - interval_starts_bohr[interval];
    const auto& c = cubic_coefficients[interval];
    const double energy = c[0] + dr * (c[1] + dr * (c[2] + dr * c[3]));
    const double derivative = c[1] + dr * (2.0 * c[2] + dr * 3.0 * c[3]);
    SkfRepulsiveValue value;
    value.energy_hartree = energy;
    value.derivative_hartree_per_bohr = derivative;
    return value;
}

double SkfData::valence_electron_count() const
{
    if (!has_atomic_data)
    {
        throw std::runtime_error("Neutral valence occupations are absent from heteronuclear SKF file: " + filename);
    }
    return reference_occupations[0] + reference_occupations[1] + reference_occupations[2];
}

SkfIntegralEvaluation SkfData::evaluate_integrals(double distance_bohr) const
{
    constexpr std::size_t interpolation_points = 8;
    constexpr double distance_fudge_bohr = 1.0;
    constexpr double derivative_step_bohr = 1.0e-5;

    if (!std::isfinite(distance_bohr) || distance_bohr < 0.0)
    {
        throw std::invalid_argument("SKF integral distance must be finite and non-negative");
    }
    if (hamiltonian.size() != overlap.size() || hamiltonian.size() < interpolation_points
        || !(grid_spacing_bohr > 0.0))
    {
        throw std::runtime_error("SKF data are incomplete or too short for DFTB+ 8-point interpolation: " + filename);
    }

    SkfIntegralEvaluation result;
    const std::size_t n_grid = hamiltonian.size();
    const double grid_end_bohr = static_cast<double>(n_grid) * grid_spacing_bohr;
    const double cutoff_bohr = grid_end_bohr + distance_fudge_bohr;
    if (distance_bohr >= cutoff_bohr)
    {
        return result;
    }

    struct PolynomialEvaluation
    {
        std::array<double, 10> value{};
        std::array<double, 10> first{};
        std::array<double, 10> second{};
    };
    const auto interpolate = [this](const std::vector<std::array<double, 10>>& table,
                                    std::size_t first_row,
                                    double distance) {
        PolynomialEvaluation evaluation;
        for (std::size_t column = 0; column < 10; ++column)
        {
            std::array<double, interpolation_points> coefficient{};
            for (std::size_t point = 0; point < interpolation_points; ++point)
            {
                coefficient[point] = table[first_row + point][column];
            }
            for (std::size_t order = 1; order < interpolation_points; ++order)
            {
                for (std::size_t point = interpolation_points - 1; point >= order; --point)
                {
                    coefficient[point] = (coefficient[point] - coefficient[point - 1])
                                         / static_cast<double>(order);
                }
            }

            // Newton form on the dimensionless, unit-spaced grid. The divided
            // differences above are mathematically equivalent to DFTB+'s
            // Neville interpolation but also give analytic first/second derivatives.
            const double t = distance / grid_spacing_bohr - static_cast<double>(first_row + 1);
            double value = coefficient[interpolation_points - 1];
            double first = 0.0;
            double second = 0.0;
            for (std::size_t point = interpolation_points - 1; point-- > 0;)
            {
                const double old_value = value;
                const double old_first = first;
                const double factor = t - static_cast<double>(point);
                second = second * factor + 2.0 * old_first;
                first = old_first * factor + old_value;
                value = old_value * factor + coefficient[point];
            }
            evaluation.value[column] = value;
            evaluation.first[column] = first / grid_spacing_bohr;
            evaluation.second[column] = second / (grid_spacing_bohr * grid_spacing_bohr);
        }
        return evaluation;
    };

    if (distance_bohr < grid_end_bohr)
    {
        const auto grid_index = static_cast<std::size_t>(std::floor(distance_bohr / grid_spacing_bohr));
        const auto last_exclusive = std::max(interpolation_points,
                                             std::min(n_grid, grid_index + interpolation_points / 2));
        const std::size_t first_row = last_exclusive - interpolation_points;
        const auto h = interpolate(hamiltonian, first_row, distance_bohr);
        const auto s = interpolate(overlap, first_row, distance_bohr);
        result.hamiltonian = h.value;
        result.overlap = s.value;
        result.d_hamiltonian_dr = h.first;
        result.d_overlap_dr = s.first;
        result.d2_hamiltonian_dr2 = h.second;
        result.d2_overlap_dr2 = s.second;
        return result;
    }

    const std::size_t first_row = n_grid - interpolation_points;
    const auto h_at_end = interpolate(hamiltonian, first_row, grid_end_bohr);
    const auto h_before = interpolate(hamiltonian, first_row, grid_end_bohr - derivative_step_bohr);
    const auto h_after = interpolate(hamiltonian, first_row, grid_end_bohr + derivative_step_bohr);
    const auto s_at_end = interpolate(overlap, first_row, grid_end_bohr);
    const auto s_before = interpolate(overlap, first_row, grid_end_bohr - derivative_step_bohr);
    const auto s_after = interpolate(overlap, first_row, grid_end_bohr + derivative_step_bohr);
    const double xr = cutoff_bohr - distance_bohr;

    const auto smooth_to_zero = [xr](double y, double yp, double ypp, double& value, double& derivative, double& second) {
        const double dx1 = -yp;
        const double dx2 = ypp;
        const double dd = 10.0 * y - 4.0 * dx1 + 0.5 * dx2;
        const double ee = -15.0 * y + 7.0 * dx1 - dx2;
        const double ff = 6.0 * y - 3.0 * dx1 + 0.5 * dx2;
        value = ((ff * xr + ee) * xr + dd) * xr * xr * xr;
        derivative = -(3.0 * dd * xr * xr + 4.0 * ee * xr * xr * xr + 5.0 * ff * xr * xr * xr * xr);
        second = 6.0 * dd * xr + 12.0 * ee * xr * xr + 20.0 * ff * xr * xr * xr;
    };
    for (std::size_t column = 0; column < 10; ++column)
    {
        const double h_first = (h_after.value[column] - h_before.value[column]) / (2.0 * derivative_step_bohr);
        const double h_second = (h_after.value[column] + h_before.value[column] - 2.0 * h_at_end.value[column])
                                / (derivative_step_bohr * derivative_step_bohr);
        smooth_to_zero(h_at_end.value[column], h_first, h_second,
                       result.hamiltonian[column], result.d_hamiltonian_dr[column], result.d2_hamiltonian_dr2[column]);
        double s_value = 0.0;
        double s_derivative = 0.0;
        double s_second = 0.0;
        const double s_first = (s_after.value[column] - s_before.value[column]) / (2.0 * derivative_step_bohr);
        const double s_second_fit = (s_after.value[column] + s_before.value[column] - 2.0 * s_at_end.value[column])
                                / (derivative_step_bohr * derivative_step_bohr);
        smooth_to_zero(s_at_end.value[column], s_first, s_second_fit,
                       s_value, s_derivative, s_second);
        result.overlap[column] = s_value;
        result.d_overlap_dr[column] = s_derivative;
        result.d2_overlap_dr2[column] = s_second;
    }
    return result;
}

SkfData SkfData::read_legacy(const std::string& filename, bool is_homonuclear)
{
    const auto lines = read_lines(filename);
    if (lines.empty())
    {
        throw std::runtime_error("Empty SKF file: " + filename);
    }
    if (!trim(lines.front()).empty() && trim(lines.front()).front() == '@')
    {
        throw std::runtime_error(filename + ": '@' extended s/p/d/f SKF format is not supported by the legacy reader");
    }

    SkfData data;
    data.filename = filename;
    data.homonuclear = is_homonuclear;

    const auto grid = require_numbers(lines, 0, 2, filename, "grid spacing and grid point count");
    data.grid_spacing_bohr = grid[0];
    if (!(data.grid_spacing_bohr > 0.0))
    {
        throw std::runtime_error(filename + ": SKF grid spacing must be positive");
    }
    const double n_grid_as_double = grid[1];
    if (n_grid_as_double < 2.0 || std::floor(n_grid_as_double) != n_grid_as_double)
    {
        throw std::runtime_error(filename + ": invalid SKF grid point count");
    }
    data.declared_grid_points = static_cast<std::size_t>(n_grid_as_double);
    // Match DFTB+'s legacy reader, which decrements the file header count once.
    const std::size_t n_table_rows = data.declared_grid_points - 1;

    std::size_t next = 1;
    if (is_homonuclear)
    {
        const auto atom = require_numbers(lines, next++, 10, filename, "homonuclear atomic data");
        // Legacy SKF stores each shell group in reverse order: d, p, s.
        data.onsite_hartree = {atom[2], atom[1], atom[0]};
        data.hubbard_u_hartree = {atom[6], atom[5], atom[4]};
        data.reference_occupations = {atom[9], atom[8], atom[7]};
        const auto metadata = require_numbers(lines, next++, 20, filename, "homonuclear mass/repulsive metadata");
        data.mass_amu = metadata[0];
        data.has_atomic_data = true;
    }
    else
    {
        // The heteronuclear metadata line holds legacy polynomial-repulsive
        // fields; this stage only uses the Spline repulsive section below.
        (void)require_numbers(lines, next++, 20, filename, "heteronuclear repulsive metadata");
    }

    data.hamiltonian.reserve(n_table_rows);
    data.overlap.reserve(n_table_rows);
    for (std::size_t row = 0; row < n_table_rows; ++row, ++next)
    {
        const auto values = require_numbers(lines, next, 20, filename, "20-column legacy H/S integral row");
        std::array<double, 10> h{};
        std::array<double, 10> s{};
        std::copy_n(values.begin(), 10, h.begin());
        std::copy_n(values.begin() + 10, 10, s.begin());
        data.hamiltonian.push_back(h);
        data.overlap.push_back(s);
    }

    const auto spline_tag = find_spline_tag(lines, next, filename);
    parse_repulsive_spline(lines, spline_tag, filename, data.repulsive);
    data.has_repulsive_spline = true;
    return data;
}

void validate_pair_directions(const SkfData& ab, const SkfData& ba, double tolerance)
{
    if (ab.homonuclear || ba.homonuclear)
    {
        throw std::invalid_argument("Directional SKF validation requires two heteronuclear files");
    }
    if (ab.declared_grid_points != ba.declared_grid_points
        || std::abs(ab.grid_spacing_bohr - ba.grid_spacing_bohr) > tolerance)
    {
        throw std::runtime_error("Incompatible SKF grid spacing or length for " + ab.filename + " and " + ba.filename);
    }
    if (ab.hamiltonian.size() != ba.hamiltonian.size())
    {
        throw std::runtime_error("Incompatible SKF table row count for " + ab.filename + " and " + ba.filename);
    }
    if (ab.has_repulsive_spline != ba.has_repulsive_spline)
    {
        throw std::runtime_error("Only one directional SKF file contains a repulsive spline");
    }
    if (!ab.has_repulsive_spline)
    {
        return;
    }

    const auto& a = ab.repulsive;
    const auto& b = ba.repulsive;
    const auto mismatch = [tolerance](double x, double y) { return std::abs(x - y) > tolerance; };
    if (a.interval_starts_bohr.size() != b.interval_starts_bohr.size()
        || a.cubic_coefficients.size() != b.cubic_coefficients.size()
        || mismatch(a.cutoff_bohr, b.cutoff_bohr))
    {
        throw std::runtime_error("Incompatible SKF repulsive spline dimensions or cutoff for "
                                 + ab.filename + " and " + ba.filename);
    }
    for (std::size_t i = 0; i < a.interval_starts_bohr.size(); ++i)
    {
        if (mismatch(a.interval_starts_bohr[i], b.interval_starts_bohr[i])
            || mismatch(a.interval_ends_bohr[i], b.interval_ends_bohr[i]))
        {
            throw std::runtime_error("Incompatible SKF repulsive spline knots for " + ab.filename + " and " + ba.filename);
        }
    }
    for (std::size_t i = 0; i < a.cubic_coefficients.size(); ++i)
    {
        for (std::size_t j = 0; j < a.cubic_coefficients[i].size(); ++j)
        {
            if (mismatch(a.cubic_coefficients[i][j], b.cubic_coefficients[i][j]))
            {
                throw std::runtime_error("Incompatible SKF repulsive spline coefficients for "
                                         + ab.filename + " and " + ba.filename);
            }
        }
    }
    for (std::size_t i = 0; i < a.exponential_coefficients.size(); ++i)
    {
        if (mismatch(a.exponential_coefficients[i], b.exponential_coefficients[i]))
        {
            throw std::runtime_error("Incompatible SKF exponential repulsive coefficients for "
                                     + ab.filename + " and " + ba.filename);
        }
    }
    for (std::size_t i = 0; i < a.tail_coefficients.size(); ++i)
    {
        if (mismatch(a.tail_coefficients[i], b.tail_coefficients[i]))
        {
            throw std::runtime_error("Incompatible SKF repulsive tail coefficients for "
                                     + ab.filename + " and " + ba.filename);
        }
    }
}

} // namespace ModuleDFTB
