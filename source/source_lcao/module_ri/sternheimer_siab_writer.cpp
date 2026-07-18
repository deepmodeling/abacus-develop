#include "sternheimer_siab_writer.h"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <tuple>

namespace module_ri
{
namespace sternheimer_siab
{
namespace
{

const double hermitian_tolerance = 1.0e-10;

bool finite_complex(const std::complex<double>& value)
{
    return std::isfinite(value.real()) && std::isfinite(value.imag());
}

bool valid_hex(const std::string& value, const std::size_t length)
{
    if (value.size() != length)
    {
        return false;
    }
    return std::all_of(value.begin(), value.end(), [](const char character) {
        return (character >= '0' && character <= '9') || (character >= 'a' && character <= 'f')
               || (character >= 'A' && character <= 'F');
    });
}

bool valid_commit(const std::string& value)
{
    return valid_hex(value, 40) || valid_hex(value, 64);
}

bool ascii_letter(const unsigned char character)
{
    return (character >= 'A' && character <= 'Z') || (character >= 'a' && character <= 'z');
}

bool valid_element_token(const std::string& element)
{
    if (element.empty() || !ascii_letter(static_cast<unsigned char>(element[0])))
    {
        return false;
    }
    for (std::size_t index = 1; index < element.size(); ++index)
    {
        const unsigned char character = static_cast<unsigned char>(element[index]);
        if (!ascii_letter(character) && !(character >= '0' && character <= '9') && character != '_' && character != '+'
            && character != '-')
        {
            return false;
        }
    }
    return true;
}

bool continuation_byte(const unsigned char byte)
{
    return byte >= 0x80 && byte <= 0xBF;
}

bool valid_utf8(const std::string& value)
{
    std::size_t position = 0;
    while (position < value.size())
    {
        const unsigned char first = static_cast<unsigned char>(value[position]);
        if (first <= 0x7F)
        {
            ++position;
            continue;
        }
        if (first >= 0xC2 && first <= 0xDF)
        {
            if (position + 1 >= value.size() || !continuation_byte(static_cast<unsigned char>(value[position + 1])))
            {
                return false;
            }
            position += 2;
            continue;
        }
        if (first >= 0xE0 && first <= 0xEF)
        {
            if (position + 2 >= value.size())
            {
                return false;
            }
            const unsigned char second = static_cast<unsigned char>(value[position + 1]);
            const unsigned char third = static_cast<unsigned char>(value[position + 2]);
            if (!continuation_byte(third) || (first == 0xE0 && (second < 0xA0 || second > 0xBF))
                || (first == 0xED && (second < 0x80 || second > 0x9F))
                || ((first != 0xE0 && first != 0xED) && !continuation_byte(second)))
            {
                return false;
            }
            position += 3;
            continue;
        }
        if (first >= 0xF0 && first <= 0xF4)
        {
            if (position + 3 >= value.size())
            {
                return false;
            }
            const unsigned char second = static_cast<unsigned char>(value[position + 1]);
            const unsigned char third = static_cast<unsigned char>(value[position + 2]);
            const unsigned char fourth = static_cast<unsigned char>(value[position + 3]);
            if (!continuation_byte(third) || !continuation_byte(fourth)
                || (first == 0xF0 && (second < 0x90 || second > 0xBF))
                || (first == 0xF4 && (second < 0x80 || second > 0x8F))
                || ((first != 0xF0 && first != 0xF4) && !continuation_byte(second)))
            {
                return false;
            }
            position += 4;
            continue;
        }
        return false;
    }
    return true;
}

void validate_utf8_field(const std::string& value, const std::string& field)
{
    if (!valid_utf8(value))
    {
        throw std::invalid_argument("Sternheimer SIAB provenance " + field + " must be valid UTF-8");
    }
}

std::size_t validate_blocks(const std::vector<PrimitiveBlock>& blocks)
{
    std::size_t expected_offset = 0;
    for (const PrimitiveBlock& block: blocks)
    {
        if (!valid_element_token(block.element))
        {
            throw std::invalid_argument("Sternheimer SIAB primitive block element must match [A-Za-z][A-Za-z0-9_+-]*");
        }
        if (block.atom_index < 0 || block.l < 0 || block.m < -block.l || block.m > block.l)
        {
            throw std::invalid_argument("Sternheimer SIAB primitive block has invalid atom or angular indices");
        }
        if (block.n_primitive <= 0 || block.offset < 0 || static_cast<std::size_t>(block.offset) != expected_offset)
        {
            throw std::invalid_argument("Sternheimer SIAB primitive block offsets must be contiguous from zero");
        }
        expected_offset += static_cast<std::size_t>(block.n_primitive);
    }
    return expected_offset;
}

void validate_provenance(const Provenance& provenance)
{
    validate_utf8_field(provenance.abacus_commit, "abacus_commit");
    validate_utf8_field(provenance.auxiliary_basis_sha256, "auxiliary_basis_sha256");
    validate_utf8_field(provenance.kernel, "kernel");
    validate_utf8_field(provenance.orbital_sha256, "orbital_sha256");
    validate_utf8_field(provenance.pseudopotential_sha256, "pseudopotential_sha256");
    validate_utf8_field(provenance.spin_convention, "spin_convention");
    validate_utf8_field(provenance.executable_sha256, "executable_sha256");
    if (!valid_commit(provenance.abacus_commit))
    {
        throw std::invalid_argument("Sternheimer SIAB provenance requires a 40- or 64-digit hexadecimal ABACUS commit");
    }
    if (!valid_hex(provenance.auxiliary_basis_sha256, 64) || !valid_hex(provenance.orbital_sha256, 64)
        || !valid_hex(provenance.pseudopotential_sha256, 64))
    {
        throw std::invalid_argument("Sternheimer SIAB provenance hashes must be 64-digit hexadecimal SHA256 values");
    }
    if (provenance.cell_bohr.size() != 9)
    {
        throw std::invalid_argument(
            "Sternheimer SIAB provenance cell_bohr must contain complete 3x3 row-major lattice vectors");
    }
    for (const double cell_value: provenance.cell_bohr)
    {
        if (!std::isfinite(cell_value))
        {
            throw std::invalid_argument("Sternheimer SIAB provenance cell_bohr values must be finite");
        }
    }
    const long double a = provenance.cell_bohr[0];
    const long double b = provenance.cell_bohr[1];
    const long double c = provenance.cell_bohr[2];
    const long double d = provenance.cell_bohr[3];
    const long double e = provenance.cell_bohr[4];
    const long double f = provenance.cell_bohr[5];
    const long double g = provenance.cell_bohr[6];
    const long double h = provenance.cell_bohr[7];
    const long double i = provenance.cell_bohr[8];
    const long double determinant = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g);
    if (!std::isfinite(determinant) || determinant == 0.0L)
    {
        throw std::invalid_argument("Sternheimer SIAB provenance cell_bohr lattice must be nonsingular");
    }
    if (!std::isfinite(provenance.ecut_ry) || provenance.ecut_ry <= 0.0)
    {
        throw std::invalid_argument("Sternheimer SIAB provenance ecut_ry must be finite and positive");
    }
    if (provenance.kernel.empty() || provenance.spin_convention.empty())
    {
        throw std::invalid_argument("Sternheimer SIAB provenance kernel and spin convention must not be empty");
    }

    const bool has_task4 = !provenance.executable_sha256.empty() || provenance.exx_pca_thr != -1.0
                           || provenance.sternheimer_nfreq != 0 || !provenance.frequency_ha.empty()
                           || !provenance.frequency_weights_ha.empty() || provenance.mpi_ranks != 0
                           || provenance.omp_threads != 0;
    if (has_task4)
    {
        if (!valid_hex(provenance.executable_sha256, 64))
        {
            throw std::invalid_argument(
                "Sternheimer SIAB production provenance requires a 64-digit executable SHA256 value");
        }
        if (!std::isfinite(provenance.exx_pca_thr) || provenance.exx_pca_thr < 0.0
            || provenance.sternheimer_nfreq <= 0 || provenance.mpi_ranks <= 0 || provenance.omp_threads <= 0
            || provenance.frequency_ha.size() != static_cast<std::size_t>(provenance.sternheimer_nfreq)
            || provenance.frequency_weights_ha.size() != static_cast<std::size_t>(provenance.sternheimer_nfreq))
        {
            throw std::invalid_argument("Sternheimer SIAB production provenance fields are incomplete or invalid");
        }
        for (std::size_t index = 0; index != provenance.frequency_ha.size(); ++index)
        {
            if (!std::isfinite(provenance.frequency_ha[index]) || provenance.frequency_ha[index] < 0.0
                || !std::isfinite(provenance.frequency_weights_ha[index])
                || provenance.frequency_weights_ha[index] < 0.0)
            {
                throw std::invalid_argument(
                    "Sternheimer SIAB production frequencies and weights must be finite and non-negative");
            }
            if (index != 0 && provenance.frequency_ha[index] <= provenance.frequency_ha[index - 1])
            {
                throw std::invalid_argument(
                    "Sternheimer SIAB production frequencies must be strictly increasing in frequency-index order");
            }
        }
    }
}

std::vector<const ReferenceRow*> validate_and_sort_rows(const std::vector<ReferenceRow>& rows,
                                                        const std::size_t n_primitive)
{
    std::vector<const ReferenceRow*> sorted_rows;
    sorted_rows.reserve(rows.size());
    for (const ReferenceRow& row: rows)
    {
        if (row.q.size() != n_primitive)
        {
            throw std::invalid_argument("Sternheimer SIAB Q row size does not match n_primitive");
        }
        if (row.occupied_state < 0 || row.auxiliary_channel < 0 || row.frequency_index < 0)
        {
            throw std::invalid_argument("Sternheimer SIAB reference indices must be non-negative");
        }
        if (!std::isfinite(row.frequency_ha) || !std::isfinite(row.occupation) || !std::isfinite(row.frequency_weight)
            || !std::isfinite(row.norm) || row.occupation < 0.0 || row.frequency_weight < 0.0 || row.norm <= 0.0)
        {
            throw std::invalid_argument("Sternheimer SIAB reference metadata is invalid");
        }
        if (!std::all_of(row.q.begin(), row.q.end(), finite_complex))
        {
            throw std::invalid_argument("Sternheimer SIAB Q contains a non-finite value");
        }
        sorted_rows.push_back(&row);
    }

    const auto key = [](const ReferenceRow* row) {
        return std::make_tuple(row->occupied_state, row->auxiliary_channel, row->frequency_index);
    };
    std::sort(sorted_rows.begin(), sorted_rows.end(), [&key](const ReferenceRow* left, const ReferenceRow* right) {
        return key(left) < key(right);
    });
    for (std::size_t i = 1; i < sorted_rows.size(); ++i)
    {
        if (key(sorted_rows[i - 1]) == key(sorted_rows[i]))
        {
            throw std::invalid_argument("Sternheimer SIAB reference row keys must be unique");
        }
    }
    return sorted_rows;
}

void validate_overlap_s(const std::vector<std::complex<double>>& overlap_s, const std::size_t n_primitive)
{
    if (n_primitive != 0 && n_primitive > static_cast<std::size_t>(-1) / n_primitive)
    {
        throw std::invalid_argument("Sternheimer SIAB primitive overlap dimensions overflow");
    }
    if (overlap_s.size() != n_primitive * n_primitive)
    {
        throw std::invalid_argument("Sternheimer SIAB S size does not match n_primitive squared");
    }
    if (!std::all_of(overlap_s.begin(), overlap_s.end(), finite_complex))
    {
        throw std::invalid_argument("Sternheimer SIAB S contains a non-finite value");
    }
    for (std::size_t row = 0; row < n_primitive; ++row)
    {
        for (std::size_t column = 0; column < n_primitive; ++column)
        {
            if (std::abs(overlap_s[row * n_primitive + column] - std::conj(overlap_s[column * n_primitive + row]))
                > hermitian_tolerance)
            {
                throw std::invalid_argument("Sternheimer SIAB S must be Hermitian within absolute tolerance 1e-10");
            }
        }
    }
}

bool round_trips_exactly(const std::string& text, const double value)
{
    std::istringstream input(text);
    input.imbue(std::locale::classic());
    double parsed = 0.0;
    input >> parsed;
    const bool exact_subnormal
        = value != 0.0 && std::abs(value) < std::numeric_limits<double>::min() && parsed == value;
    if (input.bad() || (input.fail() && !exact_subnormal))
    {
        return false;
    }
    return input.rdbuf()->in_avail() == 0 && parsed == value
           && (value != 0.0 || std::signbit(parsed) == std::signbit(value));
}

std::string normalize_exponent(const std::string& text)
{
    const std::size_t exponent_position = text.find_first_of("eE");
    if (exponent_position == std::string::npos)
    {
        return text;
    }

    std::size_t exponent_digit = exponent_position + 1;
    bool negative_exponent = false;
    if (exponent_digit < text.size() && (text[exponent_digit] == '+' || text[exponent_digit] == '-'))
    {
        negative_exponent = text[exponent_digit] == '-';
        ++exponent_digit;
    }
    while (exponent_digit + 1 < text.size() && text[exponent_digit] == '0')
    {
        ++exponent_digit;
    }
    return text.substr(0, exponent_position) + 'e' + (negative_exponent ? "-" : "") + text.substr(exponent_digit);
}

int parse_exponent(const std::string& text, const std::size_t exponent_position)
{
    std::size_t exponent_digit = exponent_position + 1;
    bool negative_exponent = false;
    if (exponent_digit < text.size() && text[exponent_digit] == '-')
    {
        negative_exponent = true;
        ++exponent_digit;
    }
    int exponent = 0;
    for (; exponent_digit < text.size(); ++exponent_digit)
    {
        exponent = exponent * 10 + (text[exponent_digit] - '0');
    }
    return negative_exponent ? -exponent : exponent;
}

std::string expand_exponent(const std::string& text)
{
    const std::size_t exponent_position = text.find('e');
    if (exponent_position == std::string::npos)
    {
        return text;
    }

    std::size_t mantissa_begin = 0;
    std::string sign;
    if (text[0] == '-' || text[0] == '+')
    {
        sign.assign(1, text[0]);
        mantissa_begin = 1;
    }
    const std::size_t decimal_position = text.find('.', mantissa_begin);
    std::string digits = text.substr(mantissa_begin, exponent_position - mantissa_begin);
    if (decimal_position != std::string::npos && decimal_position < exponent_position)
    {
        digits.erase(decimal_position - mantissa_begin, 1);
    }

    const int exponent = parse_exponent(text, exponent_position);
    const int digits_before_decimal = decimal_position == std::string::npos
                                          ? static_cast<int>(digits.size())
                                          : static_cast<int>(decimal_position - mantissa_begin);
    const int expanded_decimal_position = digits_before_decimal + exponent;
    std::string expanded = sign;
    if (expanded_decimal_position <= 0)
    {
        expanded += "0.";
        expanded.append(static_cast<std::size_t>(-expanded_decimal_position), '0');
        expanded += digits;
    }
    else if (expanded_decimal_position >= static_cast<int>(digits.size()))
    {
        expanded += digits;
        expanded.append(static_cast<std::size_t>(expanded_decimal_position - static_cast<int>(digits.size())), '0');
    }
    else
    {
        expanded += digits.substr(0, static_cast<std::size_t>(expanded_decimal_position));
        expanded += '.';
        expanded += digits.substr(static_cast<std::size_t>(expanded_decimal_position));
    }
    return expanded;
}

std::string format_exponent(const int exponent)
{
    std::ostringstream output;
    output.imbue(std::locale::classic());
    output << exponent;
    return output.str();
}

std::string decimal_to_scientific(const std::string& decimal)
{
    std::size_t digits_begin = 0;
    std::string sign;
    if (decimal[0] == '-' || decimal[0] == '+')
    {
        sign.assign(1, decimal[0]);
        digits_begin = 1;
    }
    const std::size_t decimal_position = decimal.find('.', digits_begin);
    const int digits_before_decimal = decimal_position == std::string::npos
                                          ? static_cast<int>(decimal.size() - digits_begin)
                                          : static_cast<int>(decimal_position - digits_begin);
    std::string digits = decimal.substr(digits_begin);
    if (decimal_position != std::string::npos)
    {
        digits.erase(decimal_position - digits_begin, 1);
    }
    const std::size_t first_nonzero = digits.find_first_not_of('0');
    if (first_nonzero == std::string::npos)
    {
        return sign + "0e0";
    }
    const std::size_t last_nonzero = digits.find_last_not_of('0');
    const int exponent = digits_before_decimal - static_cast<int>(first_nonzero) - 1;
    const std::string significant = digits.substr(first_nonzero, last_nonzero - first_nonzero + 1);
    std::string scientific = sign + significant.substr(0, 1);
    if (significant.size() > 1)
    {
        scientific += '.';
        scientific += significant.substr(1);
    }
    scientific += 'e';
    scientific += format_exponent(exponent);
    return scientific;
}

void add_exact_candidate(std::vector<std::string>& candidates, const std::string& candidate, const double value)
{
    std::string compatible = candidate;
    if (compatible.find_first_of(".eE") == std::string::npos)
    {
        compatible += ".0";
    }
    if (round_trips_exactly(compatible, value)
        && std::find(candidates.begin(), candidates.end(), compatible) == candidates.end())
    {
        candidates.push_back(compatible);
    }
}

bool scientific_notation(const std::string& value)
{
    return value.find('e') != std::string::npos;
}

bool better_candidate(const std::string& left, const std::string& right)
{
    if (left.size() != right.size())
    {
        return left.size() < right.size();
    }
    if (scientific_notation(left) != scientific_notation(right))
    {
        return !scientific_notation(left);
    }
    return left < right;
}

std::string small_integral_compatibility(const double value)
{
    std::ostringstream output;
    output.imbue(std::locale::classic());
    output << std::fixed << std::setprecision(0) << value << ".0";
    return output.str();
}

std::string format_double(const double value)
{
    if (std::trunc(value) == value && std::abs(value) < 100.0)
    {
        return small_integral_compatibility(value);
    }

    std::vector<std::string> candidates;
    for (int precision = 1; precision <= std::numeric_limits<double>::max_digits10; ++precision)
    {
        std::ostringstream output;
        output.imbue(std::locale::classic());
        output << std::defaultfloat << std::setprecision(precision) << value;
        if (!output.good() || !round_trips_exactly(output.str(), value))
        {
            continue;
        }
        const std::string normalized = normalize_exponent(output.str());
        const std::string decimal = expand_exponent(normalized);
        add_exact_candidate(candidates, normalized, value);
        add_exact_candidate(candidates, decimal, value);
        add_exact_candidate(candidates, decimal_to_scientific(decimal), value);
    }
    if (candidates.empty())
    {
        throw std::runtime_error("Failed to format Sternheimer SIAB floating-point value exactly");
    }
    return *std::min_element(candidates.begin(), candidates.end(), better_candidate);
}

std::string json_string(const std::string& value)
{
    static const char hexadecimal[] = "0123456789ABCDEF";
    std::string escaped;
    escaped.reserve(value.size() + 2);
    escaped += '"';
    for (const unsigned char character: value)
    {
        switch (character)
        {
        case '"':
            escaped += "\\\"";
            break;
        case '\\':
            escaped += "\\\\";
            break;
        case '\b':
            escaped += "\\b";
            break;
        case '\f':
            escaped += "\\f";
            break;
        case '\n':
            escaped += "\\n";
            break;
        case '\r':
            escaped += "\\r";
            break;
        case '\t':
            escaped += "\\t";
            break;
        default:
            if (character < 0x20)
            {
                escaped += "\\u00";
                escaped += hexadecimal[character >> 4];
                escaped += hexadecimal[character & 0x0F];
            }
            else
            {
                escaped += static_cast<char>(character);
            }
        }
    }
    escaped += '"';
    return escaped;
}

std::string provenance_json(const Provenance& provenance)
{
    std::ostringstream output;
    output.imbue(std::locale::classic());
    output << "{\"abacus_commit\":" << json_string(provenance.abacus_commit)
           << ",\"auxiliary_basis_sha256\":" << json_string(provenance.auxiliary_basis_sha256) << ",\"cell_bohr\":[";
    for (std::size_t index = 0; index < provenance.cell_bohr.size(); ++index)
    {
        if (index != 0)
        {
            output << ',';
        }
        output << format_double(provenance.cell_bohr[index]);
    }
    output << "],\"ecut_ry\":" << format_double(provenance.ecut_ry) << ",\"kernel\":" << json_string(provenance.kernel)
           << ",\"orbital_sha256\":" << json_string(provenance.orbital_sha256)
           << ",\"pseudopotential_sha256\":" << json_string(provenance.pseudopotential_sha256)
           << ",\"spin_convention\":" << json_string(provenance.spin_convention);
    if (!provenance.executable_sha256.empty())
    {
        output << ",\"executable_sha256\":" << json_string(provenance.executable_sha256)
               << ",\"exx_pca_thr\":" << format_double(provenance.exx_pca_thr) << ",\"frequency_ha\":[";
        for (std::size_t index = 0; index != provenance.frequency_ha.size(); ++index)
        {
            if (index != 0)
            {
                output << ',';
            }
            output << format_double(provenance.frequency_ha[index]);
        }
        output << "],\"frequency_weights_ha\":[";
        for (std::size_t index = 0; index != provenance.frequency_weights_ha.size(); ++index)
        {
            if (index != 0)
            {
                output << ',';
            }
            output << format_double(provenance.frequency_weights_ha[index]);
        }
        output << "],\"mpi_ranks\":" << provenance.mpi_ranks << ",\"omp_threads\":" << provenance.omp_threads
               << ",\"sternheimer_nfreq\":" << provenance.sternheimer_nfreq;
    }
    output << '}';
    if (!output.good())
    {
        throw std::runtime_error("Failed to format Sternheimer SIAB provenance JSON");
    }
    return output.str();
}

class TemporaryFile
{
  public:
    explicit TemporaryFile(const std::string& path) : path_(path), keep_(false)
    {
    }

    ~TemporaryFile()
    {
        if (!keep_)
        {
            std::remove(path_.c_str());
        }
    }

    void keep()
    {
        keep_ = true;
    }

  private:
    const std::string path_;
    bool keep_;
};

} // namespace

void write_v1(const std::string& path,
              const double grid_volume_bohr3,
              const std::vector<PrimitiveBlock>& blocks,
              const std::vector<ReferenceRow>& rows,
              const std::vector<std::complex<double>>& overlap_s,
              const Provenance& provenance)
{
    if (path.empty())
    {
        throw std::invalid_argument("Sternheimer SIAB output path must not be empty");
    }
    if (!std::isfinite(grid_volume_bohr3) || grid_volume_bohr3 <= 0.0)
    {
        throw std::invalid_argument("Sternheimer SIAB grid volume must be finite and positive");
    }

    const std::size_t n_primitive = validate_blocks(blocks);
    const std::vector<const ReferenceRow*> sorted_rows = validate_and_sort_rows(rows, n_primitive);
    validate_overlap_s(overlap_s, n_primitive);
    validate_provenance(provenance);
    const std::string json = provenance_json(provenance);

    const std::string temporary_path = path + ".tmp";
    TemporaryFile temporary_file(temporary_path);
    std::ofstream output(temporary_path.c_str(), std::ios::binary | std::ios::trunc);
    output.imbue(std::locale::classic());
    if (!output.is_open())
    {
        throw std::runtime_error("Failed to open Sternheimer SIAB temporary output: " + temporary_path);
    }

    output << "<STERNHEIMER_SIAB_HEADER>\n"
           << "format_version 1\n"
           << "n_reference " << sorted_rows.size() << "\n"
           << "n_primitive " << n_primitive << "\n"
           << "n_blocks " << blocks.size() << "\n"
           << "grid_volume_bohr3 " << format_double(grid_volume_bohr3) << "\n"
           << "</STERNHEIMER_SIAB_HEADER>\n"
           << "<PRIMITIVE_BLOCKS>\n"
           << "# element atom_index l m n_primitive offset\n";
    for (const PrimitiveBlock& block: blocks)
    {
        output << block.element << " " << block.atom_index << " " << block.l << " " << block.m << " "
               << block.n_primitive << " " << block.offset << "\n";
    }

    output << "</PRIMITIVE_BLOCKS>\n"
           << "<REFERENCE_METADATA>\n"
           << "# occupied_state auxiliary_channel frequency_ha occupation frequency_weight norm\n";
    for (const ReferenceRow* row: sorted_rows)
    {
        output << row->occupied_state << " " << row->auxiliary_channel << " " << format_double(row->frequency_ha) << " "
               << format_double(row->occupation) << " " << format_double(row->frequency_weight) << " "
               << format_double(row->norm) << "\n";
    }

    output << "</REFERENCE_METADATA>\n"
           << "<OVERLAP_Q>\n"
           << "# row-major <Y_rho|B_e>, real imag\n";
    for (const ReferenceRow* row: sorted_rows)
    {
        for (const std::complex<double>& value: row->q)
        {
            output << format_double(value.real()) << " " << format_double(value.imag()) << "\n";
        }
    }

    output << "</OVERLAP_Q>\n"
           << "<OVERLAP_S>\n"
           << "# row-major <B_e|B_ep>, real imag\n";
    for (const std::complex<double>& value: overlap_s)
    {
        output << format_double(value.real()) << " " << format_double(value.imag()) << "\n";
    }
    output << "</OVERLAP_S>\n"
           << "<PROVENANCE_JSON>\n"
           << json << "\n"
           << "</PROVENANCE_JSON>\n";

    output.flush();
    if (!output.good())
    {
        throw std::runtime_error("Failed to flush Sternheimer SIAB temporary output: " + temporary_path);
    }
    output.close();
    if (output.fail())
    {
        throw std::runtime_error("Failed to close Sternheimer SIAB temporary output: " + temporary_path);
    }
    if (std::rename(temporary_path.c_str(), path.c_str()) != 0)
    {
        const std::string reason = std::strerror(errno);
        throw std::runtime_error("Failed to replace Sternheimer SIAB output " + path + ": " + reason);
    }
    temporary_file.keep();
}

} // namespace sternheimer_siab
} // namespace module_ri
