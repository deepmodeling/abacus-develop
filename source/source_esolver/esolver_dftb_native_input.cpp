#include "esolver_dftb_native_input.h"

#include "source_cell/unitcell.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace ModuleESolver
{
namespace NativeDftb
{
namespace
{
struct KptMesh
{
    std::array<int, 3> dimensions;
    std::array<double, 3> offset;
    bool monkhorst_pack;
};

struct BandPathHeader
{
    int segments;
    int intervals;
};

struct BandSegment
{
    std::string start_label;
    std::string end_label;
    std::array<double, 3> start;
    std::array<double, 3> end;
};

std::string lowercase(std::string value)
{
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char character) {
        return static_cast<char>(std::tolower(character));
    });
    return value;
}

template <typename T>
T parse_kpt_value(const std::vector<std::string>& tokens, std::size_t* position, const std::string& description)
{
    if (*position >= tokens.size()) throw std::runtime_error("Missing " + description + " in ABACUS KPT file");
    std::istringstream parser(tokens[(*position)++]);
    T value{};
    std::string extra;
    if (!(parser >> value) || (parser >> extra))
        throw std::runtime_error("Invalid " + description + " in ABACUS KPT file");
    return value;
}

bool parse_path_and_species_option(const std::string& key,
                                   std::istringstream& row,
                                   NativeDftbConfig& config,
                                   const std::string& filename,
                                   int line_number)
{
    if (key == "skf_dir")
    {
        if (!(row >> config.skf_directory)) throw std::runtime_error("Missing skf_dir value in " + filename);
    }
    else if (key == "band_path_file")
    {
        if (!(row >> config.band_path_file)) throw std::runtime_error("Missing band_path_file value");
    }
    else if (key == "hubbard_deriv")
    {
        std::string species;
        double value = 0.0;
        if (!(row >> species >> value))
            throw std::runtime_error("Expected 'hubbard_deriv Element value' at line " + std::to_string(line_number));
        if (!config.hubbard_derivatives.insert(std::make_pair(species, value)).second)
            throw std::runtime_error("Duplicate Hubbard derivative for species " + species);
    }
    else
    {
        return false;
    }
    return true;
}

bool parse_thermal_option(const std::string& key, std::istringstream& row, NativeDftbConfig& config)
{
    if (key == "temperature_kelvin")
    {
        if (!(row >> config.temperature_kelvin)) throw std::runtime_error("Invalid temperature_kelvin value");
    }
    else if (key == "scc_tolerance")
    {
        if (!(row >> config.scc_tolerance)) throw std::runtime_error("Invalid scc_tolerance value");
    }
    else
    {
        return false;
    }
    return true;
}

bool parse_mixing_option(const std::string& key, std::istringstream& row, NativeDftbConfig& config)
{
    if (key == "max_scc_iterations")
    {
        if (!(row >> config.maximum_scc_iterations)) throw std::runtime_error("Invalid max_scc_iterations value");
    }
    else if (key == "mixing_parameter")
    {
        if (!(row >> config.mixing_parameter)) throw std::runtime_error("Invalid mixing_parameter value");
    }
    else if (key == "mixing_method")
    {
        if (!(row >> config.mixing_method)) throw std::runtime_error("Missing mixing_method value");
        config.mixing_method = lowercase(config.mixing_method);
    }
    else if (key == "mixing_history")
    {
        if (!(row >> config.mixing_history)) throw std::runtime_error("Invalid mixing_history value");
    }
    else
    {
        return false;
    }
    return true;
}

bool parse_scc_option(const std::string& key, std::istringstream& row, NativeDftbConfig& config)
{
    return parse_thermal_option(key, row, config) || parse_mixing_option(key, row, config);
}

bool parse_broyden_option(const std::string& key, std::istringstream& row, NativeDftbConfig& config)
{
    if (key == "broyden_inverse_jacobi_weight")
    {
        if (!(row >> config.broyden_inverse_jacobi_weight))
            throw std::runtime_error("Invalid broyden_inverse_jacobi_weight value");
    }
    else if (key == "broyden_minimal_weight")
    {
        if (!(row >> config.broyden_minimal_weight)) throw std::runtime_error("Invalid broyden_minimal_weight value");
    }
    else if (key == "broyden_maximal_weight")
    {
        if (!(row >> config.broyden_maximal_weight)) throw std::runtime_error("Invalid broyden_maximal_weight value");
    }
    else if (key == "broyden_weight_factor")
    {
        if (!(row >> config.broyden_weight_factor)) throw std::runtime_error("Invalid broyden_weight_factor value");
    }
    else
    {
        return false;
    }
    return true;
}

bool parse_output_precision_option(const std::string& key, std::istringstream& row, NativeDftbConfig& config)
{
    if (key != "output_precision") return false;
    if (!(row >> config.output_precision)) throw std::runtime_error("Invalid output_precision value");
    return true;
}

bool parse_third_order_option(const std::string& key, std::istringstream& row, NativeDftbConfig& config)
{
    if (key != "third_order") return false;
    std::string value;
    if (!(row >> value)) throw std::runtime_error("Invalid third_order value");
    value = lowercase(value);
    if (value == "yes" || value == "true" || value == "1" || value == "on")
    {
        config.third_order = true;
    }
    else if (value == "no" || value == "false" || value == "0" || value == "off")
    {
        config.third_order = false;
    }
    else
    {
        throw std::runtime_error("third_order must be yes/no, true/false, or 1/0");
    }
    return true;
}

bool parse_output_option(const std::string& key, std::istringstream& row, NativeDftbConfig& config)
{
    return parse_output_precision_option(key, row, config) || parse_third_order_option(key, row, config);
}

void parse_config_option(const std::string& key,
                         std::istringstream& row,
                         NativeDftbConfig& config,
                         const std::string& filename,
                         int line_number)
{
    const bool parsed = parse_path_and_species_option(key, row, config, filename, line_number)
                        || parse_scc_option(key, row, config) || parse_broyden_option(key, row, config)
                        || parse_output_option(key, row, config);
    if (!parsed)
        throw std::runtime_error("Unknown native DFTB configuration key at line "
                                 + std::to_string(line_number) + ": " + key);
}

void validate_scc_controls(const NativeDftbConfig& config)
{
    if (!(config.temperature_kelvin >= 0.0) || !(config.scc_tolerance > 0.0)
        || config.maximum_scc_iterations <= 0
        || !(config.mixing_parameter > 0.0 && config.mixing_parameter <= 1.0)
        || (config.mixing_method != "linear" && config.mixing_method != "pulay"
            && config.mixing_method != "broyden")
        || config.mixing_history < 2 || config.mixing_history > 20)
        throw std::runtime_error("Invalid native DFTB temperature or SCC controls");
}

void validate_broyden_controls(const NativeDftbConfig& config)
{
    if (!(config.broyden_inverse_jacobi_weight > 0.0) || !(config.broyden_minimal_weight > 0.0)
        || !(config.broyden_maximal_weight >= config.broyden_minimal_weight)
        || !(config.broyden_weight_factor > 0.0) || !std::isfinite(config.broyden_inverse_jacobi_weight)
        || !std::isfinite(config.broyden_minimal_weight) || !std::isfinite(config.broyden_maximal_weight)
        || !std::isfinite(config.broyden_weight_factor))
        throw std::runtime_error("Invalid native DFTB Broyden controls");
}

void validate_config(const NativeDftbConfig& config)
{
    if (config.skf_directory.empty()) throw std::runtime_error("Native DFTB config must define skf_dir");
    validate_scc_controls(config);
    validate_broyden_controls(config);
    if (config.output_precision < 1 || config.output_precision > 17)
        throw std::runtime_error("Invalid native DFTB output precision");
}

std::vector<std::string> read_kpt_tokens(const std::string& filename)
{
    std::ifstream input(filename.c_str());
    if (!input) throw std::runtime_error("Cannot open ABACUS KPT file: " + filename);
    std::vector<std::string> tokens;
    std::string line;
    while (std::getline(input, line))
    {
        const std::size_t hash_comment = line.find('#');
        const std::size_t slash_comment = line.find("//");
        std::size_t comment = std::string::npos;
        if (hash_comment != std::string::npos) comment = hash_comment;
        if (slash_comment != std::string::npos) comment = std::min(comment, slash_comment);
        if (comment != std::string::npos) line.erase(comment);
        std::istringstream row(line);
        std::string token;
        while (row >> token) tokens.push_back(token);
    }
    return tokens;
}

std::size_t find_kpt_header(const std::vector<std::string>& tokens)
{
    std::size_t header = 0;
    while (header < tokens.size() && tokens[header] != "K_POINTS" && tokens[header] != "KPOINTS"
           && tokens[header] != "K")
        ++header;
    if (header == tokens.size()) throw std::runtime_error("Native DFTB expects a valid ABACUS KPT file");
    return header;
}

void validate_kpt_tail(const std::vector<std::string>& tokens, std::size_t position)
{
    if (position != tokens.size()) throw std::runtime_error("Unexpected trailing data in ABACUS KPT file");
}

void validate_mesh(const std::array<int, 3>& dimensions)
{
    if (dimensions[0] <= 0 || dimensions[1] <= 0 || dimensions[2] <= 0)
        throw std::runtime_error("KPT mesh dimensions must be positive");
    long long mesh_size = dimensions[0];
    if (dimensions[1] > 100000 / mesh_size) throw std::runtime_error("ABACUS KPT mesh exceeds the 100000-point limit");
    mesh_size *= dimensions[1];
    if (dimensions[2] > 100000 / mesh_size) throw std::runtime_error("ABACUS KPT mesh exceeds the 100000-point limit");
}

KptMesh read_kpt_mesh(const std::vector<std::string>& tokens, std::size_t* position)
{
    const std::string mode = lowercase(parse_kpt_value<std::string>(tokens, position, "mesh mode"));
    KptMesh mesh;
    mesh.offset = {{0.0, 0.0, 0.0}};
    for (int axis = 0; axis < 3; ++axis)
        mesh.dimensions[axis] = parse_kpt_value<int>(tokens, position, "mesh dimension");
    if (*position != tokens.size() && tokens.size() - *position != 3)
        throw std::runtime_error("An automatic ABACUS KPT mesh accepts either zero or three offsets");
    if (*position < tokens.size())
        for (int axis = 0; axis < 3; ++axis)
            mesh.offset[axis] = parse_kpt_value<double>(tokens, position, "mesh offset");
    for (int axis = 0; axis < 3; ++axis)
        if (!std::isfinite(mesh.offset[axis])) throw std::runtime_error("KPT mesh offsets must be finite");
    if (mode != "gamma" && mode != "monkhorst-pack" && mode != "mp")
        throw std::runtime_error("Native DFTB supports Gamma or Monkhorst-Pack KPT meshes");
    mesh.monkhorst_pack = mode != "gamma";
    validate_mesh(mesh.dimensions);
    return mesh;
}

std::vector<ModuleDFTB::DftbWeightedKPoint> make_mesh_kpoints(const KptMesh& mesh)
{
    const double weight = 1.0 / static_cast<double>(mesh.dimensions[0] * mesh.dimensions[1] * mesh.dimensions[2]);
    std::vector<ModuleDFTB::DftbWeightedKPoint> result;
    for (int ix = 1; ix <= mesh.dimensions[0]; ++ix)
    {
        for (int iy = 1; iy <= mesh.dimensions[1]; ++iy)
        {
            for (int iz = 1; iz <= mesh.dimensions[2]; ++iz)
            {
                ModuleDFTB::DftbWeightedKPoint point;
                const int indices[3] = {ix, iy, iz};
                for (int axis = 0; axis < 3; ++axis)
                {
                    const double n = static_cast<double>(indices[axis]);
                    const double dim = static_cast<double>(mesh.dimensions[axis]);
                    if (mesh.monkhorst_pack)
                        point.fractional[axis] = (mesh.offset[axis] + 2.0 * n - dim - 1.0) / (2.0 * dim);
                    else
                        point.fractional[axis] = (mesh.offset[axis] + n - 1.0) / dim;
                }
                point.weight = weight;
                result.push_back(point);
            }
        }
    }
    return result;
}

bool is_line_kpt_mode(const std::string& mode)
{
    return mode == "line" || mode == "line_direct" || mode == "l" || mode == "line_cartesian";
}

bool read_direct_coordinate_mode(const std::string& mode)
{
    if (mode == "direct" || mode == "d") return true;
    if (mode == "cartesian" || mode == "c") return false;
    if (is_line_kpt_mode(mode))
        throw std::runtime_error(
            "Native DFTB SCC needs weighted integration k-points; ABACUS line-mode points are band paths.\n"
            "Use a Gamma/MP mesh or weighted Direct/Cartesian list in KPT, and set band_path_file in dftb_native.in.");
    throw std::runtime_error("Explicit native DFTB KPT coordinates must be Direct/D or Cartesian/C");
}

std::array<double, 3> read_explicit_coordinate(const std::vector<std::string>& tokens, std::size_t* position)
{
    std::array<double, 3> coordinate;
    for (int axis = 0; axis < 3; ++axis)
    {
        coordinate[axis] = parse_kpt_value<double>(tokens, position, "explicit k-point coordinate");
        if (!std::isfinite(coordinate[axis])) throw std::runtime_error("KPT coordinates must be finite");
    }
    return coordinate;
}

ModuleDFTB::DftbWeightedKPoint read_explicit_kpoint(const std::vector<std::string>& tokens,
                                                    std::size_t* position,
                                                    bool direct,
                                                    const UnitCell& ucell)
{
    ModuleDFTB::DftbWeightedKPoint point;
    const std::array<double, 3> coordinate = read_explicit_coordinate(tokens, position);
    point.weight = parse_kpt_value<double>(tokens, position, "explicit k-point weight");
    if (!(point.weight >= 0.0) || !std::isfinite(point.weight))
        throw std::runtime_error("Malformed explicit k-point list in KPT file");
    if (direct)
    {
        point.fractional = coordinate;
    }
    else
    {
        // ABACUS Cartesian KPT coordinates are expressed in 2pi/a units.
        // Inverting k_cart = k_direct * G gives k_direct = k_cart * latvec^T.
        point.fractional[0] = coordinate[0] * ucell.latvec.e11 + coordinate[1] * ucell.latvec.e12
                              + coordinate[2] * ucell.latvec.e13;
        point.fractional[1] = coordinate[0] * ucell.latvec.e21 + coordinate[1] * ucell.latvec.e22
                              + coordinate[2] * ucell.latvec.e23;
        point.fractional[2] = coordinate[0] * ucell.latvec.e31 + coordinate[1] * ucell.latvec.e32
                              + coordinate[2] * ucell.latvec.e33;
    }
    return point;
}

void normalize_kpoint_weights(std::vector<ModuleDFTB::DftbWeightedKPoint>* points, double weight_sum)
{
    if (!(weight_sum > 0.0)) throw std::runtime_error("KPT weights must have a positive sum");
    for (std::size_t i = 0; i < points->size(); ++i) (*points)[i].weight /= weight_sum;
}

std::vector<ModuleDFTB::DftbWeightedKPoint> read_explicit_kpoints(const std::vector<std::string>& tokens,
                                                                   std::size_t* position,
                                                                   int count,
                                                                   const UnitCell& ucell)
{
    const std::string mode = lowercase(parse_kpt_value<std::string>(tokens, position, "explicit coordinate mode"));
    const bool direct = read_direct_coordinate_mode(mode);
    std::vector<ModuleDFTB::DftbWeightedKPoint> result;
    double weight_sum = 0.0;
    for (int i = 0; i < count; ++i)
    {
        const ModuleDFTB::DftbWeightedKPoint point = read_explicit_kpoint(tokens, position, direct, ucell);
        weight_sum += point.weight;
        result.push_back(point);
    }
    validate_kpt_tail(tokens, *position);
    normalize_kpoint_weights(&result, weight_sum);
    return result;
}

std::vector<std::string> read_band_path_lines(const std::string& filename)
{
    std::ifstream input(filename.c_str());
    if (!input) throw std::runtime_error("Cannot open native DFTB band path: " + filename);
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(input, line))
    {
        const std::size_t comment = line.find('#');
        if (comment != std::string::npos) line.erase(comment);
        std::istringstream row(line);
        std::string first;
        if (row >> first) lines.push_back(line);
    }
    if (lines.empty()) throw std::runtime_error("Empty native DFTB band path: " + filename);
    return lines;
}

BandPathHeader read_band_path_header(const std::vector<std::string>& lines)
{
    BandPathHeader header;
    std::string extra;
    std::istringstream row(lines.front());
    if (!(row >> header.segments >> header.intervals) || (row >> extra) || header.segments <= 0
        || header.intervals <= 0 || lines.size() != static_cast<std::size_t>(header.segments + 1))
        throw std::runtime_error("Band path header must be '<segments> <intervals>' followed by one line per segment");
    return header;
}

BandSegment read_band_segment(const std::string& line,
                              int segment_index,
                              const std::array<double, 3>& previous_end)
{
    BandSegment segment;
    std::string extra;
    std::istringstream row(line);
    if (!(row >> segment.start_label >> segment.start[0] >> segment.start[1] >> segment.start[2]
              >> segment.end_label >> segment.end[0] >> segment.end[1] >> segment.end[2])
        || (row >> extra))
        throw std::runtime_error("Each band path segment must be '<label> kx ky kz <label> kx ky kz'");
    for (int axis = 0; axis < 3; ++axis)
    {
        if (!std::isfinite(segment.start[axis]) || !std::isfinite(segment.end[axis]))
            throw std::runtime_error("Band path coordinates must be finite");
        if (segment_index > 0 && std::abs(segment.start[axis] - previous_end[axis]) > 1.0e-10)
            throw std::runtime_error("Band path segments must be connected in the listed order");
    }
    return segment;
}

void append_band_segment(std::vector<ModuleDFTB::DftbBandKPoint>* points,
                         const BandSegment& segment,
                         int segment_index,
                         int intervals)
{
    if (segment_index == 0)
    {
        ModuleDFTB::DftbBandKPoint point;
        point.fractional = segment.start;
        point.label = segment.start_label;
        points->push_back(point);
    }
    for (int step = 1; step <= intervals; ++step)
    {
        const double t = static_cast<double>(step) / static_cast<double>(intervals);
        ModuleDFTB::DftbBandKPoint point;
        for (int axis = 0; axis < 3; ++axis)
            point.fractional[axis] = segment.start[axis] + t * (segment.end[axis] - segment.start[axis]);
        if (step == intervals) point.label = segment.end_label;
        points->push_back(point);
    }
}
} // namespace

NativeDftbConfig read_native_config(const std::string& filename)
{
    std::ifstream input(filename.c_str());
    if (!input) throw std::runtime_error("Cannot open native DFTB configuration: " + filename);
    NativeDftbConfig config;
    std::string line;
    int line_number = 0;
    while (std::getline(input, line))
    {
        ++line_number;
        const std::size_t comment = line.find('#');
        if (comment != std::string::npos) line.erase(comment);
        std::istringstream row(line);
        std::string key;
        if (!(row >> key)) continue;
        parse_config_option(key, row, config, filename, line_number);
    }
    validate_config(config);
    return config;
}

std::vector<ModuleDFTB::DftbWeightedKPoint> read_native_kpoints(const std::string& filename, const UnitCell& ucell)
{
    const std::vector<std::string> tokens = read_kpt_tokens(filename);
    const std::size_t header = find_kpt_header(tokens);
    std::size_t position = header + 1;
    const int count = parse_kpt_value<int>(tokens, &position, "k-point count");
    if (count < 0 || count > 100000)
        throw std::runtime_error("ABACUS KPT k-point count is outside the supported range");
    if (count == 0)
    {
        const KptMesh mesh = read_kpt_mesh(tokens, &position);
        validate_kpt_tail(tokens, position);
        return make_mesh_kpoints(mesh);
    }
    return read_explicit_kpoints(tokens, &position, count, ucell);
}

std::vector<ModuleDFTB::DftbBandKPoint> read_native_band_path(const std::string& filename)
{
    const std::vector<std::string> lines = read_band_path_lines(filename);
    const BandPathHeader header = read_band_path_header(lines);
    std::vector<ModuleDFTB::DftbBandKPoint> points;
    std::array<double, 3> previous_end{{0.0, 0.0, 0.0}};
    for (int segment_index = 0; segment_index < header.segments; ++segment_index)
    {
        const BandSegment segment = read_band_segment(lines[static_cast<std::size_t>(segment_index + 1)],
                                                       segment_index,
                                                       previous_end);
        append_band_segment(&points, segment, segment_index, header.intervals);
        previous_end = segment.end;
    }
    return points;
}
} // namespace NativeDftb
} // namespace ModuleESolver
