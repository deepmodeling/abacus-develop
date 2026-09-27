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

}

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
            if (!(row >> species >> value)) throw std::runtime_error("Expected 'hubbard_deriv Element value' at line " + std::to_string(line_number));
            if (!config.hubbard_derivatives.insert(std::make_pair(species, value)).second)
                throw std::runtime_error("Duplicate Hubbard derivative for species " + species);
        }
        else if (key == "temperature_kelvin")
        {
            if (!(row >> config.temperature_kelvin)) throw std::runtime_error("Invalid temperature_kelvin value");
        }
        else if (key == "scc_tolerance")
        {
            if (!(row >> config.scc_tolerance)) throw std::runtime_error("Invalid scc_tolerance value");
        }
        else if (key == "max_scc_iterations")
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
        else if (key == "broyden_inverse_jacobi_weight")
        {
            if (!(row >> config.broyden_inverse_jacobi_weight)) throw std::runtime_error("Invalid broyden_inverse_jacobi_weight value");
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
        else if (key == "output_precision")
        {
            if (!(row >> config.output_precision)) throw std::runtime_error("Invalid output_precision value");
        }
        else if (key == "third_order")
        {
            std::string value;
            if (!(row >> value)) throw std::runtime_error("Invalid third_order value");
            value = lowercase(value);
            if (value == "yes" || value == "true" || value == "1" || value == "on") config.third_order = true;
            else if (value == "no" || value == "false" || value == "0" || value == "off") config.third_order = false;
            else throw std::runtime_error("third_order must be yes/no, true/false, or 1/0");
        }
        else
        {
            throw std::runtime_error("Unknown native DFTB configuration key at line " + std::to_string(line_number) + ": " + key);
        }
    }
    if (config.skf_directory.empty()) throw std::runtime_error("Native DFTB config must define skf_dir");
    if (!(config.temperature_kelvin >= 0.0) || !(config.scc_tolerance > 0.0)
        || config.maximum_scc_iterations <= 0 || !(config.mixing_parameter > 0.0 && config.mixing_parameter <= 1.0)
        || (config.mixing_method != "linear" && config.mixing_method != "pulay"
            && config.mixing_method != "broyden")
        || config.mixing_history < 2 || config.mixing_history > 20
        || !(config.broyden_inverse_jacobi_weight > 0.0)
        || !(config.broyden_minimal_weight > 0.0)
        || !(config.broyden_maximal_weight >= config.broyden_minimal_weight)
        || !(config.broyden_weight_factor > 0.0)
        || !std::isfinite(config.broyden_inverse_jacobi_weight)
        || !std::isfinite(config.broyden_minimal_weight)
        || !std::isfinite(config.broyden_maximal_weight)
        || !std::isfinite(config.broyden_weight_factor)
        || config.output_precision < 1 || config.output_precision > 17)
        throw std::runtime_error("Invalid native DFTB temperature or SCC controls");
    return config;
}

std::vector<ModuleDFTB::DftbWeightedKPoint> read_native_kpoints(const std::string& filename, const UnitCell& ucell)
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
    std::size_t header = 0;
    while (header < tokens.size() && tokens[header] != "K_POINTS" && tokens[header] != "KPOINTS" && tokens[header] != "K")
        ++header;
    if (header == tokens.size())
        throw std::runtime_error("Native DFTB expects a valid ABACUS KPT file");
    std::size_t position = header + 1;
    const int count = parse_kpt_value<int>(tokens, &position, "k-point count");
    if (count < 0 || count > 100000) throw std::runtime_error("ABACUS KPT k-point count is outside the supported range");
    std::vector<ModuleDFTB::DftbWeightedKPoint> result;
    if (count == 0)
    {
        const std::string mode = lowercase(parse_kpt_value<std::string>(tokens, &position, "mesh mode"));
        std::array<int, 3> mesh{};
        std::array<double, 3> offset{{0.0, 0.0, 0.0}};
        for (int axis = 0; axis < 3; ++axis)
            mesh[axis] = parse_kpt_value<int>(tokens, &position, "mesh dimension");
        if (position != tokens.size() && tokens.size() - position != 3)
            throw std::runtime_error("An automatic ABACUS KPT mesh accepts either zero or three offsets");
        if (position < tokens.size())
            for (int axis = 0; axis < 3; ++axis)
                offset[axis] = parse_kpt_value<double>(tokens, &position, "mesh offset");
        for (const double value : offset)
            if (!std::isfinite(value)) throw std::runtime_error("KPT mesh offsets must be finite");
        if (mode != "gamma" && mode != "monkhorst-pack" && mode != "mp")
            throw std::runtime_error("Native DFTB supports Gamma or Monkhorst-Pack KPT meshes");
        if (mesh[0] <= 0 || mesh[1] <= 0 || mesh[2] <= 0)
            throw std::runtime_error("KPT mesh dimensions must be positive");
        long long mesh_size = mesh[0];
        if (mesh[1] > 100000 / mesh_size) throw std::runtime_error("ABACUS KPT mesh exceeds the 100000-point limit");
        mesh_size *= mesh[1];
        if (mesh[2] > 100000 / mesh_size) throw std::runtime_error("ABACUS KPT mesh exceeds the 100000-point limit");
        mesh_size *= mesh[2];
        const bool monkhorst_pack = mode != "gamma";
        const double weight = 1.0 / static_cast<double>(mesh[0] * mesh[1] * mesh[2]);
        for (int ix = 1; ix <= mesh[0]; ++ix)
        {
            for (int iy = 1; iy <= mesh[1]; ++iy)
            {
                for (int iz = 1; iz <= mesh[2]; ++iz)
                {
                    ModuleDFTB::DftbWeightedKPoint point;
                    const int indices[3] = {ix, iy, iz};
                    for (int axis = 0; axis < 3; ++axis)
                    {
                        const double n = static_cast<double>(indices[axis]);
                        const double dim = static_cast<double>(mesh[axis]);
                        point.fractional[axis] = monkhorst_pack
                            ? (offset[axis] + 2.0 * n - dim - 1.0) / (2.0 * dim)
                            : (offset[axis] + n - 1.0) / dim;
                    }
                    point.weight = weight;
                    result.push_back(point);
                }
            }
        }
    }
    else
    {
        const std::string mode = lowercase(parse_kpt_value<std::string>(tokens, &position, "explicit coordinate mode"));
        const bool cartesian = mode == "cartesian" || mode == "c";
        const bool direct = mode == "direct" || mode == "d";
        if (!direct && !cartesian)
        {
            if (mode == "line" || mode == "line_direct" || mode == "l" || mode == "line_cartesian")
                throw std::runtime_error("Native DFTB SCC needs weighted integration k-points; ABACUS line-mode points are band paths.\n"
                                         "Use a Gamma/MP mesh or weighted Direct/Cartesian list in KPT, and set band_path_file in dftb_native.in.");
            throw std::runtime_error("Explicit native DFTB KPT coordinates must be Direct/D or Cartesian/C");
        }
        double weight_sum = 0.0;
        for (int i = 0; i < count; ++i)
        {
            ModuleDFTB::DftbWeightedKPoint point;
            std::array<double, 3> coordinate{};
            for (int axis = 0; axis < 3; ++axis)
                coordinate[axis] = parse_kpt_value<double>(tokens, &position, "explicit k-point coordinate");
            point.weight = parse_kpt_value<double>(tokens, &position, "explicit k-point weight");
            if (!(point.weight >= 0.0) || !std::isfinite(point.weight))
                throw std::runtime_error("Malformed explicit k-point list in KPT file");
            for (int axis = 0; axis < 3; ++axis)
                if (!std::isfinite(coordinate[axis])) throw std::runtime_error("KPT coordinates must be finite");
            if (direct)
            {
                point.fractional = coordinate;
            }
            else
            {
                // ABACUS Cartesian KPT coordinates are expressed in 2pi/a units.
                // Inverting its direct-to-Cartesian relation k_cart = k_direct * G
                // gives k_direct = k_cart * latvec^T (latvec is dimensionless here).
                point.fractional[0] = coordinate[0] * ucell.latvec.e11 + coordinate[1] * ucell.latvec.e12
                                      + coordinate[2] * ucell.latvec.e13;
                point.fractional[1] = coordinate[0] * ucell.latvec.e21 + coordinate[1] * ucell.latvec.e22
                                      + coordinate[2] * ucell.latvec.e23;
                point.fractional[2] = coordinate[0] * ucell.latvec.e31 + coordinate[1] * ucell.latvec.e32
                                      + coordinate[2] * ucell.latvec.e33;
            }
            weight_sum += point.weight;
            result.push_back(point);
        }
        if (position != tokens.size()) throw std::runtime_error("Unexpected trailing data in ABACUS KPT file");
        if (!(weight_sum > 0.0)) throw std::runtime_error("KPT weights must have a positive sum");
        for (std::size_t i = 0; i < result.size(); ++i) result[i].weight /= weight_sum;
    }
    if (position != tokens.size()) throw std::runtime_error("Unexpected trailing data in ABACUS KPT file");
    return result;
}

std::vector<ModuleDFTB::DftbBandKPoint> read_native_band_path(const std::string& filename)
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
    int segments = 0;
    int intervals = 0;
    std::istringstream header(lines.front());
    std::string extra;
    if (!(header >> segments >> intervals) || (header >> extra) || segments <= 0 || intervals <= 0
        || lines.size() != static_cast<std::size_t>(segments + 1))
        throw std::runtime_error("Band path header must be '<segments> <intervals>' followed by one line per segment");

    std::vector<ModuleDFTB::DftbBandKPoint> points;
    std::array<double, 3> previous_end{{0.0, 0.0, 0.0}};
    for (int segment = 0; segment < segments; ++segment)
    {
        std::istringstream row(lines[static_cast<std::size_t>(segment + 1)]);
        std::string start_label;
        std::string end_label;
        std::array<double, 3> start{};
        std::array<double, 3> end{};
        if (!(row >> start_label >> start[0] >> start[1] >> start[2]
                  >> end_label >> end[0] >> end[1] >> end[2]) || (row >> extra))
            throw std::runtime_error("Each band path segment must be '<label> kx ky kz <label> kx ky kz'");
        for (int axis = 0; axis < 3; ++axis)
        {
            if (!std::isfinite(start[axis]) || !std::isfinite(end[axis]))
                throw std::runtime_error("Band path coordinates must be finite");
            if (segment > 0 && std::abs(start[axis] - previous_end[axis]) > 1.0e-10)
                throw std::runtime_error("Band path segments must be connected in the listed order");
        }
        if (segment == 0)
        {
            ModuleDFTB::DftbBandKPoint point;
            point.fractional = start;
            point.label = start_label;
            points.push_back(point);
        }
        for (int step = 1; step <= intervals; ++step)
        {
            const double t = static_cast<double>(step) / static_cast<double>(intervals);
            ModuleDFTB::DftbBandKPoint point;
            for (int axis = 0; axis < 3; ++axis) point.fractional[axis] = start[axis] + t * (end[axis] - start[axis]);
            if (step == intervals) point.label = end_label;
            points.push_back(point);
        }
        previous_end = end;
    }
    return points;
}
} // namespace NativeDftb
} // namespace ModuleESolver
