#include "esolver_dftb_native.h"

#include "source_base/constants.h"
#include "source_base/global_variable.h"
#include "source_base/parallel_common.h"
#include "source_base/parallel_reduce.h"
#include "source_base/tool_quit.h"
#include "source_cell/unitcell.h"
#include "source_io/module_output/output_log.h"
#include "source_io/module_parameter/input_parameter.h"
#include "source_io/module_parameter/parameter.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace ModuleESolver
{
namespace
{
// Use the precise conversion used by DFTB+ / GEN geometry handling. ABACUS's
// legacy ANGSTROM_AU constant is rounded to 1.889727 and shifts short-range
// pair energies enough to be visible in this regression.
constexpr double angstrom_to_bohr = 1.8897261254578281;
const double kelvin_to_hartree = 3.166811563e-6;

struct NativeDftbConfig
{
    std::string skf_directory;
    std::string band_path_file;
    std::map<std::string, double> hubbard_derivatives;
    double temperature_kelvin = 0.0;
    double scc_tolerance = 1.0e-6;
    int maximum_scc_iterations = 200;
    double mixing_parameter = 0.2;
    std::string mixing_method = "linear";
    int mixing_history = 6;
    int output_precision = 12;
    bool third_order = true;
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
        || (config.mixing_method != "linear" && config.mixing_method != "pulay")
        || config.mixing_history < 2 || config.mixing_history > 20
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
} // namespace

void ESolver_DFTBNative::before_all_runners(BaseCell& cell, const Input_para& inp)
{
    this->inp_ = &inp;
    if (cell.kind() != BaseCell::Kind::unit_cell)
        ModuleBase::WARNING_QUIT("ESolver_DFTBNative", "Native DFTB requires a periodic UnitCell.");
    const UnitCell& ucell = static_cast<const UnitCell&>(cell);
    if (inp.basis_type != "dftb")
        ModuleBase::WARNING_QUIT("ESolver_DFTBNative", "Set basis_type=dftb when using esolver_type=dftbnative.");
    if (inp.calculation != "scf")
        ModuleBase::WARNING_QUIT("ESolver_DFTBNative", "Native DFTB currently supports calculation=scf only; use band_path_file for frozen-SCC bands.");
    if (inp.nspin != 1 || inp.noncolin)
        ModuleBase::WARNING_QUIT("ESolver_DFTBNative", "Native DFTB currently supports spin-degenerate calculations only.");

    const int rank = Parallel_Common::get_rank();
    std::ostream& running_log = GlobalV::ofs_running;
    int local_status = 0;
    char error_message[1024] = {};
    try
    {
        this->load_model(ucell, inp);
    }
    catch (const std::exception& error)
    {
        local_status = 1;
        std::strncpy(error_message, error.what(), sizeof(error_message) - 1);
    }
    int all_ranks_succeeded = local_status == 0 ? 1 : 0;
    Parallel_Reduce::reduce_min(all_ranks_succeeded);
    const int global_status = all_ranks_succeeded == 0 ? 1 : 0;
    if (global_status != 0)
    {
        int error_rank = local_status != 0 ? rank : std::numeric_limits<int>::max();
        Parallel_Reduce::reduce_min(error_rank);
        std::vector<int> error_codes(sizeof(error_message), 0);
        if (rank == error_rank)
            for (std::size_t i = 0; i < sizeof(error_message); ++i)
                error_codes[i] = static_cast<unsigned char>(error_message[i]);
        Parallel_Reduce::reduce_all(error_codes.data(), static_cast<int>(error_codes.size()));
        for (std::size_t i = 0; i < sizeof(error_message); ++i)
            error_message[i] = static_cast<char>(error_codes[i]);
        ModuleBase::WARNING_QUIT("ESolver_DFTBNative", error_message);
    }
    if (rank == 0)
    {
        running_log << " Native periodic DFTB model loaded from " << inp.dftb_native_input
                    << " (no DFTB+ runtime, UPF, or ABACUS orbital files)." << std::endl;
    }
}
void ESolver_DFTBNative::load_model(const UnitCell& ucell, const Input_para& inp)
{
    const NativeDftbConfig config = read_native_config(inp.dftb_native_input);
    const std::size_t ntype = static_cast<std::size_t>(ucell.ntype);
    if (ntype == 0 || ucell.nat == 0) throw std::runtime_error("Native DFTB requires at least one species and atom");
    std::vector<std::string> labels(ntype);
    for (std::size_t it = 0; it < ntype; ++it) labels[it] = ucell.atoms[it].label;

    this->skfiles_.resize(ntype * ntype);
    const std::string separator = config.skf_directory[config.skf_directory.size() - 1] == '/' ? "" : "/";
    for (std::size_t a = 0; a < ntype; ++a)
    {
        for (std::size_t b = 0; b < ntype; ++b)
        {
            const std::string path = config.skf_directory + separator + labels[a] + "-" + labels[b] + ".skf";
            this->skfiles_[a * ntype + b] = ModuleDFTB::SkfData::read_legacy(path, a == b);
        }
    }
    for (std::size_t a = 0; a < ntype; ++a)
        for (std::size_t b = a + 1; b < ntype; ++b)
            ModuleDFTB::validate_pair_directions(this->skfiles_[a * ntype + b], this->skfiles_[b * ntype + a]);

    this->template_ = ModuleDFTB::DftbPeriodicInput();
    this->template_.atoms.resize(static_cast<std::size_t>(ucell.nat));
    this->template_.pair_parameters.clear();
    for (std::size_t a = 0; a < ntype; ++a)
    {
        for (std::size_t b = a; b < ntype; ++b)
        {
            ModuleDFTB::DftbPairParameters pair;
            pair.species_a = a;
            pair.species_b = b;
            pair.ab = &this->skfiles_[a * ntype + b];
            pair.ba = &this->skfiles_[b * ntype + a];
            this->template_.pair_parameters.push_back(pair);
        }
    }
    this->template_.hubbard_derivative.resize(ntype, 0.0);
    for (std::size_t it = 0; it < ntype; ++it)
    {
        const std::map<std::string, double>::const_iterator derivative = config.hubbard_derivatives.find(labels[it]);
        if (config.third_order && derivative == config.hubbard_derivatives.end())
            throw std::runtime_error("Missing Hubbard derivative for DFTB3 species " + labels[it]);
        if (derivative != config.hubbard_derivatives.end()) this->template_.hubbard_derivative[it] = derivative->second;
    }
    this->template_.kpoints = read_native_kpoints(inp.kpoint_file, ucell);
    if (!config.band_path_file.empty())
        this->template_.band_kpoints = read_native_band_path(config.band_path_file);
    this->template_.thermal_energy_hartree = config.temperature_kelvin * kelvin_to_hartree;
    this->template_.scc_tolerance = config.scc_tolerance;
    this->template_.maximum_scc_iterations = config.maximum_scc_iterations;
    this->template_.mixing_parameter = config.mixing_parameter;
    this->template_.mixing_method = config.mixing_method;
    this->template_.mixing_history = config.mixing_history;
    this->template_.third_order = config.third_order;
    this->output_precision_ = config.output_precision;
    this->template_.total_electrons = 0.0;

    std::size_t atom_index = 0;
    for (std::size_t it = 0; it < ntype; ++it)
    {
        for (int ia = 0; ia < ucell.atoms[it].na; ++ia)
        {
            ModuleDFTB::DftbSpAtom& atom = this->template_.atoms[atom_index++];
            atom.species = it;
            atom.homonuclear_data = &this->skfiles_[it * ntype + it];
            this->template_.total_electrons += atom.homonuclear_data->valence_electron_count();
        }
    }
    if (atom_index != this->template_.atoms.size())
        throw std::runtime_error("ABACUS atom ordering/total does not match native DFTB species data");
}

ModuleDFTB::DftbPeriodicInput ESolver_DFTBNative::make_geometry(const UnitCell& ucell) const
{
    if (ucell.nat != static_cast<int>(this->template_.atoms.size()))
        throw std::runtime_error("Atom count changed after native DFTB initialization");
    ModuleDFTB::DftbPeriodicInput input = this->template_;
    const double length_bohr = ucell.lat0_angstrom * angstrom_to_bohr;
    input.lattice_bohr[0] = {{ucell.latvec.e11 * length_bohr, ucell.latvec.e12 * length_bohr, ucell.latvec.e13 * length_bohr}};
    input.lattice_bohr[1] = {{ucell.latvec.e21 * length_bohr, ucell.latvec.e22 * length_bohr, ucell.latvec.e23 * length_bohr}};
    input.lattice_bohr[2] = {{ucell.latvec.e31 * length_bohr, ucell.latvec.e32 * length_bohr, ucell.latvec.e33 * length_bohr}};
    std::size_t atom_index = 0;
    for (int it = 0; it < ucell.ntype; ++it)
    {
        for (int ia = 0; ia < ucell.atoms[it].na; ++ia)
        {
            const ModuleBase::Vector3<double>& tau = ucell.atoms[it].tau[ia];
            input.atoms[atom_index].position_bohr = {{tau.x * length_bohr, tau.y * length_bohr, tau.z * length_bohr}};
            ++atom_index;
        }
    }
    return input;
}

void ESolver_DFTBNative::runner(BaseCell& cell, const int istep)
{
    static_cast<void>(istep);
    if (cell.kind() != BaseCell::Kind::unit_cell)
        ModuleBase::WARNING_QUIT("ESolver_DFTBNative", "Native DFTB requires a periodic UnitCell.");
    const UnitCell& ucell = static_cast<const UnitCell&>(cell);
    const ModuleDFTB::DftbPeriodicInput input = this->make_geometry(ucell);
    const int rank = Parallel_Common::get_rank();
    std::ostream& running_log = GlobalV::ofs_running;
    std::string output_dir = PARAM.globalv.global_out_dir;
    if (!output_dir.empty() && output_dir.back() != '/') output_dir += '/';
    const std::string dftb_log_file = output_dir + "dftb.log";
    std::ofstream dftb_log;
    if (rank == 0)
    {
        dftb_log.open(dftb_log_file.c_str());
        if (dftb_log)
        {
            dftb_log << "# Native periodic DFTB SCC calculation log\n"
                     << "# Module: ABACUS source_dftb\n"
                     << "# Method: " << (input.third_order ? "DFTB3" : "DFTB2") << " / SCC\n"
                     << "# Settings file: " << this->inp_->dftb_native_input << "\n"
                     << "# Spin-degenerate occupations, total electrons: " << input.total_electrons << "\n"
                     << "# K-points: " << input.kpoints.size() << ", temperature: "
                     << input.thermal_energy_hartree / kelvin_to_hartree << " K\n"
                     << "# SCC threshold: max |delta q| <= " << input.scc_tolerance << " e\n"
                     << "# Mixer: " << input.mixing_method << ", beta=" << input.mixing_parameter
                     << ", history=" << input.mixing_history << "\n"
                     << "# Each row evaluates the DFTB energy functional on the output Mulliken charges of that iteration.\n"
                     << "# Diff_electronic is the change in that electronic energy; SCC_error is max |delta q|.\n"
                     << "# The converged variational free energy and component decomposition follow the iteration table.\n"
                     << "# iSCC E_electronic(Ha) Diff_electronic(Ha) SCC_error(e) net_electron_excess(e)"
                     << " F_band(Ha) E_Fermi(Ha) mixer\n";
            dftb_log << std::scientific << std::setprecision(this->output_precision_);
        }
    }
    int log_open_failure = (rank == 0 && !dftb_log) ? 1 : 0;
    int all_ranks_opened_log = log_open_failure == 0 ? 1 : 0;
    Parallel_Reduce::reduce_min(all_ranks_opened_log);
    log_open_failure = all_ranks_opened_log == 0 ? 1 : 0;
    if (log_open_failure != 0)
        ModuleBase::WARNING_QUIT("ESolver_DFTBNative", "Cannot open native DFTB iteration log: " + dftb_log_file);

    const std::function<void(const ModuleDFTB::DftbSccIteration&)> iteration_logger =
        [this, &dftb_log, &running_log, rank](const ModuleDFTB::DftbSccIteration& iteration) {
            if (rank != 0) return;
            running_log << std::setprecision(std::max(1, this->inp_->out_ndigits))
                                 << " Native DFTB SCC iteration " << std::setw(4) << iteration.iteration
                                 << ": max |delta q|=" << iteration.maximum_charge_residual << " e"
                                 << ", net excess=" << iteration.net_electron_excess << " e"
                                 << ", E_electronic=" << iteration.electronic_energy_hartree << " Ha"
                                 << ", delta E=";
            if (iteration.has_previous_energy)
                running_log << iteration.electronic_energy_change_hartree << " Ha";
            else
                running_log << "n/a";
            running_log << ", F_band=" << iteration.band_free_energy_hartree << " Ha"
                        << ", E_Fermi=" << iteration.fermi_energy_hartree << " Ha"
                        << ", mixer=" << iteration.mixing_step << std::endl;
            running_log.flush();
            dftb_log << std::setw(6) << iteration.iteration << ' ' << iteration.electronic_energy_hartree << ' ';
            if (iteration.has_previous_energy) dftb_log << iteration.electronic_energy_change_hartree;
            else dftb_log << "nan";
            dftb_log << ' ' << iteration.maximum_charge_residual << ' ' << iteration.net_electron_excess << ' '
                     << iteration.band_free_energy_hartree << ' ' << iteration.fermi_energy_hartree << ' '
                     << iteration.mixing_step << std::endl;
            dftb_log.flush();
        };
    this->result_ = ModuleDFTB::solve_periodic_dftb(input, iteration_logger);
    this->conv_esolver = this->result_.converged;
    this->energy_ry_ = 2.0 * this->result_.total_free_energy_hartree;

    if (rank != 0) return;

    running_log << std::setprecision(16)
                         << "\n Native DFTB SCC converged: " << std::boolalpha << this->result_.converged
                         << ", iterations=" << this->result_.scc_iterations
                         << ", max |delta q|=" << this->result_.maximum_charge_residual << " e\n"
                         << " Native DFTB Fermi level: " << this->result_.fermi_energy_hartree << " Ha\n"
                         << " Native DFTB band free energy: " << this->result_.band_free_energy_hartree << " Ha\n"
                         << " Native DFTB H0 energy: " << this->result_.h0_energy_hartree << " Ha\n"
                         << " Native DFTB SCC energy: " << this->result_.scc_energy_hartree << " Ha\n"
                         << " Native DFTB third-order energy: " << this->result_.third_order_energy_hartree << " Ha\n"
                         << " Native DFTB repulsive energy: " << this->result_.repulsive_energy_hartree << " Ha\n"
                         << " Native DFTB total free energy: " << this->result_.total_free_energy_hartree << " Ha\n"
                         << " Native DFTB total free energy: " << this->result_.total_free_energy_hartree * 27.211386245988 << " eV\n"
                         << " Native DFTB component sum: "
                         << this->result_.h0_energy_hartree + this->result_.scc_energy_hartree
                            + this->result_.third_order_energy_hartree + this->result_.repulsive_energy_hartree
                         << " Ha\n";
    dftb_log << "# SCC converged: " << std::boolalpha << this->result_.converged
             << ", iterations: " << this->result_.scc_iterations
             << ", maximum charge residual: " << this->result_.maximum_charge_residual << " e\n"
             << "# Fermi energy: " << this->result_.fermi_energy_hartree << " Ha\n"
             << "# Band energy: " << this->result_.band_energy_hartree << " Ha\n"
             << "# Band free energy: " << this->result_.band_free_energy_hartree << " Ha\n"
             << "# H0 energy: " << this->result_.h0_energy_hartree << " Ha\n"
             << "# SCC energy: " << this->result_.scc_energy_hartree << " Ha\n"
             << "# Third-order energy: " << this->result_.third_order_energy_hartree << " Ha\n"
             << "# Repulsive energy: " << this->result_.repulsive_energy_hartree << " Ha\n"
             << "# Total free energy: " << this->result_.total_free_energy_hartree << " Ha / "
             << this->result_.total_free_energy_hartree * 27.211386245988 << " eV\n";
    double total_electron_excess = 0.0;
    std::size_t atom_index = 0;
    running_log << " Native DFTB Mulliken populations and electron excess charges:\n"
                << "  atom  element  population(e)  electron_excess(e)\n";
    for (int it = 0; it < ucell.ntype; ++it)
    {
        const auto& species = ucell.atoms[it];
        double species_excess = 0.0;
        for (int ia = 0; ia < species.na; ++ia, ++atom_index)
        {
            const double charge = this->result_.electron_excess_charges[atom_index];
            const double population = this->result_.atomic_electron_populations[atom_index];
            species_excess += charge;
            total_electron_excess += charge;
            running_log << "  " << std::setw(4) << atom_index + 1 << "  " << std::setw(7) << species.label
                        << "  " << std::setw(18) << population << "  " << charge << "\n";
        }
        if (species.na > 0)
            running_log << " Native DFTB mean electron excess " << species.label << ": "
                        << species_excess / species.na << " e\n";
    }
    running_log << " Native DFTB net electron excess: " << total_electron_excess << " e\n"
                << " #TOTAL ENERGY# " << this->energy_ry_ * ModuleBase::Ry_to_eV << " eV (native DFTB)\n";
    dftb_log << "# Net electron excess: " << total_electron_excess << " e\n"
             << "# Atomic Mulliken populations and excess charges\n"
             << "# atom element population(e) electron_excess(e)\n";
    atom_index = 0;
    for (int it = 0; it < ucell.ntype; ++it)
    {
        const auto& species = ucell.atoms[it];
        for (int ia = 0; ia < species.na; ++ia, ++atom_index)
            dftb_log << atom_index + 1 << ' ' << species.label << ' '
                     << this->result_.atomic_electron_populations[atom_index] << ' '
                     << this->result_.electron_excess_charges[atom_index] << '\n';
    }
    dftb_log.flush();

    const std::string eigen_file = output_dir + "eig_occ.txt";
    std::ofstream eigenvalues(eigen_file.c_str());
    if (!eigenvalues)
        ModuleBase::WARNING_QUIT("ESolver_DFTBNative", "Cannot write native DFTB eigenvalue file: " + eigen_file);
    eigenvalues << std::setprecision(this->output_precision_)
                << "# K-point eigenvalues and occupations from the converged native DFTB SCC solution\n"
                << "# Module: Native DFTB eigensolver\n"
                << "# Units: energy in eV; occupations are electrons per spin-degenerate state; k-point weights sum to 1\n"
                << "# k_index weight kx(Direct) ky(Direct) kz(Direct) band_index energy(eV) occupation(e)\n";
    for (std::size_t ik = 0; ik < this->result_.kpoint_eigenvalues.size(); ++ik)
    {
        const auto& point = this->result_.kpoint_eigenvalues[ik];
        for (std::size_t band = 0; band < point.eigenvalues_hartree.size(); ++band)
            eigenvalues << ik + 1 << ' ' << point.weight << ' '
                        << point.fractional[0] << ' ' << point.fractional[1] << ' ' << point.fractional[2] << ' '
                        << band + 1 << ' ' << point.eigenvalues_hartree[band] * 27.211386245988 << ' '
                        << point.occupations[band] << '\n';
    }
    eigenvalues.close();
    running_log << " Native DFTB eigenvalues and occupations: wrote " << eigen_file << "\n";

    const std::string mulliken_file = output_dir + "mulliken.txt";
    std::ofstream mulliken(mulliken_file.c_str());
    if (!mulliken)
        ModuleBase::WARNING_QUIT("ESolver_DFTBNative", "Cannot write native DFTB Mulliken file: " + mulliken_file);
    mulliken << std::setprecision(this->output_precision_)
             << "# Atomic Mulliken populations from the converged native DFTB SCC solution\n"
             << "# Module: Native DFTB Mulliken analysis\n"
             << "# Units: population and electron excess in e; excess = population - SKF neutral valence\n"
             << "# atom element population(e) electron_excess(e)\n";
    atom_index = 0;
    for (int it = 0; it < ucell.ntype; ++it)
    {
        const auto& species = ucell.atoms[it];
        for (int ia = 0; ia < species.na; ++ia, ++atom_index)
            mulliken << atom_index + 1 << ' ' << species.label << ' '
                     << this->result_.atomic_electron_populations[atom_index] << ' '
                     << this->result_.electron_excess_charges[atom_index] << '\n';
    }
    mulliken.close();
    running_log << " Native DFTB Mulliken populations: wrote " << mulliken_file << "\n";

    if (!this->result_.band_structure.empty())
    {
        const std::string band_file = output_dir + "band.txt";
        std::ofstream bands(band_file.c_str());
        if (!bands)
            ModuleBase::WARNING_QUIT("ESolver_DFTBNative", "Cannot write native DFTB band file: " + band_file);
        bands << std::setprecision(this->output_precision_)
              << "# Band energies from the converged native DFTB SCC potential (frozen-potential path solve)\n"
              << "# Module: Native DFTB band structure\n"
              << "# Units: k_distance in bohr^-1, eigenvalues and E-Ef in eV; coordinates are fractional reciprocal\n"
              << "# Fermi energy: " << this->result_.fermi_energy_hartree * 27.211386245988 << " eV\n"
              << "# k_index k_distance(bohr^-1) kx ky kz label band_index energy(eV) E-Ef(eV)\n";
        for (std::size_t ik = 0; ik < this->result_.band_structure.size(); ++ik)
        {
            const auto& point = this->result_.band_structure[ik];
            for (std::size_t band = 0; band < point.eigenvalues_hartree.size(); ++band)
            {
                const double energy = point.eigenvalues_hartree[band];
                bands << ik + 1 << ' ' << point.distance_inverse_bohr << ' '
                      << point.fractional[0] << ' ' << point.fractional[1] << ' ' << point.fractional[2] << ' '
                      << (point.label.empty() ? "-" : point.label) << ' ' << band + 1 << ' '
                      << energy * 27.211386245988 << ' '
                      << (energy - this->result_.fermi_energy_hartree) * 27.211386245988 << '\n';
            }
        }
        bands.close();
        running_log << " Native DFTB frozen-potential band structure: "
                    << this->result_.band_structure.size() << " k-points, "
                    << (this->result_.band_structure.front().eigenvalues_hartree.size())
                    << " bands; wrote " << band_file << "\n";
    }
    dftb_log.close();
}

double ESolver_DFTBNative::cal_energy()
{
    return this->energy_ry_;
}

void ESolver_DFTBNative::cal_force(BaseCell& cell, ModuleBase::matrix& force)
{
    static_cast<void>(cell);
    static_cast<void>(force);
    ModuleBase::WARNING_QUIT("ESolver_DFTBNative", "Native DFTB analytic forces are not implemented.");
}

void ESolver_DFTBNative::cal_stress(BaseCell& cell, ModuleBase::matrix& stress)
{
    static_cast<void>(cell);
    static_cast<void>(stress);
    ModuleBase::WARNING_QUIT("ESolver_DFTBNative", "Native DFTB stress derivatives are not implemented.");
}

void ESolver_DFTBNative::after_all_runners(BaseCell& cell)
{
    static_cast<void>(cell);
    GlobalV::ofs_running << std::setprecision(16)
                         << "\n --------------------------------------------" << std::endl
                         << " !FINAL_ETOT_IS " << this->energy_ry_ * ModuleBase::Ry_to_eV
                         << " eV (native DFTB)" << std::endl;
}
} // namespace ModuleESolver
