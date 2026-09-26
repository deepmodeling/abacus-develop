#include "esolver_dftb_native.h"

#include "source_base/constants.h"
#include "source_base/global_variable.h"
#include "source_base/tool_quit.h"
#include "source_cell/unitcell.h"
#include "source_io/module_output/output_log.h"
#include "source_io/module_parameter/input_parameter.h"
#include "source_io/module_parameter/parameter.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#ifdef __MPI
#include <mpi.h>
#endif
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
    bool third_order = true;
};

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
        else if (key == "third_order")
        {
            std::string value;
            if (!(row >> value)) throw std::runtime_error("Invalid third_order value");
            std::transform(value.begin(), value.end(), value.begin(), ::tolower);
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
        || config.maximum_scc_iterations <= 0 || !(config.mixing_parameter > 0.0 && config.mixing_parameter <= 1.0))
        throw std::runtime_error("Invalid native DFTB temperature or SCC controls");
    return config;
}

std::vector<ModuleDFTB::DftbWeightedKPoint> read_native_kpoints(const std::string& filename)
{
    std::ifstream input(filename.c_str());
    if (!input) throw std::runtime_error("Cannot open ABACUS KPT file: " + filename);
    std::string header;
    int count = -1;
    if (!(input >> header >> count) || header != "K_POINTS" || count < 0)
        throw std::runtime_error("Native DFTB expects a valid ABACUS KPT file");
    std::vector<ModuleDFTB::DftbWeightedKPoint> result;
    if (count == 0)
    {
        std::string mode;
        std::array<int, 3> mesh{};
        std::array<double, 3> offset{{0.0, 0.0, 0.0}};
        if (!(input >> mode >> mesh[0] >> mesh[1] >> mesh[2]))
            throw std::runtime_error("Malformed automatic mesh in KPT file");
        input >> offset[0] >> offset[1] >> offset[2];
        if (!input) input.clear();
        if (mode != "Gamma" && mode != "Monkhorst-Pack" && mode != "MP" && mode != "mp")
            throw std::runtime_error("Native DFTB supports Gamma or Monkhorst-Pack KPT meshes");
        if (mesh[0] <= 0 || mesh[1] <= 0 || mesh[2] <= 0)
            throw std::runtime_error("KPT mesh dimensions must be positive");
        const bool monkhorst_pack = mode != "Gamma";
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
        double weight_sum = 0.0;
        for (int i = 0; i < count; ++i)
        {
            ModuleDFTB::DftbWeightedKPoint point;
            if (!(input >> point.fractional[0] >> point.fractional[1] >> point.fractional[2] >> point.weight)
                || !(point.weight > 0.0))
                throw std::runtime_error("Malformed explicit k-point list in KPT file");
            weight_sum += point.weight;
            result.push_back(point);
        }
        if (!(weight_sum > 0.0)) throw std::runtime_error("KPT weights must have a positive sum");
        for (std::size_t i = 0; i < result.size(); ++i) result[i].weight /= weight_sum;
    }
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
    if (cell.kind() != BaseCell::Kind::unitcell)
        ModuleBase::WARNING_QUIT("ESolver_DFTBNative", "Native DFTB requires a periodic UnitCell.");
    const UnitCell& ucell = static_cast<const UnitCell&>(cell);
    if (inp.nspin != 1 || inp.noncolin)
        ModuleBase::WARNING_QUIT("ESolver_DFTBNative", "Native DFTB currently supports spin-degenerate calculations only.");

    int rank = 0;
#ifdef __MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
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
#ifdef __MPI
    int global_status = 0;
    MPI_Allreduce(&local_status, &global_status, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (global_status != 0)
    {
        int error_rank = local_status != 0 ? rank : std::numeric_limits<int>::max();
        int first_error_rank = 0;
        MPI_Allreduce(&error_rank, &first_error_rank, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
        MPI_Bcast(error_message, static_cast<int>(sizeof(error_message)), MPI_CHAR,
                  first_error_rank, MPI_COMM_WORLD);
        ModuleBase::WARNING_QUIT("ESolver_DFTBNative", error_message);
    }
#else
    if (local_status != 0) ModuleBase::WARNING_QUIT("ESolver_DFTBNative", error_message);
#endif
    if (rank == 0)
    {
        GlobalV::ofs_running << " Native periodic DFTB model loaded from " << inp.dftb_native_input
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
    this->template_.kpoints = read_native_kpoints(inp.kpoint_file);
    if (!config.band_path_file.empty())
        this->template_.band_kpoints = read_native_band_path(config.band_path_file);
    this->template_.thermal_energy_hartree = config.temperature_kelvin * kelvin_to_hartree;
    this->template_.scc_tolerance = config.scc_tolerance;
    this->template_.maximum_scc_iterations = config.maximum_scc_iterations;
    this->template_.mixing_parameter = config.mixing_parameter;
    this->template_.third_order = config.third_order;
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
    if (cell.kind() != BaseCell::Kind::unitcell)
        ModuleBase::WARNING_QUIT("ESolver_DFTBNative", "Native DFTB requires a periodic UnitCell.");
    const UnitCell& ucell = static_cast<const UnitCell&>(cell);
    const ModuleDFTB::DftbPeriodicInput input = this->make_geometry(ucell);
    this->result_ = ModuleDFTB::solve_periodic_dftb(input);
    this->conv_esolver = this->result_.converged;
    this->energy_ry_ = 2.0 * this->result_.total_free_energy_hartree;

    int mpi_rank = 0;
#ifdef __MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
#endif
    if (mpi_rank != 0) return;

    GlobalV::ofs_running << std::setprecision(16)
                         << "\n Native DFTB SCC iteration history (iteration, max |delta q| in e):\n";
    for (std::size_t i = 0; i < this->result_.scc_residual_history.size(); ++i)
        GlobalV::ofs_running << "  " << std::setw(4) << i + 1 << "  "
                             << this->result_.scc_residual_history[i] << "\n";
    GlobalV::ofs_running << " Native DFTB SCC converged: " << std::boolalpha << this->result_.converged
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
    double total_electron_excess = 0.0;
    std::size_t atom_index = 0;
    GlobalV::ofs_running << " Native DFTB Mulliken populations and electron excess charges:\n"
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
            GlobalV::ofs_running << "  " << std::setw(4) << atom_index + 1 << "  " << std::setw(7) << species.label
                                 << "  " << std::setw(18) << population << "  " << charge << "\n";
        }
        if (species.na > 0)
            GlobalV::ofs_running << " Native DFTB mean electron excess " << species.label << ": "
                                 << species_excess / species.na << " e\n";
    }
    GlobalV::ofs_running << " Native DFTB net electron excess: " << total_electron_excess << " e\n"
                         << " #TOTAL ENERGY# " << this->energy_ry_ * ModuleBase::Ry_to_eV << " eV (native DFTB)\n";

    if (!this->result_.band_structure.empty())
    {
        std::string output_dir = PARAM.globalv.global_out_dir;
        if (!output_dir.empty() && output_dir.back() != '/') output_dir += '/';
        const std::string band_file = output_dir + "band_structure.dat";
        std::ofstream bands(band_file.c_str());
        if (!bands)
            ModuleBase::WARNING_QUIT("ESolver_DFTBNative", "Cannot write native DFTB band file: " + band_file);
        bands << std::setprecision(16)
              << "# Native DFTB bands from the converged SCC potential (non-self-consistent path solve)\n"
              << "# Fermi energy: " << this->result_.fermi_energy_hartree << " Ha\n"
              << "# Columns: k_index k_distance(bohr^-1) kx ky kz label band_index energy(Ha) energy(eV) E-Ef(eV)\n";
        for (std::size_t ik = 0; ik < this->result_.band_structure.size(); ++ik)
        {
            const auto& point = this->result_.band_structure[ik];
            for (std::size_t band = 0; band < point.eigenvalues_hartree.size(); ++band)
            {
                const double energy = point.eigenvalues_hartree[band];
                bands << ik + 1 << ' ' << point.distance_inverse_bohr << ' '
                      << point.fractional[0] << ' ' << point.fractional[1] << ' ' << point.fractional[2] << ' '
                      << (point.label.empty() ? "-" : point.label) << ' ' << band + 1 << ' '
                      << energy << ' ' << energy * 27.211386245988 << ' '
                      << (energy - this->result_.fermi_energy_hartree) * 27.211386245988 << '\n';
            }
        }
        bands.close();
        GlobalV::ofs_running << " Native DFTB frozen-potential band structure: "
                             << this->result_.band_structure.size() << " k-points, "
                             << (this->result_.band_structure.front().eigenvalues_hartree.size())
                             << " bands; wrote " << band_file << "\n";
    }
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
