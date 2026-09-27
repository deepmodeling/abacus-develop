#include "esolver_dftb_native.h"
#include "esolver_dftb_native_input.h"

#include "source_base/constants.h"
#include "source_base/parallel_common.h"
#include "source_base/parallel_reduce.h"
#include "source_base/tool_quit.h"
#include "source_cell/unitcell.h"
#include "source_io/module_output/output_log.h"
#include "source_io/module_parameter/input_parameter.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
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
} // namespace

void ESolver_DFTBNative::before_all_runners(BaseCell& cell, const Input_para& inp)
{
    this->inp_ = &inp;
    if (cell.kind() != BaseCell::Kind::unitcell)
        ModuleBase::WARNING_QUIT("ESolver_DFTBNative", "Native DFTB requires a periodic UnitCell.");
    const UnitCell& ucell = static_cast<const UnitCell&>(cell);
    if (inp.basis_type != "dftb")
        ModuleBase::WARNING_QUIT("ESolver_DFTBNative", "Set basis_type=dftb when using esolver_type=dftbnative.");
    if (inp.calculation != "scf")
        ModuleBase::WARNING_QUIT("ESolver_DFTBNative", "Native DFTB currently supports calculation=scf only; use band_path_file for frozen-SCC bands.");
    if (inp.nspin != 1 || inp.noncolin)
        ModuleBase::WARNING_QUIT("ESolver_DFTBNative", "Native DFTB currently supports spin-degenerate calculations only.");

    const int rank = Parallel_Common::get_rank();
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
        std::cout << " Experimental native periodic DFTB model loaded from " << inp.dftb_native_input
                    << " (no DFTB+ runtime, UPF, or ABACUS orbital files)." << std::endl;
    }
}
void ESolver_DFTBNative::load_model(const UnitCell& ucell, const Input_para& inp)
{
    const NativeDftb::NativeDftbConfig config = NativeDftb::read_native_config(inp.dftb_native_input);
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
            ModuleDFTB::validate_pair_directions(this->skfiles_[a * ntype + b],
                                                 this->skfiles_[b * ntype + a],
                                                 1.0e-8);

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
    this->template_.kpoints = NativeDftb::read_native_kpoints(inp.kpoint_file, ucell);
    if (!config.band_path_file.empty())
        this->template_.band_kpoints = NativeDftb::read_native_band_path(config.band_path_file);
    this->template_.thermal_energy_hartree = config.temperature_kelvin * kelvin_to_hartree;
    this->template_.scc_tolerance = config.scc_tolerance;
    this->template_.maximum_scc_iterations = config.maximum_scc_iterations;
    this->template_.mixing_parameter = config.mixing_parameter;
    this->template_.mixing_method = config.mixing_method;
    this->template_.mixing_history = config.mixing_history;
    this->template_.broyden_inverse_jacobi_weight = config.broyden_inverse_jacobi_weight;
    this->template_.broyden_minimal_weight = config.broyden_minimal_weight;
    this->template_.broyden_maximal_weight = config.broyden_maximal_weight;
    this->template_.broyden_weight_factor = config.broyden_weight_factor;
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
    if (cell.kind() != BaseCell::Kind::unitcell)
        ModuleBase::WARNING_QUIT("ESolver_DFTBNative", "Native DFTB requires a periodic UnitCell.");
    const UnitCell& ucell = static_cast<const UnitCell&>(cell);
    const ModuleDFTB::DftbPeriodicInput input = this->make_geometry(ucell);
    const int rank = Parallel_Common::get_rank();
    std::ostream& running_log = std::cout;
    std::string output_dir = this->output_dir_;
    if (!output_dir.empty() && output_dir.back() != '/') output_dir += '/';
    const std::string dftb_log_file = output_dir + "dftb.log";
    std::ofstream dftb_log;
    if (rank == 0)
    {
        dftb_log.open(dftb_log_file.c_str());
        if (dftb_log)
        {
            dftb_log << std::scientific << std::setprecision(this->output_precision_);
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
                     << "# Broyden weights: inverse-Jacobi=" << input.broyden_inverse_jacobi_weight
                     << ", minimum=" << input.broyden_minimal_weight
                     << ", maximum=" << input.broyden_maximal_weight
                     << ", factor=" << input.broyden_weight_factor << "\n"
                     << "# Each row evaluates the DFTB energy functional on the output Mulliken charges of that iteration.\n"
                     << "# Diff_electronic is the change in that electronic energy; SCC_error is max |delta q|.\n"
                     << "# The converged variational free energy and component decomposition follow the iteration table.\n"
                     << "# MPI ranks: " << Parallel_Common::get_size() << "\n"
                     << "# Slater-Koster files (ordered species pair, A-B, B-A):\n";
            for (const auto& pair : input.pair_parameters)
            {
                dftb_log << "#   " << ucell.atoms[pair.species_a].label << '-'
                         << ucell.atoms[pair.species_b].label << "  " << pair.ab->filename
                         << "  " << pair.ba->filename << "\n";
            }
            dftb_log << "# Species data (label, neutral valence, Hubbard derivative Ha/e, onsite s/p Ha, Hubbard U s/p Ha):\n";
            for (int it = 0; it < ucell.ntype; ++it)
            {
                const ModuleDFTB::SkfData& data = this->skfiles_[static_cast<std::size_t>(it) * ucell.ntype + it];
                dftb_log << "#   " << ucell.atoms[it].label << ' ' << data.valence_electron_count() << ' '
                         << input.hubbard_derivative[static_cast<std::size_t>(it)] << ' '
                         << data.onsite_hartree[0] << ' ' << data.onsite_hartree[1] << ' '
                         << data.hubbard_u_hartree[0] << ' ' << data.hubbard_u_hartree[1] << "\n";
            }
            dftb_log << "# Lattice vectors (Bohr; row vectors):\n";
            for (const auto& vector : input.lattice_bohr)
                dftb_log << "#   " << vector[0] << ' ' << vector[1] << ' ' << vector[2] << "\n";
            dftb_log << "# Atom positions (index, element, x, y, z in Bohr):\n";
            for (std::size_t atom = 0; atom < input.atoms.size(); ++atom)
            {
                const auto& position = input.atoms[atom].position_bohr;
                dftb_log << "#   " << atom + 1 << ' ' << ucell.atoms[input.atoms[atom].species].label << ' '
                         << position[0] << ' ' << position[1] << ' ' << position[2] << "\n";
            }
            dftb_log << "# SCC integration k-points (index, direct coordinates, normalized weight):\n";
            for (std::size_t ik = 0; ik < input.kpoints.size(); ++ik)
            {
                const auto& point = input.kpoints[ik];
                dftb_log << "#   " << ik + 1 << ' ' << point.fractional[0] << ' ' << point.fractional[1] << ' '
                         << point.fractional[2] << ' ' << point.weight << "\n";
            }
            dftb_log << "# Frozen-SCC band path points: " << input.band_kpoints.size() << "\n"
                     << "# iSCC E_electronic(Ha) Diff_electronic(Ha) SCC_error(e) net_electron_excess(e)"
                     << " F_band(Ha) E_Fermi(Ha) mixer\n";
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

    this->write_result_outputs(ucell, input, output_dir, running_log, dftb_log);
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

} // namespace ModuleESolver
