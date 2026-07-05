#include "source_lcao/module_ri/sternheimer_abacus_st_smoke.h"

#include "source_basis/module_ao/ORB_read.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_cell/unitcell.h"
#include "source_estate/elecstate.h"
#include "source_estate/module_pot/potential_new.h"
#include "source_hamilt/module_xc/exx_info.h"
#include "source_io/module_parameter/parameter.h"
#include "source_lcao/module_ri/conv_coulomb_pot_k.h"
#include "source_lcao/module_ri/exx_abfs-construct_orbs.h"
#include "source_lcao/module_ri/sternheimer_abfs_perturbation.h"
#include "source_lcao/module_ri/sternheimer_delta.h"
#include "source_lcao/module_ri/sternheimer_fd_solver.h"
#include "source_lcao/module_ri/sternheimer_rpa.h"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <map>
#include <stdexcept>
#include <sstream>
#include <string>
#include <utility>

namespace ModuleRI
{
namespace
{

constexpr const char* kSmokeEnv = "ABACUS_STERNHEIMER_FD_ST_SMOKE";
constexpr const char* kOutputEnv = "ABACUS_STERNHEIMER_FD_ST_OUT";
constexpr const char* kBandsEnv = "ABACUS_STERNHEIMER_FD_ST_BANDS";
constexpr const char* kChannelsEnv = "ABACUS_STERNHEIMER_FD_ST_CHANNELS";
constexpr const char* kMaxDenseEnv = "ABACUS_STERNHEIMER_FD_ST_MAX_DENSE";
constexpr const char* kLanczosSubspaceEnv = "ABACUS_STERNHEIMER_FD_ST_LANCZOS_SUBSPACE";
constexpr const char* kOmegaEnv = "ABACUS_STERNHEIMER_FD_ST_OMEGA";
constexpr const char* kSolverToleranceEnv = "ABACUS_STERNHEIMER_FD_ST_SOLVER_TOL";
constexpr const char* kSolverMaxIterEnv = "ABACUS_STERNHEIMER_FD_ST_MAX_ITER";
constexpr const char* kPCAThresholdEnv = "ABACUS_STERNHEIMER_FD_ST_PCA_THRESHOLD";
constexpr const char* kCCPRmeshTimesEnv = "ABACUS_STERNHEIMER_FD_ST_CCP_RMESH_TIMES";
constexpr const char* kOrbitalDirEnv = "ABACUS_STERNHEIMER_FD_ST_ORBITAL_DIR";
constexpr const char* kOrbitalFilesEnv = "ABACUS_STERNHEIMER_FD_ST_ORBITAL_FILES";

std::string lower_string(std::string value)
{
    std::transform(value.begin(), value.end(), value.begin(), [](const unsigned char ch) {
        return static_cast<char>(std::tolower(ch));
    });
    return value;
}

bool env_is_true(const char* name)
{
    const char* raw = std::getenv(name);
    if (raw == nullptr)
    {
        return false;
    }
    const std::string value = lower_string(raw);
    return !(value.empty() || value == "0" || value == "false" || value == "off" || value == "no");
}

int positive_int_from_env(const char* name, const int default_value)
{
    const char* raw = std::getenv(name);
    if (raw == nullptr)
    {
        return default_value;
    }
    std::size_t parsed = 0;
    const int value = std::stoi(raw, &parsed);
    if (raw[parsed] != '\0' || value <= 0)
    {
        throw std::invalid_argument(std::string("Invalid positive integer in ") + name + ".");
    }
    return value;
}

double positive_double_from_env(const char* name, const double default_value)
{
    const char* raw = std::getenv(name);
    if (raw == nullptr)
    {
        return default_value;
    }
    std::size_t parsed = 0;
    const double value = std::stod(raw, &parsed);
    if (raw[parsed] != '\0' || value <= 0.0)
    {
        throw std::invalid_argument(std::string("Invalid positive floating-point value in ") + name + ".");
    }
    return value;
}

double nonnegative_double_from_env(const char* name, const double default_value)
{
    const char* raw = std::getenv(name);
    if (raw == nullptr)
    {
        return default_value;
    }
    std::size_t parsed = 0;
    const double value = std::stod(raw, &parsed);
    if (raw[parsed] != '\0' || value < 0.0)
    {
        throw std::invalid_argument(std::string("Invalid non-negative floating-point value in ") + name + ".");
    }
    return value;
}

std::string default_report_path(const std::string& output_dir)
{
    const char* explicit_path = std::getenv(kOutputEnv);
    if (explicit_path != nullptr && explicit_path[0] != '\0')
    {
        return explicit_path;
    }
    if (output_dir.empty())
    {
        return "STERNHEIMER_FD_ST.dat";
    }
    if (output_dir.back() == '/')
    {
        return output_dir + "STERNHEIMER_FD_ST.dat";
    }
    return output_dir + "/STERNHEIMER_FD_ST.dat";
}

std::string join_output_path(const std::string& output_dir, const std::string& filename)
{
    if (output_dir.empty())
    {
        return filename;
    }
    if (output_dir.back() == '/')
    {
        return output_dir + filename;
    }
    return output_dir + "/" + filename;
}

std::string chi0_status_path(const std::string& output_dir)
{
    return join_output_path(output_dir, "STERNHEIMER_CHI0.dat");
}

std::string chi0_v1_filename(const int iq, const int ifrequency)
{
    std::ostringstream out;
    out << "v1_sternheimer_chi0_iq_" << iq << "_ifreq_" << ifrequency << "_rank" << GlobalV::MY_RANK << ".dat";
    return out.str();
}

int occupied_band_count(const elecstate::ElecState& elec_state, const int k_index)
{
    int count = 0;
    for (int ib = 0; ib != elec_state.wg.nc; ++ib)
    {
        if (elec_state.wg(k_index, ib) > 1.0e-8)
        {
            count = ib + 1;
        }
    }
    return count;
}

std::vector<double> eigenvalues_ry_from_elec_state(const elecstate::ElecState& elec_state, const int k_index)
{
    std::vector<double> eigenvalues;
    eigenvalues.reserve(static_cast<std::size_t>(elec_state.ekb.nc));
    for (int ib = 0; ib != elec_state.ekb.nc; ++ib)
    {
        eigenvalues.push_back(elec_state.ekb(k_index, ib));
    }
    return eigenvalues;
}

std::vector<double> occupations_from_elec_state(const elecstate::ElecState& elec_state, const int k_index)
{
    std::vector<double> occupations;
    occupations.reserve(static_cast<std::size_t>(elec_state.wg.nc));
    for (int ib = 0; ib != elec_state.wg.nc; ++ib)
    {
        occupations.push_back(elec_state.wg(k_index, ib));
    }
    return occupations;
}

std::vector<std::string> split_orbital_file_list(const std::string& raw)
{
    std::vector<std::string> files;
    std::string current;
    for (const char ch: raw)
    {
        if (ch == ',' || ch == ':')
        {
            if (!current.empty())
            {
                files.push_back(current);
                current.clear();
            }
        }
        else
        {
            current.push_back(ch);
        }
    }
    if (!current.empty())
    {
        files.push_back(current);
    }
    return files;
}

std::string join_path(const std::string& dir, const std::string& file)
{
    if (file.empty())
    {
        return file;
    }
    if (!file.empty() && file.front() == '/')
    {
        return file;
    }
    if (dir.empty() || dir.back() == '/')
    {
        return dir + file;
    }
    return dir + "/" + file;
}

void validate_orbital_files(const std::string& orbital_dir, const std::vector<std::string>& orbital_files)
{
    if (orbital_files.empty())
    {
        throw std::runtime_error(
            "No numerical orbital files are available. Set NUMERICAL_ORBITAL in STRU or "
            "ABACUS_STERNHEIMER_FD_ST_ORBITAL_FILES.");
    }
    for (const std::string& file: orbital_files)
    {
        const std::string path = join_path(orbital_dir, file);
        std::ifstream in(path.c_str());
        if (!in)
        {
            throw std::runtime_error("Cannot open Sternheimer ST orbital file: " + path + ".");
        }
    }
}

std::vector<std::string> orbital_files_from_env_or_cell(const UnitCell& ucell)
{
    const char* raw = std::getenv(kOrbitalFilesEnv);
    if (raw != nullptr && raw[0] != '\0')
    {
        return split_orbital_file_list(raw);
    }
    return ucell.orbital_fn;
}

std::string orbital_dir_from_env_or_input()
{
    const char* raw = std::getenv(kOrbitalDirEnv);
    if (raw != nullptr)
    {
        return raw;
    }
    return PARAM.inp.orbital_dir;
}

void read_sternheimer_orbitals(const UnitCell& ucell, LCAO_Orbitals& orb)
{
    std::vector<std::string> orbital_files = orbital_files_from_env_or_cell(ucell);
    if (static_cast<int>(orbital_files.size()) != ucell.ntype)
    {
        throw std::runtime_error("Sternheimer ST requires one numerical orbital file per atom type.");
    }

    std::string orbital_dir = orbital_dir_from_env_or_input();
    const bool all_absolute = std::all_of(orbital_files.begin(), orbital_files.end(), [](const std::string& file) {
        return !file.empty() && file.front() == '/';
    });
    if (all_absolute)
    {
        orbital_dir.clear();
    }
    validate_orbital_files(orbital_dir, orbital_files);

    const double lcao_ecut = PARAM.inp.lcao_ecut > 0.0 ? PARAM.inp.lcao_ecut : PARAM.inp.ecutwfc;
    if (lcao_ecut <= 0.0)
    {
        throw std::runtime_error("Sternheimer ST requires positive lcao_ecut or ecutwfc to read orbitals.");
    }

    orb.init(GlobalV::ofs_running,
             ucell.ntype,
             orbital_dir,
             orbital_files.data(),
             ucell.descriptor_file,
             PARAM.inp.lmaxmax,
             lcao_ecut,
             PARAM.inp.lcao_dk,
             PARAM.inp.lcao_dr,
             PARAM.inp.lcao_rmax,
             PARAM.globalv.deepks_setorb,
             PARAM.inp.out_mat_r,
             PARAM.inp.cal_force,
             GlobalV::MY_RANK);
}

std::map<Conv_Coulomb_Pot_K::Coulomb_Type, std::vector<std::map<std::string, std::string>>>
make_fock_hartree_coulomb_param()
{
    return {{Conv_Coulomb_Pot_K::Coulomb_Type::Fock, {{{"alpha", "1"}, {"singularity_correction", "limits"}}}}};
}

std::vector<SternheimerABFGridChannel> build_abfs_ccp_grid_channels(const UnitCell& ucell,
                                                                    const SternheimerFDHamiltonian::Grid& grid,
                                                                    const int max_channels,
                                                                    const double pca_threshold,
                                                                    const double ccp_rmesh_times)
{
    LCAO_Orbitals orb;
    read_sternheimer_orbitals(ucell, orb);
    auto lcaos = Exx_Abfs::Construct_Orbs::change_orbs(orb, GlobalC::exx_info.info_ri.kmesh_times);
    Exx_Abfs::Construct_Orbs::filter_empty_orbs(lcaos);

    auto abfs = Exx_Abfs::Construct_Orbs::abfs_same_atom(ucell,
                                                         orb,
                                                         lcaos,
                                                         GlobalC::exx_info.info_ri.kmesh_times,
                                                         pca_threshold);
    Exx_Abfs::Construct_Orbs::filter_empty_orbs(abfs);
    const auto abfs_ccp = Conv_Coulomb_Pot_K::cal_orbs_ccp(abfs, make_fock_hartree_coulomb_param(), ccp_rmesh_times);
    const auto radials_by_type = make_sternheimer_radial_perturbations_from_orbitals(abfs_ccp);

    std::vector<int> atom_types;
    std::vector<ModuleBase::Vector3<double>> atom_positions;
    atom_types.reserve(ucell.nat);
    atom_positions.reserve(ucell.nat);
    for (int it = 0; it != ucell.ntype; ++it)
    {
        const Atom& atom = ucell.atoms[it];
        for (int ia = 0; ia != atom.na; ++ia)
        {
            atom_types.push_back(it);
            atom_positions.push_back(atom.tau[ia] * ucell.lat0);
        }
    }

    return sample_sternheimer_abf_grid_channels(radials_by_type, atom_types, atom_positions, grid, max_channels);
}

std::vector<SternheimerFDHamiltonian::Vector> build_lcao_candidate_grid_orbitals(
    const UnitCell& ucell,
    const SternheimerFDHamiltonian::Grid& grid)
{
    LCAO_Orbitals orb;
    read_sternheimer_orbitals(ucell, orb);
    auto lcaos = Exx_Abfs::Construct_Orbs::change_orbs(orb, GlobalC::exx_info.info_ri.kmesh_times);
    Exx_Abfs::Construct_Orbs::filter_empty_orbs(lcaos);
    const auto radials_by_type = make_sternheimer_radial_perturbations_from_orbitals(lcaos);

    std::vector<int> atom_types;
    std::vector<ModuleBase::Vector3<double>> atom_positions;
    atom_types.reserve(ucell.nat);
    atom_positions.reserve(ucell.nat);
    for (int it = 0; it != ucell.ntype; ++it)
    {
        const Atom& atom = ucell.atoms[it];
        for (int ia = 0; ia != atom.na; ++ia)
        {
            atom_types.push_back(it);
            atom_positions.push_back(atom.tau[ia] * ucell.lat0);
        }
    }

    const std::vector<SternheimerABFGridChannel> channels
        = sample_sternheimer_abf_grid_channels(radials_by_type, atom_types, atom_positions, grid, -1);
    std::vector<SternheimerFDHamiltonian::Vector> candidates;
    candidates.reserve(channels.size());
    for (const SternheimerABFGridChannel& channel: channels)
    {
        SternheimerFDHamiltonian::Vector candidate(
            channel.potential_r.size(),
            SternheimerFDHamiltonian::Complex(0.0, 0.0));
        for (std::size_t ir = 0; ir != channel.potential_r.size(); ++ir)
        {
            candidate[ir] = SternheimerFDHamiltonian::Complex(channel.potential_r[ir], 0.0);
        }
        candidates.push_back(std::move(candidate));
    }
    return candidates;
}

std::vector<std::vector<double>> collect_channel_potentials(const std::vector<SternheimerABFGridChannel>& channels)
{
    std::vector<std::vector<double>> potentials;
    potentials.reserve(channels.size());
    for (const SternheimerABFGridChannel& channel: channels)
    {
        potentials.push_back(channel.potential_r);
    }
    return potentials;
}

SternheimerRPA::Chi0V1Metadata make_chi0_v1_metadata(const UnitCell& ucell,
                                                     const std::vector<SternheimerABFGridChannel>& channels,
                                                     const int ifrequency,
                                                     const double omega_ha,
                                                     const double weight_ha)
{
    SternheimerRPA::Chi0V1Metadata metadata;
    metadata.iq = 1;
    metadata.ifrequency = ifrequency;
    metadata.omega = omega_ha;
    metadata.weight = weight_ha;
    metadata.atom_naux.assign(static_cast<std::size_t>(ucell.nat), 0);
    for (const SternheimerABFGridChannel& channel: channels)
    {
        if (channel.atom_index < 0 || channel.atom_index >= ucell.nat)
        {
            throw std::runtime_error("Sternheimer chi0 found an ABFS channel with invalid atom index.");
        }
        if (channel.atom_local_index < 0)
        {
            throw std::runtime_error("Sternheimer chi0 found an ABFS channel with invalid local index.");
        }
        int& atom_naux = metadata.atom_naux[static_cast<std::size_t>(channel.atom_index)];
        atom_naux = std::max(atom_naux, channel.atom_local_index + 1);
    }
    return metadata;
}

std::vector<SternheimerRPA::AuxiliaryChannel> make_chi0_auxiliary_channels(
    const std::vector<SternheimerABFGridChannel>& channels)
{
    std::vector<SternheimerRPA::AuxiliaryChannel> auxiliary_channels;
    auxiliary_channels.reserve(channels.size());
    for (const SternheimerABFGridChannel& channel: channels)
    {
        SternheimerRPA::AuxiliaryChannel auxiliary_channel;
        auxiliary_channel.channel_index = channel.channel_index;
        auxiliary_channel.atom_index = channel.atom_index;
        auxiliary_channel.atom_local_index = channel.atom_local_index;
        auxiliary_channels.push_back(auxiliary_channel);
    }
    return auxiliary_channels;
}

void write_chi0_index_file(const std::string& filename,
                           const std::vector<std::pair<std::string, SternheimerRPA::Chi0V1Metadata>>& entries)
{
    std::ofstream out(filename.c_str(), std::ios::out | std::ios::trunc);
    if (!out)
    {
        throw std::runtime_error("Failed to open Sternheimer chi0 index file: " + filename);
    }
    out << std::setprecision(16);
    out << "# ABACUS Sternheimer chi0 v1 index\n";
    out << "format_version 1\n";
    out << "marker " << SternheimerRPA::chi0_v1_marker() << '\n';
    out << "nq 1\n";
    out << "nfreq " << entries.size() << '\n';
    out << "iq ifrequency omega_Ha weight_Ha filename\n";
    for (const auto& entry: entries)
    {
        const SternheimerRPA::Chi0V1Metadata& metadata = entry.second;
        out << metadata.iq << ' ' << metadata.ifrequency << ' ' << metadata.omega << ' ' << metadata.weight << ' '
            << entry.first << '\n';
    }
}

SternheimerFDZeroOrderStates solve_fd_zero_order_auto(const SternheimerFDHamiltonian& hamiltonian,
                                                      const int num_bands,
                                                      const double volume_element,
                                                      const int max_dense_size,
                                                      const int lanczos_max_subspace_size)
{
    if (hamiltonian.grid().size() <= max_dense_size)
    {
        return solve_sternheimer_fd_zero_order_dense(hamiltonian, num_bands, volume_element, max_dense_size);
    }

    SternheimerFDLanczosOptions options;
    options.max_subspace_size = lanczos_max_subspace_size;
    options.residual_tolerance = 1.0e-8;
    options.initial_seed = 1;
    return solve_sternheimer_fd_zero_order_lanczos(hamiltonian, num_bands, volume_element, options);
}

std::vector<SternheimerFDHamiltonian::Vector> occupied_wavefunctions_from_states(
    const SternheimerFDZeroOrderStates& states,
    const elecstate::ElecState& elec_state,
    const int k_index)
{
    std::vector<SternheimerFDHamiltonian::Vector> occupied;
    for (int ib = 0; ib != static_cast<int>(states.wavefunctions.size()); ++ib)
    {
        if (elec_state.wg(k_index, ib) > 1.0e-8)
        {
            occupied.push_back(states.wavefunctions[ib]);
        }
    }
    return occupied;
}

void write_failure_report(std::ofstream& out, const std::string& reason)
{
    out << "# ABACUS Sternheimer FD linear-response smoke test\n";
    out << "status failed\n";
    out << "reason " << reason << '\n';
}

} // namespace

bool sternheimer_abacus_st_smoke_enabled()
{
    return env_is_true(kSmokeEnv);
}

void run_sternheimer_abacus_st_smoke(const elecstate::Potential& potential,
                                     const ModulePW::PW_Basis& pw_basis,
                                     const UnitCell& ucell,
                                     const elecstate::ElecState& elec_state,
                                     const std::string& output_dir)
{
    if (!sternheimer_abacus_st_smoke_enabled() || GlobalV::MY_RANK != 0)
    {
        return;
    }

    const std::string report_path = default_report_path(output_dir);
    std::ofstream out(report_path.c_str(), std::ios::out | std::ios::trunc);
    if (!out)
    {
        GlobalV::ofs_running << " Sternheimer FD ST smoke: failed to open " << report_path << std::endl;
        return;
    }

    try
    {
        if (GlobalV::NPROC != 1)
        {
            throw std::runtime_error(
                "The current FD ST smoke test requires a single MPI rank so pw_basis.nrxx is the full grid.");
        }
        if (elec_state.ekb.nc <= 0 || elec_state.wg.nc <= 0)
        {
            throw std::runtime_error("ABACUS DFT eigenvalues or occupations are not available.");
        }

        const int requested_bands = positive_int_from_env(kBandsEnv, 1);
        const int num_bands = std::min(requested_bands, elec_state.ekb.nc);
        const int max_channels = positive_int_from_env(kChannelsEnv, 1);
        const int max_dense_size = positive_int_from_env(kMaxDenseEnv, 4096);
        const int lanczos_max_subspace_size = positive_int_from_env(kLanczosSubspaceEnv, 320);
        const double omega = nonnegative_double_from_env(kOmegaEnv, 0.5);
        const double solver_tolerance = positive_double_from_env(kSolverToleranceEnv, 1.0e-8);
        const int solver_max_iter = positive_int_from_env(kSolverMaxIterEnv, 300);
        const double pca_threshold = nonnegative_double_from_env(kPCAThresholdEnv, PARAM.inp.exx_pca_threshold);
        const double ccp_rmesh_times = positive_double_from_env(kCCPRmeshTimesEnv, PARAM.inp.rpa_ccp_rmesh_times);

        SternheimerABACUSSTSmokeResult result;
        result.grid_data = make_sternheimer_fd_grid(pw_basis);
        result.omega = omega;
        result.pca_threshold = pca_threshold;
        result.ccp_rmesh_times = ccp_rmesh_times;
        result.perturbation_source = "abfs_ccp";

        const SternheimerFDHamiltonian hamiltonian = make_sternheimer_fd_hamiltonian(potential, pw_basis, ucell, 0, 1.0);
        const SternheimerFDZeroOrderStates states = solve_fd_zero_order_auto(hamiltonian,
                                                                             num_bands,
                                                                             result.grid_data.volume_element,
                                                                             max_dense_size,
                                                                             lanczos_max_subspace_size);
        const std::vector<SternheimerFDHamiltonian::Vector> occupied
            = occupied_wavefunctions_from_states(states, elec_state, 0);
        if (occupied.empty())
        {
            throw std::runtime_error("No occupied FD zero-order states are available for Sternheimer ST.");
        }

        const std::vector<SternheimerABFGridChannel> channels
            = build_abfs_ccp_grid_channels(ucell, result.grid_data.grid, max_channels, pca_threshold, ccp_rmesh_times);
        result.num_available_channels = static_cast<int>(channels.size());
        if (channels.empty())
        {
            throw std::runtime_error("No ABFS CCP perturbation channels were generated.");
        }

        auto dot = [&result](const SternheimerFDHamiltonian::Vector& lhs,
                             const SternheimerFDHamiltonian::Vector& rhs) {
            return sternheimer_fd_grid_dot(lhs, rhs, result.grid_data.volume_element);
        };

        SternheimerRPA::SolverOptions solver_options;
        solver_options.max_iter = solver_max_iter;
        solver_options.residual_tol = solver_tolerance;

        for (int ib = 0; ib != static_cast<int>(states.wavefunctions.size()); ++ib)
        {
            const double occupation = elec_state.wg(0, ib);
            if (occupation <= 1.0e-8)
            {
                continue;
            }

            for (const SternheimerABFGridChannel& channel: channels)
            {
                SternheimerFDHamiltonian::Vector rhs;
                SternheimerRPA::build_rhs_from_hartree_perturbation(channel.potential_r,
                                                                    states.wavefunctions[ib],
                                                                    rhs);
                SternheimerFDHamiltonian::Vector projected_rhs = rhs;
                SternheimerRPA::project_out_subspace(occupied, dot, projected_rhs);

                const SternheimerFDLinearResponse response
                    = solve_sternheimer_fd_linear_response(hamiltonian,
                                                           occupied,
                                                           states.eigenvalues[ib],
                                                           rhs,
                                                           omega,
                                                           result.grid_data.volume_element,
                                                           solver_options);

                SternheimerABACUSSTChannelResult channel_result;
                channel_result.band_index = ib;
                channel_result.channel_index = channel.channel_index;
                channel_result.atom_index = channel.atom_index;
                channel_result.angular_momentum = channel.angular_momentum;
                channel_result.radial_index = channel.radial_index;
                channel_result.magnetic_index = channel.magnetic_index;
                channel_result.fd_eigenvalue = states.eigenvalues[ib];
                channel_result.occupation = occupation;
                channel_result.rhs_norm = sternheimer_fd_grid_norm(rhs, result.grid_data.volume_element);
                channel_result.projected_rhs_norm
                    = sternheimer_fd_grid_norm(projected_rhs, result.grid_data.volume_element);
                channel_result.solver_converged = response.solver.converged;
                channel_result.solver_iterations = response.solver.iterations;
                channel_result.solver_relative_residual = response.solver.relative_residual;
                channel_result.equation_residual_norm = response.residual_norm;
                channel_result.polarizability = SternheimerRPA::accumulate_polarizability_grid_element(
                    channel.potential_r,
                    states.wavefunctions[ib],
                    response.delta_wavefunction,
                    result.grid_data.volume_element);
                result.channels.push_back(channel_result);
            }
        }

        out << format_sternheimer_abacus_st_report(result);
        GlobalV::ofs_running << " Sternheimer FD ST smoke report: " << report_path << std::endl;
    }
    catch (const std::exception& error)
    {
        write_failure_report(out, error.what());
        GlobalV::ofs_running << " Sternheimer FD ST smoke failed: " << error.what() << std::endl;
        GlobalV::ofs_running << " Sternheimer FD ST smoke report: " << report_path << std::endl;
    }
}

void run_sternheimer_abacus_chi0_output(const elecstate::Potential& potential,
                                        const ModulePW::PW_Basis& pw_basis,
                                        const UnitCell& ucell,
                                        const elecstate::ElecState& elec_state,
                                        const std::string& output_dir)
{
    if (!PARAM.inp.out_sternheimer_librpa || GlobalV::MY_RANK != 0)
    {
        return;
    }

    const std::string status_path = chi0_status_path(output_dir);
    std::ofstream out(status_path.c_str(), std::ios::out | std::ios::trunc);
    if (!out)
    {
        GlobalV::ofs_running << " Sternheimer chi0 output: failed to open " << status_path << std::endl;
        return;
    }
    out << std::setprecision(16);
    out << "# ABACUS Sternheimer chi0 output for LibRPA\n";

    try
    {
        if (PARAM.inp.out_librpa_reader_version != 1)
        {
            throw std::runtime_error("out_sternheimer_librpa requires out_librpa_reader_version=1.");
        }
        if (GlobalV::NPROC != 1)
        {
            throw std::runtime_error(
                "The current Sternheimer chi0 output requires a single MPI rank so pw_basis.nrxx is the full grid.");
        }
        if (elec_state.ekb.nc <= 0 || elec_state.wg.nc <= 0)
        {
            throw std::runtime_error("ABACUS DFT eigenvalues or occupations are not available.");
        }

        const int occupied_count = occupied_band_count(elec_state, 0);
        if (occupied_count <= 0)
        {
            throw std::runtime_error("No occupied DFT bands are available for Sternheimer chi0 output.");
        }

        const int requested_bands = positive_int_from_env(kBandsEnv, occupied_count);
        const int num_bands = std::min(requested_bands, elec_state.ekb.nc);
        if (num_bands < occupied_count)
        {
            throw std::runtime_error("Sternheimer chi0 output requires all occupied bands; use the smoke test for "
                                     "band-limited debugging.");
        }

        const int max_dense_size = positive_int_from_env(kMaxDenseEnv, 4096);
        const int lanczos_max_subspace_size = positive_int_from_env(kLanczosSubspaceEnv, 320);
        const double solver_tolerance = positive_double_from_env(kSolverToleranceEnv, 1.0e-8);
        const int solver_max_iter = positive_int_from_env(kSolverMaxIterEnv, 300);
        const double pca_threshold = nonnegative_double_from_env(kPCAThresholdEnv, PARAM.inp.exx_pca_threshold);
        const double ccp_rmesh_times = positive_double_from_env(kCCPRmeshTimesEnv, PARAM.inp.rpa_ccp_rmesh_times);
        const int nfreq = PARAM.inp.sternheimer_nfreq;
        const std::string frequency_grid_file = PARAM.inp.sternheimer_frequency_grid_file;
        const bool use_delta_sternheimer = PARAM.inp.sternheimer_delta;

        const std::vector<double> eigenvalues_ry = eigenvalues_ry_from_elec_state(elec_state, 0);
        const std::vector<double> occupations = occupations_from_elec_state(elec_state, 0);
        const SternheimerRPA::TransitionEnergyWindow transition_window
            = SternheimerRPA::transition_energy_window_from_eigenvalues_ry(eigenvalues_ry, occupations);
        const bool use_frequency_grid_file = !frequency_grid_file.empty();
        const std::string frequency_grid_source = use_frequency_grid_file ? "file" : "greenx_minimax";
        const SternheimerRPA::FrequencyGrid frequency_grid
            = use_frequency_grid_file
                  ? SternheimerRPA::read_frequency_grid_file(frequency_grid_file, nfreq)
                  : SternheimerRPA::generate_greenx_minimax_frequency_grid(nfreq,
                                                                           transition_window.emin_ha,
                                                                           transition_window.emax_ha);

        const SternheimerABACUSFDGridData grid_data = make_sternheimer_fd_grid(pw_basis);
        const SternheimerFDHamiltonian hamiltonian = make_sternheimer_fd_hamiltonian(potential, pw_basis, ucell, 0, 1.0);
        const SternheimerFDZeroOrderStates states = solve_fd_zero_order_auto(hamiltonian,
                                                                             num_bands,
                                                                             grid_data.volume_element,
                                                                             max_dense_size,
                                                                             lanczos_max_subspace_size);
        const std::vector<SternheimerFDHamiltonian::Vector> occupied
            = occupied_wavefunctions_from_states(states, elec_state, 0);
        if (occupied.empty())
        {
            throw std::runtime_error("No occupied FD zero-order states are available for Sternheimer chi0 output.");
        }

        const std::vector<SternheimerABFGridChannel> channels
            = build_abfs_ccp_grid_channels(ucell, grid_data.grid, -1, pca_threshold, ccp_rmesh_times);
        if (channels.empty())
        {
            throw std::runtime_error("No ABFS CCP perturbation channels were generated.");
        }

        const std::vector<std::vector<double>> potentials = collect_channel_potentials(channels);
        const int num_channels = static_cast<int>(channels.size());

        SternheimerDeltaSubspace delta_subspace;
        if (use_delta_sternheimer)
        {
            const std::vector<SternheimerFDHamiltonian::Vector> candidate_orbitals
                = build_lcao_candidate_grid_orbitals(ucell, grid_data.grid);
            if (candidate_orbitals.empty())
            {
                throw std::runtime_error("Sternheimer delta mode found no sampled LCAO candidate orbitals.");
            }

            SternheimerDeltaSubspaceOptions delta_options;
            delta_options.max_virtual_states = PARAM.inp.sternheimer_delta_max_states;
            delta_options.norm_tolerance = PARAM.inp.sternheimer_delta_norm_tol;
            delta_subspace = build_delta_sternheimer_subspace(hamiltonian,
                                                              occupied,
                                                              candidate_orbitals,
                                                              grid_data.volume_element,
                                                              delta_options);
            if (delta_subspace.virtual_states.empty())
            {
                throw std::runtime_error("Sternheimer delta mode produced no fixed virtual states.");
            }
        }

        SternheimerRPA::SolverOptions solver_options;
        solver_options.max_iter = solver_max_iter;
        solver_options.residual_tol = solver_tolerance;

        bool all_converged = true;
        int solved_equations = 0;
        double max_solver_relative_residual = 0.0;
        double max_equation_residual_norm = 0.0;
        std::vector<std::pair<std::string, SternheimerRPA::Chi0V1Metadata>> index_entries;
        index_entries.reserve(frequency_grid.omega_ha.size());

        for (int ifrequency = 0; ifrequency != nfreq; ++ifrequency)
        {
            std::vector<SternheimerRPA::Complex> chi0_branch(
                static_cast<std::size_t>(num_channels) * static_cast<std::size_t>(num_channels),
                SternheimerRPA::Complex(0.0, 0.0));

            const double omega_ha = frequency_grid.omega_ha[static_cast<std::size_t>(ifrequency)];
            const double omega_ry = 2.0 * omega_ha;

            for (int ib = 0; ib != static_cast<int>(states.wavefunctions.size()); ++ib)
            {
                const double occupation = elec_state.wg(0, ib);
                if (occupation <= 1.0e-8)
                {
                    continue;
                }

                for (int ichannel = 0; ichannel != num_channels; ++ichannel)
                {
                    SternheimerFDHamiltonian::Vector rhs;
                    SternheimerRPA::build_rhs_from_hartree_perturbation(channels[static_cast<std::size_t>(ichannel)]
                                                                            .potential_r,
                                                                        states.wavefunctions[ib],
                                                                        rhs);
                    SternheimerFDHamiltonian::Vector delta_wavefunction;
                    SternheimerRPA::SolverResult solver_result;
                    double equation_residual_norm = 0.0;
                    if (use_delta_sternheimer)
                    {
                        const std::vector<SternheimerFDHamiltonian::Complex> perturbation_matrix_elements
                            = delta_sternheimer_perturbation_matrix_elements(
                                delta_subspace.virtual_states,
                                channels[static_cast<std::size_t>(ichannel)].potential_r,
                                states.wavefunctions[ib],
                                grid_data.volume_element);
                        const SternheimerDeltaLinearResponse response
                            = solve_delta_sternheimer_linear_response(hamiltonian,
                                                                      occupied,
                                                                      states.eigenvalues[ib],
                                                                      rhs,
                                                                      delta_subspace.virtual_states,
                                                                      perturbation_matrix_elements,
                                                                      omega_ry,
                                                                      grid_data.volume_element,
                                                                      solver_options);
                        delta_wavefunction = response.response.reconstructed_wavefunction;
                        solver_result = response.solver;
                        equation_residual_norm = response.residual_norm;
                    }
                    else
                    {
                        const SternheimerFDLinearResponse response
                            = solve_sternheimer_fd_linear_response(hamiltonian,
                                                                   occupied,
                                                                   states.eigenvalues[ib],
                                                                   rhs,
                                                                   omega_ry,
                                                                   grid_data.volume_element,
                                                                   solver_options);
                        delta_wavefunction = response.delta_wavefunction;
                        solver_result = response.solver;
                        equation_residual_norm = response.residual_norm;
                    }
                    SternheimerRPA::accumulate_chi0_branch_column(potentials,
                                                                  states.wavefunctions[ib],
                                                                  delta_wavefunction,
                                                                  grid_data.volume_element,
                                                                  occupation,
                                                                  ichannel,
                                                                  chi0_branch);
                    all_converged = all_converged && solver_result.converged;
                    ++solved_equations;
                    max_solver_relative_residual
                        = std::max(max_solver_relative_residual, solver_result.relative_residual);
                    max_equation_residual_norm = std::max(max_equation_residual_norm, equation_residual_norm);
                }
            }

            const std::vector<SternheimerRPA::Complex> chi0
                = SternheimerRPA::symmetrize_chi0_imaginary_frequency(chi0_branch, num_channels);
            const SternheimerRPA::Chi0V1Metadata metadata
                = make_chi0_v1_metadata(ucell,
                                        channels,
                                        ifrequency + 1,
                                        omega_ha,
                                        frequency_grid.weights_ha[static_cast<std::size_t>(ifrequency)]);
            const std::vector<SternheimerRPA::AuxiliaryChannel> auxiliary_channels
                = make_chi0_auxiliary_channels(channels);
            const std::string data_file = chi0_v1_filename(metadata.iq, metadata.ifrequency);
            SternheimerRPA::write_chi0_v1_file(data_file, metadata, auxiliary_channels, chi0);
            index_entries.push_back({data_file, metadata});
            GlobalV::ofs_running << " Sternheimer chi0 v1 output: " << data_file << std::endl;
        }

        write_chi0_index_file("v1_sternheimer_chi0_index.dat", index_entries);

        const int grid_size = grid_data.grid.nx * grid_data.grid.ny * grid_data.grid.nz;
        out << "status success\n";
        out << "format v1\n";
        out << "data_files " << index_entries.size() << '\n';
        out << "index_file v1_sternheimer_chi0_index.dat\n";
        out << "grid " << grid_data.grid.nx << ' ' << grid_data.grid.ny << ' ' << grid_data.grid.nz << " size "
            << grid_size << " dV " << grid_data.volume_element << '\n';
        out << "nfreq " << nfreq << '\n';
        out << "frequency_grid_source " << frequency_grid_source << '\n';
        if (use_frequency_grid_file)
        {
            out << "frequency_grid_file " << frequency_grid_file << '\n';
        }
        out << "transition_window_Ha " << transition_window.emin_ha << ' ' << transition_window.emax_ha << '\n';
        out << "ifrequency omega_Ha weight_Ha omega_Ry data_file\n";
        for (const auto& entry: index_entries)
        {
            const SternheimerRPA::Chi0V1Metadata& metadata = entry.second;
            out << metadata.ifrequency << ' ' << metadata.omega << ' ' << metadata.weight << ' '
                << 2.0 * metadata.omega << ' ' << entry.first << '\n';
        }
        out << "pca_threshold " << pca_threshold << '\n';
        out << "ccp_rmesh_times " << ccp_rmesh_times << '\n';
        out << "occupied_bands " << occupied.size() << '\n';
        out << "abfs_channels " << num_channels << '\n';
        out << "sternheimer_delta " << (use_delta_sternheimer ? "yes" : "no") << '\n';
        if (use_delta_sternheimer)
        {
            out << "sternheimer_delta_max_states " << PARAM.inp.sternheimer_delta_max_states << '\n';
            out << "sternheimer_delta_norm_tol " << PARAM.inp.sternheimer_delta_norm_tol << '\n';
            out << "sternheimer_delta_virtual_states " << delta_subspace.virtual_states.size() << '\n';
            out << "sternheimer_delta_accepted_candidates " << delta_subspace.accepted_candidates << '\n';
            out << "sternheimer_delta_discarded_candidates " << delta_subspace.discarded_candidates << '\n';
        }
        out << "solved_equations " << solved_equations << '\n';
        out << "all_converged " << (all_converged ? "yes" : "no") << '\n';
        out << "max_solver_relative_residual " << max_solver_relative_residual << '\n';
        out << "max_equation_residual_norm " << max_equation_residual_norm << '\n';
        GlobalV::ofs_running << " Sternheimer chi0 status: " << status_path << std::endl;
    }
    catch (const std::exception& error)
    {
        out << "status failed\n";
        out << "reason " << error.what() << '\n';
        GlobalV::ofs_running << " Sternheimer chi0 output failed: " << error.what() << std::endl;
        GlobalV::ofs_running << " Sternheimer chi0 status: " << status_path << std::endl;
    }
}

} // namespace ModuleRI
