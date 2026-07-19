#include "source_lcao/module_ri/sternheimer_abacus_st_smoke.h"
#include "source_lcao/module_ri/singular_value.h"

#include "source_basis/module_ao/ORB_read.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_cell/unitcell.h"
#include "source_estate/elecstate.h"
#include "source_estate/module_pot/potential_new.h"
#include "source_hamilt/module_xc/exx_info.h"
#include "source_io/module_parameter/parameter.h"
#include "source_lcao/module_gint/gint_atom.h"
#include "source_lcao/module_ri/conv_coulomb_pot_k.h"
#include "source_lcao/module_ri/exx_abfs-construct_orbs.h"
#include "source_lcao/module_ri/sternheimer_abfs_perturbation.h"
#include "source_lcao/module_ri/sternheimer_channel_parallel.h"
#include "source_lcao/module_ri/sternheimer_delta.h"
#include "source_lcao/module_ri/sternheimer_abacus_fd_nonlocal.h"
#include "source_lcao/module_ri/sternheimer_fd_solver.h"
#include "source_lcao/module_ri/sternheimer_periodic_solver.h"
#include "source_lcao/module_ri/sternheimer_rpa.h"

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <stdexcept>
#include <sstream>
#include <string>
#include <utility>

#ifdef __MPI
#include <mpi.h>
#endif

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
constexpr const char* kFrequencyRankShiftEnv = "ABACUS_STERNHEIMER_FD_ST_FREQ_RANK_SHIFT";
constexpr const char* kDeltaComponentDiagnosticEnv = "ABACUS_STERNHEIMER_DELTA_COMPONENT_DIAG";
constexpr const char* kLCAOSOSDiagnosticEnv = "ABACUS_STERNHEIMER_LCAO_SOS_DIAG";
constexpr const char* kDeltaABlockModeEnv = "ABACUS_STERNHEIMER_DELTA_A_BLOCK";
constexpr double kHartreeToRydberg = 2.0;

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

SternheimerDeltaABlockMode delta_a_block_mode_from_env()
{
    const char* raw = std::getenv(kDeltaABlockModeEnv);
    if (raw == nullptr || raw[0] == '\0')
    {
        return SternheimerDeltaABlockMode::ReferenceValueGradient;
    }
    return parse_sternheimer_delta_a_block_mode(lower_string(raw));
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

int int_from_env(const char* name, const int default_value)
{
    const char* raw = std::getenv(name);
    if (raw == nullptr)
    {
        return default_value;
    }
    std::size_t parsed = 0;
    const int value = std::stoi(raw, &parsed);
    if (raw[parsed] != '\0')
    {
        throw std::invalid_argument(std::string("Invalid integer in ") + name + ".");
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

std::string chi0_v1_filename(const int iq, const int ifrequency, const int rank = GlobalV::MY_RANK)
{
    std::ostringstream out;
    out << "v1_sternheimer_chi0_iq_" << iq << "_ifreq_" << ifrequency << "_rank" << rank << ".dat";
    return out.str();
}

std::string delta_component_v1_filename(const std::string& component,
                                        const int iq,
                                        const int ifrequency,
                                        const int rank = GlobalV::MY_RANK)
{
    std::ostringstream out;
    out << "v1_sternheimer_delta_" << component << "_iq_" << iq << "_ifreq_" << ifrequency
        << "_rank" << rank << ".dat";
    return out.str();
}

std::string lcao_sos_v1_filename(const int iq, const int ifrequency, const int rank = GlobalV::MY_RANK)
{
    std::ostringstream out;
    out << "v1_sternheimer_lcao_sos_iq_" << iq << "_ifreq_" << ifrequency << "_rank" << rank << ".dat";
    return out.str();
}

std::string chi0_progress_filename(const int rank = GlobalV::MY_RANK)
{
    std::ostringstream out;
    out << "STERNHEIMER_CHI0_PROGRESS_rank" << rank << ".dat";
    return out.str();
}

std::string chi0_failure_filename(const int rank = GlobalV::MY_RANK)
{
    std::ostringstream out;
    out << "STERNHEIMER_CHI0_FAILURE_rank" << rank << ".dat";
    return out.str();
}

double elapsed_seconds_since(const std::chrono::steady_clock::time_point& start)
{
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}

void reset_chi0_progress_file()
{
    std::ofstream out(chi0_progress_filename().c_str());
    out << std::setprecision(16);
    out << "event rank ifrequency owner_rank band channel solved_equations converged iterations "
           "solver_relative_residual equation_residual_norm elapsed_s note\n";
}

void append_chi0_progress_event(const std::string& event,
                                const int ifrequency,
                                const int owner_rank,
                                const int band,
                                const int channel,
                                const int solved_equations,
                                const SternheimerRPA::SolverResult* solver_result,
                                const double equation_residual_norm,
                                const double elapsed_seconds,
                                const std::string& note)
{
    std::ofstream out(chi0_progress_filename().c_str(), std::ios::app);
    out << std::setprecision(16);
    out << event << ' ' << GlobalV::MY_RANK << ' ' << ifrequency << ' ' << owner_rank << ' ' << band << ' '
        << channel << ' ' << solved_equations << ' ';
    if (solver_result == nullptr)
    {
        out << "- -1 -1 ";
    }
    else
    {
        out << (solver_result->converged ? "yes" : "no") << ' ' << solver_result->iterations << ' '
            << solver_result->relative_residual << ' ';
    }
    out << equation_residual_norm << ' ' << elapsed_seconds << ' ' << note << '\n';
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

SternheimerRPA::TransitionEnergyWindow transition_window_from_all_kpoints(
    const elecstate::ElecState& elec_state,
    const std::vector<SternheimerLCAOOccupiedKPoint>& records)
{
    SternheimerRPA::TransitionEnergyWindow combined;
    combined.emin_ha = std::numeric_limits<double>::max();
    combined.emax_ha = 0.0;
    for (const SternheimerLCAOOccupiedKPoint& record: records)
    {
        const SternheimerRPA::TransitionEnergyWindow window
            = SternheimerRPA::transition_energy_window_from_eigenvalues_ry(
                eigenvalues_ry_from_elec_state(elec_state, record.local_k_index),
                occupations_from_elec_state(elec_state, record.local_k_index));
        combined.emin_ha = std::min(combined.emin_ha, window.emin_ha);
        combined.emax_ha = std::max(combined.emax_ha, window.emax_ha);
    }
    return combined;
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

using SternheimerOrbitalSet = std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>;

struct SternheimerABFSInput
{
    std::vector<std::vector<SternheimerRadialPerturbation>> radials_by_type;
    std::vector<int> atom_types;
    std::vector<ModuleBase::Vector3<double>> atom_positions;
};

SternheimerOrbitalSet build_sternheimer_abfs(const UnitCell& ucell, const double pca_threshold)
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
    return abfs;
}

SternheimerABFSInput make_sternheimer_abfs_input(const UnitCell& ucell,
                                                 const SternheimerOrbitalSet& orbitals)
{
    SternheimerABFSInput input;
    input.radials_by_type = make_sternheimer_radial_perturbations_from_orbitals(orbitals);
    input.atom_types.reserve(ucell.nat);
    input.atom_positions.reserve(ucell.nat);
    for (int it = 0; it != ucell.ntype; ++it)
    {
        const Atom& atom = ucell.atoms[it];
        for (int ia = 0; ia != atom.na; ++ia)
        {
            input.atom_types.push_back(it);
            input.atom_positions.push_back(atom.tau[ia] * ucell.lat0);
        }
    }
    return input;
}

SternheimerABFSInput build_abfs_ccp_input(const UnitCell& ucell,
                                          const double pca_threshold,
                                          const double ccp_rmesh_times)
{
    const SternheimerOrbitalSet abfs = build_sternheimer_abfs(ucell, pca_threshold);
    const SternheimerOrbitalSet abfs_ccp
        = Conv_Coulomb_Pot_K::cal_orbs_ccp(abfs, make_fock_hartree_coulomb_param(), ccp_rmesh_times);
    return make_sternheimer_abfs_input(ucell, abfs_ccp);
}

SternheimerABFSInput build_abfs_density_input(const UnitCell& ucell, const double pca_threshold)
{
    return make_sternheimer_abfs_input(ucell, build_sternheimer_abfs(ucell, pca_threshold));
}

std::vector<SternheimerABFGridChannel> build_abfs_ccp_grid_channels(const UnitCell& ucell,
                                                                    const SternheimerFDHamiltonian::Grid& grid,
                                                                    const int max_channels,
                                                                    const double pca_threshold,
                                                                    const double ccp_rmesh_times)
{
    const SternheimerABFSInput input = build_abfs_ccp_input(ucell, pca_threshold, ccp_rmesh_times);
    return sample_sternheimer_abf_grid_channels(
        input.radials_by_type, input.atom_types, input.atom_positions, grid, max_channels);
}

struct SternheimerPeriodicABFGridData
{
    std::vector<SternheimerABFBlochGridChannel> densities;
    std::vector<SternheimerABFBlochGridChannel> potentials;
};

SternheimerPeriodicABFGridData build_abfs_full_coulomb_bloch_grid_channels(
    const UnitCell& ucell,
    const SternheimerFDHamiltonian::Grid& grid,
    const SternheimerReducedKPoint& qpoint,
    const double gamma_inverse_k2,
    const int max_channels_per_atom,
    const double pca_threshold)
{
    const SternheimerABFSInput input = build_abfs_density_input(ucell, pca_threshold);
    const std::vector<SternheimerABFBlochGridChannel> all_densities
        = sample_sternheimer_abf_bloch_grid_channels(
            input.radials_by_type, input.atom_types, input.atom_positions, grid, qpoint, -1);
    SternheimerPeriodicABFGridData result;
    result.densities = limit_sternheimer_abf_channels_per_atom(all_densities, max_channels_per_atom);
    result.potentials
        = solve_sternheimer_abf_periodic_full_coulomb(
            result.densities, grid, qpoint, gamma_inverse_k2);
    return result;
}

std::vector<std::string> find_coulomb_v1_rank_files(const int iq, const int mpi_ranks)
{
    std::vector<std::string> filenames;
    for (int rank = 0; rank != mpi_ranks; ++rank)
    {
        const std::string filename = "v1_coulomb_full_iq_" + std::to_string(iq)
                                     + "_rank" + std::to_string(rank) + ".dat";
        std::ifstream input(filename.c_str(), std::ios::binary);
        if (input.good())
        {
            filenames.push_back(filename);
        }
    }
    if (filenames.empty())
    {
        throw std::runtime_error("Periodic Gamma Sternheimer requires the full-Coulomb reader-v1 rank files.");
    }
    return filenames;
}

std::vector<int> atom_auxiliary_sizes(const std::vector<SternheimerABFBlochGridChannel>& channels,
                                      const int natom)
{
    std::vector<int> sizes(static_cast<std::size_t>(natom), 0);
    for (const SternheimerABFBlochGridChannel& channel: channels)
    {
        if (channel.atom_index < 0 || channel.atom_index >= natom || channel.atom_local_index < 0)
        {
            throw std::runtime_error("Periodic Gamma Sternheimer found invalid auxiliary channel metadata.");
        }
        int& size = sizes[static_cast<std::size_t>(channel.atom_index)];
        size = std::max(size, channel.atom_local_index + 1);
    }
    return sizes;
}

void write_sternheimer_grid_coulomb_diagnostic(
    const std::string& filename,
    const std::vector<SternheimerABFBlochGridChannel>& densities,
    const std::vector<SternheimerABFBlochGridChannel>& potentials,
    const double volume_element)
{
    if (densities.size() != potentials.size())
    {
        throw std::invalid_argument("Sternheimer grid Coulomb density/potential channel counts differ.");
    }
    std::ofstream out(filename.c_str());
    if (!out)
    {
        throw std::runtime_error("Failed to open Sternheimer grid Coulomb diagnostic: " + filename);
    }
    out << std::setprecision(16);
    out << "# row col real imag\n";
    out << "naux " << densities.size() << '\n';
    for (std::size_t row = 0; row != densities.size(); ++row)
    {
        for (std::size_t col = 0; col != potentials.size(); ++col)
        {
            if (densities[row].potential_r.size() != potentials[col].potential_r.size())
            {
                throw std::invalid_argument("Sternheimer grid Coulomb channel sizes differ.");
            }
            SternheimerRPA::Complex value(0.0, 0.0);
            for (std::size_t ir = 0; ir != densities[row].potential_r.size(); ++ir)
            {
                value += std::conj(densities[row].potential_r[ir]) * potentials[col].potential_r[ir];
            }
            value *= volume_element;
            out << row << ' ' << col << ' ' << value.real() << ' ' << value.imag() << '\n';
        }
    }
}

std::vector<SternheimerDeltaGridFunction> build_lcao_candidate_grid_functions(
    const UnitCell& ucell,
    const SternheimerFDHamiltonian::Grid& grid,
    const LCAO_Orbitals* provided_orbitals = nullptr,
    const SternheimerReducedKPoint kpoint = {0.0, 0.0, 0.0})
{
    LCAO_Orbitals loaded_orbitals;
    if (provided_orbitals == nullptr)
    {
        read_sternheimer_orbitals(ucell, loaded_orbitals);
        provided_orbitals = &loaded_orbitals;
    }
    const LCAO_Orbitals& orb = *provided_orbitals;
    const int grid_size = grid.size();
    constexpr int sample_chunk_size = 8192;

    std::vector<SternheimerDeltaGridFunction> candidates;
    int atom_index = 0;
    for (int type = 0; type != ucell.ntype; ++type)
    {
        const Atom& unitcell_atom = ucell.atoms[type];
        const Numerical_Orbital& atom_orbitals = orb.Phi[type];
        Atom sampling_atom;
        sampling_atom.nwl = atom_orbitals.getLmax();
        sampling_atom.l_nchi.resize(static_cast<std::size_t>(sampling_atom.nwl + 1), 0);
        for (int angular_momentum = 0; angular_momentum <= sampling_atom.nwl; ++angular_momentum)
        {
            sampling_atom.l_nchi[static_cast<std::size_t>(angular_momentum)]
                = atom_orbitals.getNchi(angular_momentum);
            sampling_atom.nw += (2 * angular_momentum + 1)
                                * sampling_atom.l_nchi[static_cast<std::size_t>(angular_momentum)];
        }
        sampling_atom.set_index();

        for (int ia = 0; ia != unitcell_atom.na; ++ia, ++atom_index)
        {
            const int orbital_count = sampling_atom.nw;
            const std::size_t candidate_begin = candidates.size();
            candidates.resize(candidate_begin + static_cast<std::size_t>(orbital_count));
            for (int iw = 0; iw != orbital_count; ++iw)
            {
                SternheimerDeltaGridFunction& candidate
                    = candidates[candidate_begin + static_cast<std::size_t>(iw)];
                candidate.values.assign(static_cast<std::size_t>(grid_size),
                                        SternheimerFDHamiltonian::Complex(0.0, 0.0));
                for (SternheimerFDHamiltonian::Vector& gradient: candidate.gradients)
                {
                    gradient.assign(static_cast<std::size_t>(grid_size),
                                    SternheimerFDHamiltonian::Complex(0.0, 0.0));
                }
            }

            ModuleGint::GintAtom sampler(&sampling_atom,
                                         type,
                                         ia,
                                         atom_index,
                                         ModuleGint::Vec3i(0, 0, 0),
                                         ModuleGint::Vec3i(0, 0, 0),
                                         ModuleGint::Vec3d(0.0, 0.0, 0.0),
                                         &orb.Phi[type],
                                         &ucell);
            const ModuleBase::Vector3<double> atom_position = unitcell_atom.tau[ia] * ucell.lat0;
            const std::vector<std::array<int, 3>> periodic_images
                = enumerate_delta_sternheimer_periodic_images(grid,
                                                               {atom_position.x,
                                                                atom_position.y,
                                                                atom_position.z},
                                                               atom_orbitals.getRcut());
            for (int first = 0; first < grid_size; first += sample_chunk_size)
            {
                const int chunk_size = std::min(sample_chunk_size, grid_size - first);
                std::vector<ModuleGint::Vec3d> coordinates(static_cast<std::size_t>(chunk_size));
                const std::size_t buffer_size
                    = static_cast<std::size_t>(chunk_size) * static_cast<std::size_t>(orbital_count);
                std::vector<double> values(buffer_size, 0.0);
                std::array<std::vector<double>, 3> gradients{{
                    std::vector<double>(buffer_size, 0.0),
                    std::vector<double>(buffer_size, 0.0),
                    std::vector<double>(buffer_size, 0.0),
                }};
                for (const std::array<int, 3>& image: periodic_images)
                {
                    for (int local = 0; local != chunk_size; ++local)
                    {
                        const int linear = first + local;
                        const int ix = linear / (grid.ny * grid.nz);
                        const int remainder = linear % (grid.ny * grid.nz);
                        const int iy = remainder / grid.nz;
                        const int iz = remainder % grid.nz;
                        const std::array<double, 3> position
                            = sternheimer_fd_grid_cartesian_position(grid, ix, iy, iz);
                        const std::array<double, 3> translation
                            = sternheimer_fd_grid_lattice_translation(grid, image);
                        coordinates[static_cast<std::size_t>(local)]
                            = ModuleGint::Vec3d(position[0] - atom_position.x - translation[0],
                                               position[1] - atom_position.y - translation[1],
                                               position[2] - atom_position.z - translation[2]);
                    }
                    sampler.set_phi_dphi(coordinates,
                                         orbital_count,
                                         values.data(),
                                         gradients[0].data(),
                                         gradients[1].data(),
                                         gradients[2].data());
                    accumulate_delta_sternheimer_bloch_samples(values,
                                                                gradients,
                                                                chunk_size,
                                                                orbital_count,
                                                                static_cast<std::size_t>(first),
                                                                candidate_begin,
                                                                kpoint,
                                                                image,
                                                                candidates);
                }
            }
        }
    }
    return candidates;
}

struct SternheimerSampledLCAOKPoint
{
    SternheimerFDZeroOrderStates states;
    SternheimerFDZeroOrderStates unoccupied_states;
    std::vector<SternheimerDeltaGridFunction> sampled_ao_functions;
    std::vector<SternheimerDeltaGridFunction> occupied_functions;
    std::vector<SternheimerDeltaGridFunction> occupied_projector_functions;
    double occupied_raw_norm_min = std::numeric_limits<double>::infinity();
    double occupied_raw_norm_max = 0.0;
    double unoccupied_raw_norm_min = std::numeric_limits<double>::infinity();
    double unoccupied_raw_norm_max = 0.0;
};

SternheimerSampledLCAOKPoint sample_sternheimer_lcao_kpoint(
    const UnitCell& ucell,
    const SternheimerFDHamiltonian::Grid& grid,
    const LCAO_Orbitals& orbitals,
    const SternheimerLCAOOccupiedKPoint& record,
    const double volume_element,
    const double norm_tolerance)
{
    SternheimerSampledLCAOKPoint sampled;
    sampled.sampled_ao_functions = build_lcao_candidate_grid_functions(ucell, grid, &orbitals, record.kpoint);
    if (sampled.sampled_ao_functions.empty())
    {
        throw std::runtime_error("Sternheimer periodic LCAO sampling found no AO functions.");
    }

    const std::size_t occupied_count = record.coefficients.size();
    sampled.states.eigenvalues.reserve(occupied_count);
    sampled.states.wavefunctions.reserve(occupied_count);
    sampled.states.residual_norms.reserve(occupied_count);
    sampled.occupied_functions.reserve(occupied_count);
    for (std::size_t ib = 0; ib != occupied_count; ++ib)
    {
        if (record.coefficients[ib].size() != sampled.sampled_ao_functions.size())
        {
            throw std::runtime_error(
                "Sternheimer periodic LCAO coefficient basis size does not match sampled AO functions.");
        }
        SternheimerDeltaGridFunction occupied_function
            = linear_combination_delta_sternheimer_grid_functions(sampled.sampled_ao_functions,
                                                                   record.coefficients[ib]);
        const double norm = sternheimer_fd_grid_norm(occupied_function.values, volume_element);
        if (norm <= norm_tolerance)
        {
            throw std::runtime_error("Sternheimer periodic sampled occupied function has zero norm.");
        }
        sampled.occupied_raw_norm_min = std::min(sampled.occupied_raw_norm_min, norm);
        sampled.occupied_raw_norm_max = std::max(sampled.occupied_raw_norm_max, norm);
        const SternheimerFDHamiltonian::Complex inverse_norm(1.0 / norm, 0.0);
        for (auto& value: occupied_function.values)
        {
            value *= inverse_norm;
        }
        for (auto& gradient: occupied_function.gradients)
        {
            for (auto& value: gradient)
            {
                value *= inverse_norm;
            }
        }
        sampled.states.eigenvalues.push_back(record.eigenvalues[ib]);
        sampled.states.wavefunctions.push_back(occupied_function.values);
        sampled.states.residual_norms.push_back(0.0);
        sampled.occupied_functions.push_back(std::move(occupied_function));
    }
    sampled.unoccupied_states.eigenvalues.reserve(record.unoccupied_eigenvalues.size());
    sampled.unoccupied_states.wavefunctions.reserve(record.unoccupied_coefficients.size());
    sampled.unoccupied_states.residual_norms.reserve(record.unoccupied_coefficients.size());
    for (std::size_t ib = 0; ib != record.unoccupied_coefficients.size(); ++ib)
    {
        if (record.unoccupied_coefficients[ib].size() != sampled.sampled_ao_functions.size())
        {
            throw std::runtime_error(
                "Sternheimer periodic LCAO unoccupied coefficient basis size does not match sampled AO functions.");
        }
        SternheimerDeltaGridFunction unoccupied_function
            = linear_combination_delta_sternheimer_grid_functions(sampled.sampled_ao_functions,
                                                                   record.unoccupied_coefficients[ib]);
        const double norm = sternheimer_fd_grid_norm(unoccupied_function.values, volume_element);
        if (norm <= norm_tolerance)
        {
            throw std::runtime_error("Sternheimer periodic sampled unoccupied function has zero norm.");
        }
        sampled.unoccupied_raw_norm_min = std::min(sampled.unoccupied_raw_norm_min, norm);
        sampled.unoccupied_raw_norm_max = std::max(sampled.unoccupied_raw_norm_max, norm);
        const SternheimerFDHamiltonian::Complex inverse_norm(1.0 / norm, 0.0);
        for (auto& value: unoccupied_function.values)
        {
            value *= inverse_norm;
        }
        sampled.unoccupied_states.eigenvalues.push_back(record.unoccupied_eigenvalues[ib]);
        sampled.unoccupied_states.wavefunctions.push_back(std::move(unoccupied_function.values));
        sampled.unoccupied_states.residual_norms.push_back(0.0);
    }
    sampled.occupied_projector_functions = orthonormalize_delta_sternheimer_grid_functions(
        sampled.occupied_functions, volume_element, norm_tolerance);
    if (sampled.occupied_projector_functions.size() != occupied_count)
    {
        throw std::runtime_error("Sternheimer periodic occupied projector lost a linearly dependent KS state.");
    }
    return sampled;
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

std::vector<SternheimerFDHamiltonian::Vector> collect_channel_potentials(
    const std::vector<SternheimerABFBlochGridChannel>& channels)
{
    std::vector<SternheimerFDHamiltonian::Vector> potentials;
    potentials.reserve(channels.size());
    for (const SternheimerABFBlochGridChannel& channel: channels)
    {
        potentials.push_back(channel.potential_r);
    }
    return potentials;
}

std::vector<double> scale_potential(const std::vector<double>& potential, const double factor)
{
    std::vector<double> scaled = potential;
    for (double& value: scaled)
    {
        value *= factor;
    }
    return scaled;
}

template <typename Scalar>
std::vector<std::vector<Scalar>> scale_potentials(const std::vector<std::vector<Scalar>>& potentials,
                                                  const double factor)
{
    std::vector<std::vector<Scalar>> scaled = potentials;
    for (std::vector<Scalar>& potential: scaled)
    {
        for (Scalar& value: potential)
        {
            value *= factor;
        }
    }
    return scaled;
}

template <typename Channel>
void write_abfs_channel_diagnostic(const std::string& filename, const std::vector<Channel>& channels)
{
    std::ofstream out(filename.c_str(), std::ios::out | std::ios::trunc);
    if (!out)
    {
        throw std::runtime_error("Failed to open Sternheimer ABFS channel diagnostic file: " + filename);
    }
    out << std::setprecision(16);
    out << "# ABACUS Sternheimer ABFS channel diagnostic\n";
    out << "# channel atom atom_local type l radial m label max_abs\n";
    for (const Channel& channel: channels)
    {
        out << channel.channel_index << ' '
            << channel.atom_index << ' '
            << channel.atom_local_index << ' '
            << channel.type_index << ' '
            << channel.angular_momentum << ' '
            << channel.radial_index << ' '
            << channel.magnetic_index << ' '
            << channel.label << ' '
            << channel.max_abs << '\n';
    }
}

template <typename Channel>
SternheimerRPA::Chi0V1Metadata make_chi0_v1_metadata(const UnitCell& ucell,
                                                     const std::vector<Channel>& channels,
                                                     const int iq,
                                                     const int ifrequency,
                                                     const double omega_ha,
                                                     const double weight_ha)
{
    SternheimerRPA::Chi0V1Metadata metadata;
    metadata.iq = iq;
    metadata.ifrequency = ifrequency;
    metadata.omega = omega_ha;
    metadata.weight = weight_ha;
    metadata.atom_naux.assign(static_cast<std::size_t>(ucell.nat), 0);
    for (const Channel& channel: channels)
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

template <typename Channel>
std::vector<SternheimerRPA::AuxiliaryChannel> make_chi0_auxiliary_channels(const std::vector<Channel>& channels)
{
    std::vector<SternheimerRPA::AuxiliaryChannel> auxiliary_channels;
    auxiliary_channels.reserve(channels.size());
    for (const Channel& channel: channels)
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

void broadcast_zero_order_states(SternheimerFDZeroOrderStates& states, const int grid_size, const bool enabled)
{
#ifdef __MPI
    if (!enabled || GlobalV::NPROC <= 1)
    {
        return;
    }

    int num_states = GlobalV::MY_RANK == 0 ? static_cast<int>(states.wavefunctions.size()) : 0;
    MPI_Bcast(&num_states, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (num_states < 0)
    {
        throw std::runtime_error("Invalid Sternheimer zero-order state count broadcast.");
    }

    if (GlobalV::MY_RANK != 0)
    {
        states.eigenvalues.assign(static_cast<std::size_t>(num_states), 0.0);
        states.residual_norms.assign(static_cast<std::size_t>(num_states), 0.0);
        states.wavefunctions.assign(static_cast<std::size_t>(num_states),
                                    SternheimerFDHamiltonian::Vector(
                                        static_cast<std::size_t>(grid_size),
                                        SternheimerFDHamiltonian::Complex(0.0, 0.0)));
    }

    if (num_states > 0)
    {
        MPI_Bcast(states.eigenvalues.data(), num_states, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(states.residual_norms.data(), num_states, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    for (int istate = 0; istate != num_states; ++istate)
    {
        if (static_cast<int>(states.wavefunctions[static_cast<std::size_t>(istate)].size()) != grid_size)
        {
            throw std::runtime_error("Sternheimer zero-order wavefunction size does not match the full grid.");
        }
        MPI_Bcast(reinterpret_cast<double*>(states.wavefunctions[static_cast<std::size_t>(istate)].data()),
                  2 * grid_size,
                  MPI_DOUBLE,
                  0,
                  MPI_COMM_WORLD);
    }
#else
    (void)states;
    (void)grid_size;
    (void)enabled;
#endif
}

void reduce_chi0_output_stats(bool& all_converged,
                              int& solved_equations,
                              double& max_solver_relative_residual,
                              double& max_equation_residual_norm,
                              const bool enabled)
{
#ifdef __MPI
    if (!enabled || GlobalV::NPROC <= 1)
    {
        return;
    }

    int converged_flag = all_converged ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &converged_flag, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    all_converged = converged_flag != 0;

    MPI_Allreduce(MPI_IN_PLACE, &solved_equations, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    double residuals[2] = {max_solver_relative_residual, max_equation_residual_norm};
    MPI_Allreduce(MPI_IN_PLACE, residuals, 2, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    max_solver_relative_residual = residuals[0];
    max_equation_residual_norm = residuals[1];
#else
    (void)all_converged;
    (void)solved_equations;
    (void)max_solver_relative_residual;
    (void)max_equation_residual_norm;
    (void)enabled;
#endif
}

void reduce_kpoint_parallel_complex_vector(std::vector<SternheimerRPA::Complex>& values,
                                           const bool enabled)
{
#ifdef __MPI
    if (!enabled || GlobalV::NPROC <= 1 || values.empty())
    {
        return;
    }
    MPI_Allreduce(MPI_IN_PLACE,
                  values.data(),
                  static_cast<int>(values.size()),
                  MPI_DOUBLE_COMPLEX,
                  MPI_SUM,
                  MPI_COMM_WORLD);
#else
    (void)values;
    (void)enabled;
#endif
}

void reduce_kpoint_parallel_int_vector(std::vector<int>& values, const bool enabled)
{
#ifdef __MPI
    if (!enabled || GlobalV::NPROC <= 1 || values.empty())
    {
        return;
    }
    MPI_Allreduce(MPI_IN_PLACE,
                  values.data(),
                  static_cast<int>(values.size()),
                  MPI_INT,
                  MPI_SUM,
                  MPI_COMM_WORLD);
#else
    (void)values;
    (void)enabled;
#endif
}

void reduce_kpoint_parallel_double_vector(std::vector<double>& values, const bool enabled)
{
#ifdef __MPI
    if (!enabled || GlobalV::NPROC <= 1 || values.empty())
    {
        return;
    }
    MPI_Allreduce(MPI_IN_PLACE,
                  values.data(),
                  static_cast<int>(values.size()),
                  MPI_DOUBLE,
                  MPI_SUM,
                  MPI_COMM_WORLD);
#else
    (void)values;
    (void)enabled;
#endif
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

void run_sternheimer_periodic_lcao_chi0_output(
    const elecstate::Potential& potential,
    const ModulePW::PW_Basis& pw_basis,
    const UnitCell& ucell,
    const elecstate::ElecState& elec_state,
    const LCAO_Orbitals& orbitals,
    const std::vector<SternheimerLCAOOccupiedKPoint>& occupied_kpoints,
    const std::array<int, 3>& kmesh,
    std::ofstream& out,
    const bool use_frequency_mpi,
    const std::chrono::steady_clock::time_point& chi0_start_time)
{
    if (PARAM.inp.nspin != 1)
    {
        throw std::runtime_error("The first periodic Sternheimer driver supports only nspin=1 insulators.");
    }
    if (PARAM.inp.symmetry != "-1")
    {
        throw std::runtime_error("Periodic Sternheimer output requires symmetry=-1 and the full k mesh.");
    }

    const SternheimerPeriodicResponsePlan response_plan
        = build_sternheimer_periodic_response_plan(occupied_kpoints, PARAM.inp.sternheimer_q_index);
    if (PARAM.inp.sternheimer_q_index <= 0)
    {
        throw std::runtime_error("Internal error: the periodic Sternheimer path requires a positive q index.");
    }
    validate_sternheimer_periodic_kmesh(kmesh, static_cast<int>(occupied_kpoints.size()));
    constexpr double q_tolerance = 1.0e-10;
    const bool gamma_qpoint
        = std::all_of(response_plan.qpoint.begin(), response_plan.qpoint.end(), [q_tolerance](const double coordinate) {
              return std::abs(coordinate) <= q_tolerance;
          });
    const double massidda_chi = gamma_qpoint
                                    ? Singular_Value::cal_massidda(ucell, kmesh, 2, 1.0, 5, 1.0e-4)
                                    : 0.0;
    const double gamma_inverse_k2
        = sternheimer_periodic_gamma_inverse_k2(response_plan.qpoint,
                                                 PARAM.inp.exx_singularity_correction,
                                                 massidda_chi);
    append_chi0_progress_event("periodic_plan_ready",
                               0,
                               -1,
                               -1,
                               -1,
                               0,
                               nullptr,
                               -1.0,
                               elapsed_seconds_since(chi0_start_time),
                               "kq_pairs=" + std::to_string(response_plan.kq_pairs.size()));

    const int kpoint_groups = PARAM.globalv.kpar_lcao;
    const bool use_kpoint_mpi = kpoint_groups > 1;
    const bool use_parallel_grid_mpi = use_frequency_mpi || use_kpoint_mpi;
    if (use_frequency_mpi && use_kpoint_mpi)
    {
        throw std::runtime_error(
            "Periodic Sternheimer does not yet support nested frequency MPI and k-point MPI.");
    }

    const double solver_tolerance = positive_double_from_env(kSolverToleranceEnv, 1.0e-8);
    const int solver_max_iter = positive_int_from_env(kSolverMaxIterEnv, 300);
    const double pca_threshold = nonnegative_double_from_env(kPCAThresholdEnv, PARAM.inp.exx_pca_threshold);
    const double ccp_rmesh_times = positive_double_from_env(kCCPRmeshTimesEnv, PARAM.inp.rpa_ccp_rmesh_times);
    const int max_channels = positive_int_from_env(kChannelsEnv, -1);
    const int nfreq = PARAM.inp.sternheimer_nfreq;
    const int default_frequency_rank_shift = use_frequency_mpi && GlobalV::NPROC > 1 ? 1 : 0;
    const int frequency_rank_shift = int_from_env(kFrequencyRankShiftEnv, default_frequency_rank_shift);
    const bool use_delta_sternheimer = PARAM.inp.sternheimer_delta;
    const SternheimerDeltaABlockMode delta_a_block_mode = delta_a_block_mode_from_env();
    const bool write_delta_components = use_delta_sternheimer && env_is_true(kDeltaComponentDiagnosticEnv);
    const bool write_lcao_sos = env_is_true(kLCAOSOSDiagnosticEnv);
    if (write_lcao_sos)
    {
        for (const SternheimerLCAOOccupiedKPoint& record: occupied_kpoints)
        {
            if (record.unoccupied_coefficients.empty())
            {
                throw std::runtime_error(
                    "Periodic direct LCAO-SOS diagnostic requires gathered unoccupied LCAO states.");
            }
        }
    }
    const std::string frequency_grid_file = PARAM.inp.sternheimer_frequency_grid_file;
    const SternheimerRPA::TransitionEnergyWindow transition_window
        = transition_window_from_all_kpoints(elec_state, occupied_kpoints);
    const bool use_frequency_grid_file = !frequency_grid_file.empty();
    const std::string frequency_grid_source = use_frequency_grid_file ? "file" : "greenx_minimax";
    const SternheimerRPA::FrequencyGrid frequency_grid
        = use_frequency_grid_file
              ? SternheimerRPA::read_frequency_grid_file(frequency_grid_file, nfreq)
              : SternheimerRPA::generate_greenx_minimax_frequency_grid(nfreq,
                                                                       transition_window.emin_ha,
                                                                       transition_window.emax_ha);
    append_chi0_progress_event("frequency_grid_ready",
                               0,
                               -1,
                               -1,
                               -1,
                               0,
                               nullptr,
                               -1.0,
                               elapsed_seconds_since(chi0_start_time),
                               "nfreq=" + std::to_string(nfreq));

    const SternheimerABACUSFDGridData grid_data
        = use_parallel_grid_mpi ? make_sternheimer_fd_full_grid(pw_basis) : make_sternheimer_fd_grid(pw_basis);
    const std::vector<double> kpoint_parallel_full_potential
        = use_kpoint_mpi ? copy_sternheimer_full_local_potential(potential, pw_basis, 0)
                         : std::vector<double>();
    append_chi0_progress_event("full_grid_ready",
                               0,
                               -1,
                               -1,
                               -1,
                               0,
                               nullptr,
                               -1.0,
                               elapsed_seconds_since(chi0_start_time),
                               "grid_size=" + std::to_string(grid_data.grid.size()));
    SternheimerPeriodicABFGridData periodic_abfs
        = build_abfs_full_coulomb_bloch_grid_channels(
            ucell,
            grid_data.grid,
            response_plan.qpoint,
            gamma_inverse_k2,
            max_channels,
            pca_threshold);
    if (periodic_abfs.potentials.empty())
    {
        throw std::runtime_error("No periodic ABFS full-Coulomb perturbation channels were generated.");
    }
    const int num_channels = static_cast<int>(periodic_abfs.potentials.size());
    double gamma_projection_relative_error = 0.0;
    if (gamma_qpoint)
    {
        const auto target = SternheimerRPA::read_coulomb_v1_files(
            find_coulomb_v1_rank_files(response_plan.iq, GlobalV::NPROC));
        if (target.iq != response_plan.iq
            || target.atom_naux != atom_auxiliary_sizes(periodic_abfs.densities, ucell.nat))
        {
            throw std::runtime_error("Periodic Gamma Sternheimer full-Coulomb v1 metadata do not match the grid ABFS.");
        }
        gamma_projection_relative_error
            = compare_sternheimer_periodic_coulomb_projection(periodic_abfs.densities,
                                                               periodic_abfs.potentials,
                                                               target.values,
                                                               grid_data.volume_element)
                  .relative_error;
    }
    const std::vector<SternheimerABFBlochGridChannel>& channels = periodic_abfs.potentials;
    append_chi0_progress_event("channels_ready",
                               0,
                               -1,
                               -1,
                               -1,
                               0,
                               nullptr,
                               -1.0,
                               elapsed_seconds_since(chi0_start_time),
                               "channels=" + std::to_string(num_channels));
    if (GlobalV::MY_RANK == 0)
    {
        write_abfs_channel_diagnostic("STERNHEIMER_ABFS_CHANNELS.dat", channels);
        if (sternheimer_grid_coulomb_diagnostic_enabled(num_channels))
        {
            write_sternheimer_grid_coulomb_diagnostic("STERNHEIMER_GRID_COULOMB.dat",
                                                       periodic_abfs.densities,
                                                       periodic_abfs.potentials,
                                                       grid_data.volume_element);
            if (gamma_qpoint)
            {
                const std::vector<SternheimerABFBlochGridChannel> body_potentials
                    = solve_sternheimer_abf_periodic_full_coulomb(
                        periodic_abfs.densities, grid_data.grid, response_plan.qpoint, 0.0);
                write_sternheimer_grid_coulomb_diagnostic("STERNHEIMER_GRID_COULOMB_BODY.dat",
                                                           periodic_abfs.densities,
                                                           body_potentials,
                                                           grid_data.volume_element);
            }
        }
        if (sternheimer_abfs_diag_only_enabled())
        {
            out << "status abfs_diag_only\n";
            out << "sternheimer_q_index " << response_plan.iq << '\n';
            out << "sternheimer_qpoint " << response_plan.qpoint[0] << ' ' << response_plan.qpoint[1] << ' '
                << response_plan.qpoint[2] << '\n';
            out << "grid " << grid_data.grid.nx << ' ' << grid_data.grid.ny << ' ' << grid_data.grid.nz
                << " size " << grid_data.grid.size() << " dV " << grid_data.volume_element << '\n';
            out << "abfs_channels " << num_channels << '\n';
            out << "perturbation_coulomb_kernel full_periodic_poisson\n";
            out << "periodic_kmesh " << kmesh[0] << ' ' << kmesh[1] << ' ' << kmesh[2] << '\n';
            out << "periodic_gamma_massidda_chi " << massidda_chi << '\n';
            out << "periodic_gamma_coulomb_projection diagnostic_only_physical_poisson\n";
            out << "periodic_gamma_projection_relative_error " << gamma_projection_relative_error << '\n';
            out << "periodic_gamma_limit constant_mode_only_no_headwing\n";
            out << "perturbation_ccp_rmesh_times_used no\n";
        }
    }
    if (sternheimer_abfs_diag_only_enabled())
    {
        return;
    }

    const std::vector<SternheimerFDHamiltonian::Vector> potentials = collect_channel_potentials(channels);
    const std::vector<SternheimerFDHamiltonian::Vector> perturbations_ry
        = scale_potentials(potentials, kHartreeToRydberg);
    const std::vector<SternheimerRPA::AuxiliaryChannel> auxiliary_channels
        = make_chi0_auxiliary_channels(channels);

    std::vector<int> frequency_owners(static_cast<std::size_t>(nfreq), 0);
    std::vector<std::vector<SternheimerRPA::Complex>> chi0_branches(static_cast<std::size_t>(nfreq));
    std::vector<std::vector<SternheimerRPA::Complex>> delta_sos_branches(static_cast<std::size_t>(nfreq));
    std::vector<std::vector<SternheimerRPA::Complex>> delta_pulay_branches(static_cast<std::size_t>(nfreq));
    std::vector<std::vector<SternheimerRPA::Complex>> delta_out_branches(static_cast<std::size_t>(nfreq));
    std::vector<std::vector<SternheimerRPA::Complex>> lcao_sos_branches(static_cast<std::size_t>(nfreq));
    for (int ifrequency = 0; ifrequency != nfreq; ++ifrequency)
    {
        const int owner_rank = use_frequency_mpi
                                   ? SternheimerRPA::frequency_owner_rank(ifrequency,
                                                                          GlobalV::NPROC,
                                                                          frequency_rank_shift)
                                   : 0;
        frequency_owners[static_cast<std::size_t>(ifrequency)] = owner_rank;
        if (use_kpoint_mpi || owner_rank == GlobalV::MY_RANK)
        {
            chi0_branches[static_cast<std::size_t>(ifrequency)].assign(
                static_cast<std::size_t>(num_channels) * static_cast<std::size_t>(num_channels),
                SternheimerRPA::Complex(0.0, 0.0));
            if (write_delta_components)
            {
                delta_sos_branches[static_cast<std::size_t>(ifrequency)]
                    = chi0_branches[static_cast<std::size_t>(ifrequency)];
                delta_pulay_branches[static_cast<std::size_t>(ifrequency)]
                    = chi0_branches[static_cast<std::size_t>(ifrequency)];
                delta_out_branches[static_cast<std::size_t>(ifrequency)]
                    = chi0_branches[static_cast<std::size_t>(ifrequency)];
            }
            if (write_lcao_sos)
            {
                lcao_sos_branches[static_cast<std::size_t>(ifrequency)]
                    = chi0_branches[static_cast<std::size_t>(ifrequency)];
            }
        }
    }

    SternheimerRPA::SolverOptions solver_options;
    solver_options.max_iter = solver_max_iter;
    solver_options.residual_tol = solver_tolerance;
    SternheimerDeltaSubspaceOptions delta_options;
    delta_options.max_virtual_states = PARAM.inp.sternheimer_delta_max_states;
    delta_options.norm_tolerance = PARAM.inp.sternheimer_delta_norm_tol;

    bool all_converged = true;
    int solved_equations = 0;
    double max_solver_relative_residual = 0.0;
    double max_equation_residual_norm = 0.0;
    double max_full_grid_equation_residual_norm = 0.0;
    std::vector<int> target_projector_dimensions(response_plan.kq_pairs.size(), 0);
    std::vector<int> target_delta_dimensions(response_plan.kq_pairs.size(), 0);
    std::vector<double> target_delta_eigenvalue_min(response_plan.kq_pairs.size(), 0.0);
    std::vector<double> target_delta_eigenvalue_max(response_plan.kq_pairs.size(), 0.0);
    std::vector<double> target_delta_grid_hamiltonian_relative_difference(response_plan.kq_pairs.size(), 0.0);
    std::vector<double> target_delta_grid_hamiltonian_max_abs_difference(response_plan.kq_pairs.size(), 0.0);
    std::vector<double> target_lcao_unoccupied_min(response_plan.kq_pairs.size(), 0.0);
    std::vector<double> target_lcao_unoccupied_max(response_plan.kq_pairs.size(), 0.0);
    std::vector<double> target_lcao_occupied_raw_norm_min(response_plan.kq_pairs.size(), 0.0);
    std::vector<double> target_lcao_occupied_raw_norm_max(response_plan.kq_pairs.size(), 0.0);
    std::vector<double> target_lcao_unoccupied_raw_norm_min(response_plan.kq_pairs.size(), 0.0);
    std::vector<double> target_lcao_unoccupied_raw_norm_max(response_plan.kq_pairs.size(), 0.0);
    std::vector<double> target_lcao_occ_unocc_overlap_max(response_plan.kq_pairs.size(), 0.0);

    const std::vector<std::size_t> owned_pair_indices
        = sternheimer_owned_kq_pair_indices(response_plan,
                                            use_kpoint_mpi ? GlobalV::MY_RANK : 0,
                                            use_kpoint_mpi ? kpoint_groups : 1);
    for (const std::size_t pair_index: owned_pair_indices)
    {
        const SternheimerKQPair& pair = response_plan.kq_pairs[pair_index];
        const int source_record_index
            = response_plan.record_index_by_global_k[static_cast<std::size_t>(pair.source_index)];
        const int target_record_index
            = response_plan.record_index_by_global_k[static_cast<std::size_t>(pair.target_index)];
        const SternheimerLCAOOccupiedKPoint& source_record
            = occupied_kpoints[static_cast<std::size_t>(source_record_index)];
        const SternheimerLCAOOccupiedKPoint& target_record
            = occupied_kpoints[static_cast<std::size_t>(target_record_index)];

        append_chi0_progress_event("kpair_start",
                                   0,
                                   -1,
                                   -1,
                                   -1,
                                   solved_equations,
                                   nullptr,
                                   -1.0,
                                   elapsed_seconds_since(chi0_start_time),
                                   "source_k=" + std::to_string(pair.source_index + 1)
                                       + ",target_k=" + std::to_string(pair.target_index + 1));

        const SternheimerSampledLCAOKPoint source
            = sample_sternheimer_lcao_kpoint(ucell,
                                             grid_data.grid,
                                             orbitals,
                                             source_record,
                                             grid_data.volume_element,
                                             PARAM.inp.sternheimer_delta_norm_tol);
        const SternheimerSampledLCAOKPoint target
            = sample_sternheimer_lcao_kpoint(ucell,
                                             grid_data.grid,
                                             orbitals,
                                             target_record,
                                             grid_data.volume_element,
                                             PARAM.inp.sternheimer_delta_norm_tol);
        const SternheimerFDHamiltonian hamiltonian = [&]() {
            if (use_kpoint_mpi)
            {
                SternheimerABACUSFDGridData target_grid_data = grid_data;
                target_grid_data.grid.kpoint = target_record.kpoint;
                auto nonlocal_projector = make_sternheimer_fd_nonlocal_projector_from_unitcell(
                    ucell, target_grid_data.grid, target_grid_data.volume_element);
                return make_sternheimer_fd_hamiltonian_from_local_potential(target_grid_data,
                                                                             kpoint_parallel_full_potential,
                                                                             1.0,
                                                                             std::move(nonlocal_projector));
            }
            return use_frequency_mpi
                       ? make_sternheimer_fd_full_hamiltonian(
                             potential, pw_basis, ucell, 0, 1.0, target_record.kpoint)
                       : make_sternheimer_fd_hamiltonian(
                             potential, pw_basis, ucell, 0, 1.0, target_record.kpoint);
        }();

        std::vector<SternheimerFDHamiltonian::Vector> target_occupied_projector;
        target_occupied_projector.reserve(target.occupied_projector_functions.size());
        for (const SternheimerDeltaGridFunction& function: target.occupied_projector_functions)
        {
            target_occupied_projector.push_back(function.values);
        }
        SternheimerDeltaSubspace delta_subspace;
        target_projector_dimensions[pair_index] = static_cast<int>(target_occupied_projector.size());
        if (use_delta_sternheimer)
        {
            SternheimerDeltaSubspaceOptions pair_delta_options = delta_options;
            pair_delta_options.max_virtual_states = sternheimer_delta_virtual_state_limit(
                delta_options.max_virtual_states,
                static_cast<int>(target.sampled_ao_functions.size()),
                static_cast<int>(target_occupied_projector.size()));
            delta_subspace = build_delta_sternheimer_subspace_by_mode(hamiltonian,
                                                                      target.occupied_projector_functions,
                                                                      target.sampled_ao_functions,
                                                                      grid_data.volume_element,
                                                                      pair_delta_options,
                                                                      delta_a_block_mode);
            if (delta_subspace.virtual_states.empty())
            {
                throw std::runtime_error("Periodic Sternheimer produced no target-sector Delta virtual states.");
            }
            target_delta_dimensions[pair_index] = static_cast<int>(delta_subspace.virtual_states.size());
            const auto delta_minmax = std::minmax_element(
                delta_subspace.virtual_states.begin(),
                delta_subspace.virtual_states.end(),
                [](const SternheimerDeltaVirtualState& lhs, const SternheimerDeltaVirtualState& rhs) {
                    return lhs.eigenvalue < rhs.eigenvalue;
                });
            target_delta_eigenvalue_min[pair_index] = delta_minmax.first->eigenvalue;
            target_delta_eigenvalue_max[pair_index] = delta_minmax.second->eigenvalue;
            target_delta_grid_hamiltonian_relative_difference[pair_index]
                = delta_subspace.full_grid_hamiltonian_relative_difference;
            target_delta_grid_hamiltonian_max_abs_difference[pair_index]
                = delta_subspace.full_grid_hamiltonian_max_abs_difference;
        }

        const int target_occupied_count = static_cast<int>(target.states.eigenvalues.size());
        if (target_occupied_count >= elec_state.ekb.nc)
        {
            throw std::runtime_error("Periodic Sternheimer target k point has no LCAO unoccupied states.");
        }
        double lcao_unoccupied_min = std::numeric_limits<double>::infinity();
        double lcao_unoccupied_max = -std::numeric_limits<double>::infinity();
        for (int ib = target_occupied_count; ib != elec_state.ekb.nc; ++ib)
        {
            const double eigenvalue = elec_state.ekb(target_record.local_k_index, ib);
            lcao_unoccupied_min = std::min(lcao_unoccupied_min, eigenvalue);
            lcao_unoccupied_max = std::max(lcao_unoccupied_max, eigenvalue);
        }
        target_lcao_unoccupied_min[pair_index] = lcao_unoccupied_min;
        target_lcao_unoccupied_max[pair_index] = lcao_unoccupied_max;
        target_lcao_occupied_raw_norm_min[pair_index] = target.occupied_raw_norm_min;
        target_lcao_occupied_raw_norm_max[pair_index] = target.occupied_raw_norm_max;

        std::vector<SternheimerDeltaVirtualState> lcao_virtual_states;
        if (write_lcao_sos)
        {
            lcao_virtual_states.reserve(target.unoccupied_states.wavefunctions.size());
            for (std::size_t ia = 0; ia != target.unoccupied_states.wavefunctions.size(); ++ia)
            {
                SternheimerDeltaVirtualState state;
                state.orbital = target.unoccupied_states.wavefunctions[ia];
                state.residual.assign(state.orbital.size(), SternheimerFDHamiltonian::Complex(0.0, 0.0));
                state.eigenvalue = target.unoccupied_states.eigenvalues[ia];
                lcao_virtual_states.push_back(std::move(state));
            }
            target_lcao_unoccupied_raw_norm_min[pair_index] = target.unoccupied_raw_norm_min;
            target_lcao_unoccupied_raw_norm_max[pair_index] = target.unoccupied_raw_norm_max;
            double max_overlap = 0.0;
            for (const auto& occupied: target.states.wavefunctions)
            {
                for (const auto& virtual_state: lcao_virtual_states)
                {
                    max_overlap = std::max(
                        max_overlap,
                        std::abs(sternheimer_fd_grid_dot(occupied,
                                                        virtual_state.orbital,
                                                        grid_data.volume_element)));
                }
            }
            target_lcao_occ_unocc_overlap_max[pair_index] = max_overlap;
        }

        for (int ifrequency = 0; ifrequency != nfreq; ++ifrequency)
        {
            const int owner_rank = frequency_owners[static_cast<std::size_t>(ifrequency)];
            if (!use_kpoint_mpi && owner_rank != GlobalV::MY_RANK)
            {
                continue;
            }
            const double omega_ry = 2.0 * frequency_grid.omega_ha[static_cast<std::size_t>(ifrequency)];
            std::vector<SternheimerRPA::Complex>& chi0_branch
                = chi0_branches[static_cast<std::size_t>(ifrequency)];
            std::vector<SternheimerRPA::Complex>* delta_sos_branch
                = write_delta_components ? &delta_sos_branches[static_cast<std::size_t>(ifrequency)] : nullptr;
            std::vector<SternheimerRPA::Complex>* delta_pulay_branch
                = write_delta_components ? &delta_pulay_branches[static_cast<std::size_t>(ifrequency)] : nullptr;
            std::vector<SternheimerRPA::Complex>* delta_out_branch
                = write_delta_components ? &delta_out_branches[static_cast<std::size_t>(ifrequency)] : nullptr;
            std::vector<SternheimerRPA::Complex>* lcao_sos_branch
                = write_lcao_sos ? &lcao_sos_branches[static_cast<std::size_t>(ifrequency)] : nullptr;

            for (int ib = 0; ib != static_cast<int>(source.states.wavefunctions.size()); ++ib)
            {
                const double occupation = sternheimer_lcao_weighted_occupation(source_record, ib);
                struct PeriodicChannelEquationResult
                {
                    SternheimerRPA::SolverResult solver;
                    double equation_residual_norm = 0.0;
                    double full_grid_equation_residual_norm = 0.0;
                };
                const std::vector<PeriodicChannelEquationResult> channel_results
                    = run_sternheimer_channel_tasks<PeriodicChannelEquationResult>(
                        num_channels,
                        [&](const int ichannel) {
                            const std::size_t channel_index = static_cast<std::size_t>(ichannel);
                            SternheimerFDHamiltonian::Vector rhs;
                            SternheimerRPA::build_rhs_from_hartree_perturbation(
                                perturbations_ry[channel_index],
                                source.states.wavefunctions[ib],
                                rhs);
                            std::vector<SternheimerFDHamiltonian::Complex> perturbation_matrix_elements;
                            if (use_delta_sternheimer)
                            {
                                perturbation_matrix_elements = delta_sternheimer_perturbation_matrix_elements(
                                    delta_subspace.virtual_states,
                                    perturbations_ry[channel_index],
                                    source.states.wavefunctions[ib],
                                    grid_data.volume_element);
                            }
                            const SternheimerPeriodicLinearResponse response
                                = solve_sternheimer_periodic_linear_response(use_delta_sternheimer,
                                                                             hamiltonian,
                                                                             target_occupied_projector,
                                                                             source.states.eigenvalues[ib],
                                                                             rhs,
                                                                             delta_subspace.virtual_states,
                                                                             perturbation_matrix_elements,
                                                                             omega_ry,
                                                                             grid_data.volume_element,
                                                                             solver_options);
                            SternheimerRPA::accumulate_chi0_branch_column(
                                potentials,
                                source.states.wavefunctions[ib],
                                response.wavefunction,
                                grid_data.volume_element,
                                occupation,
                                ichannel,
                                chi0_branch);
                            if (write_delta_components)
                            {
                                SternheimerRPA::accumulate_chi0_branch_column(
                                    potentials,
                                    source.states.wavefunctions[ib],
                                    response.delta_components.in_sos_wavefunction,
                                    grid_data.volume_element,
                                    occupation,
                                    ichannel,
                                    *delta_sos_branch);
                                SternheimerRPA::accumulate_chi0_branch_column(
                                    potentials,
                                    source.states.wavefunctions[ib],
                                    response.delta_components.in_pulay_wavefunction,
                                    grid_data.volume_element,
                                    occupation,
                                    ichannel,
                                    *delta_pulay_branch);
                                SternheimerRPA::accumulate_chi0_branch_column(
                                    potentials,
                                    source.states.wavefunctions[ib],
                                    response.delta_components.out_wavefunction,
                                    grid_data.volume_element,
                                    occupation,
                                    ichannel,
                                    *delta_out_branch);
                            }
                            if (write_lcao_sos)
                            {
                                const auto lcao_matrix_elements = delta_sternheimer_perturbation_matrix_elements(
                                    lcao_virtual_states,
                                    perturbations_ry[channel_index],
                                    source.states.wavefunctions[ib],
                                    grid_data.volume_element);
                                const auto lcao_response = build_delta_sternheimer_sos_wavefunction(
                                    lcao_virtual_states,
                                    lcao_matrix_elements,
                                    source.states.eigenvalues[ib],
                                    omega_ry);
                                SternheimerRPA::accumulate_chi0_branch_column(
                                    potentials,
                                    source.states.wavefunctions[ib],
                                    lcao_response,
                                    grid_data.volume_element,
                                    occupation,
                                    ichannel,
                                    *lcao_sos_branch);
                            }
                            PeriodicChannelEquationResult result;
                            result.solver = response.solver;
                            result.equation_residual_norm = response.residual_norm;
                            result.full_grid_equation_residual_norm
                                = response.full_grid_equation_residual_norm;
                            return result;
                        });

                for (int ichannel = 0; ichannel != num_channels; ++ichannel)
                {
                    const PeriodicChannelEquationResult& result
                        = channel_results[static_cast<std::size_t>(ichannel)];
                    all_converged = all_converged && result.solver.converged;
                    ++solved_equations;
                    max_solver_relative_residual
                        = std::max(max_solver_relative_residual, result.solver.relative_residual);
                    max_equation_residual_norm
                        = std::max(max_equation_residual_norm, result.equation_residual_norm);
                    max_full_grid_equation_residual_norm
                        = std::max(max_full_grid_equation_residual_norm,
                                   result.full_grid_equation_residual_norm);
                    append_chi0_progress_event("equation",
                                               ifrequency + 1,
                                               owner_rank,
                                               ib,
                                               ichannel,
                                               solved_equations,
                                               &result.solver,
                                               result.equation_residual_norm,
                                               elapsed_seconds_since(chi0_start_time),
                                               "source_k=" + std::to_string(pair.source_index + 1)
                                                   + ",target_k=" + std::to_string(pair.target_index + 1));
                }
            }
        }
        append_chi0_progress_event("kpair_finish",
                                   0,
                                   -1,
                                   -1,
                                   -1,
                                   solved_equations,
                                   nullptr,
                                   -1.0,
                                   elapsed_seconds_since(chi0_start_time),
                                   std::string("mode=") + (use_delta_sternheimer ? "delta" : "standard")
                                       + ",delta_dim=" + std::to_string(delta_subspace.virtual_states.size()));
    }

    if (use_kpoint_mpi)
    {
        for (int ifrequency = 0; ifrequency != nfreq; ++ifrequency)
        {
            reduce_kpoint_parallel_complex_vector(chi0_branches[static_cast<std::size_t>(ifrequency)], true);
            if (write_delta_components)
            {
                reduce_kpoint_parallel_complex_vector(
                    delta_sos_branches[static_cast<std::size_t>(ifrequency)], true);
                reduce_kpoint_parallel_complex_vector(
                    delta_pulay_branches[static_cast<std::size_t>(ifrequency)], true);
                reduce_kpoint_parallel_complex_vector(
                    delta_out_branches[static_cast<std::size_t>(ifrequency)], true);
            }
            if (write_lcao_sos)
            {
                reduce_kpoint_parallel_complex_vector(
                    lcao_sos_branches[static_cast<std::size_t>(ifrequency)], true);
            }
        }
        reduce_kpoint_parallel_int_vector(target_projector_dimensions, true);
        reduce_kpoint_parallel_int_vector(target_delta_dimensions, true);
        reduce_kpoint_parallel_double_vector(target_delta_eigenvalue_min, true);
        reduce_kpoint_parallel_double_vector(target_delta_eigenvalue_max, true);
        reduce_kpoint_parallel_double_vector(target_delta_grid_hamiltonian_relative_difference, true);
        reduce_kpoint_parallel_double_vector(target_delta_grid_hamiltonian_max_abs_difference, true);
        reduce_kpoint_parallel_double_vector(target_lcao_unoccupied_min, true);
        reduce_kpoint_parallel_double_vector(target_lcao_unoccupied_max, true);
        reduce_kpoint_parallel_double_vector(target_lcao_occupied_raw_norm_min, true);
        reduce_kpoint_parallel_double_vector(target_lcao_occupied_raw_norm_max, true);
        reduce_kpoint_parallel_double_vector(target_lcao_unoccupied_raw_norm_min, true);
        reduce_kpoint_parallel_double_vector(target_lcao_unoccupied_raw_norm_max, true);
        reduce_kpoint_parallel_double_vector(target_lcao_occ_unocc_overlap_max, true);
    }
    if (!write_lcao_sos)
    {
        std::fill(target_lcao_unoccupied_raw_norm_min.begin(), target_lcao_unoccupied_raw_norm_min.end(), -1.0);
        std::fill(target_lcao_unoccupied_raw_norm_max.begin(), target_lcao_unoccupied_raw_norm_max.end(), -1.0);
        std::fill(target_lcao_occ_unocc_overlap_max.begin(), target_lcao_occ_unocc_overlap_max.end(), -1.0);
    }

    for (int ifrequency = 0; ifrequency != nfreq; ++ifrequency)
    {
        if ((use_kpoint_mpi && GlobalV::MY_RANK != 0)
            || (!use_kpoint_mpi
                && frequency_owners[static_cast<std::size_t>(ifrequency)] != GlobalV::MY_RANK))
        {
            continue;
        }
        const std::vector<SternheimerRPA::Complex> chi0
            = SternheimerRPA::symmetrize_chi0_imaginary_frequency(
                chi0_branches[static_cast<std::size_t>(ifrequency)], num_channels);
        const SternheimerRPA::Chi0V1Metadata metadata
            = make_chi0_v1_metadata(ucell,
                                    channels,
                                    response_plan.iq,
                                    ifrequency + 1,
                                    frequency_grid.omega_ha[static_cast<std::size_t>(ifrequency)],
                                    frequency_grid.weights_ha[static_cast<std::size_t>(ifrequency)]);
        const std::string data_file = chi0_v1_filename(metadata.iq, metadata.ifrequency);
        SternheimerRPA::write_chi0_v1_file(data_file, metadata, auxiliary_channels, chi0);
        GlobalV::ofs_running << " Sternheimer periodic chi0 v1 output: " << data_file << std::endl;
        if (write_delta_components)
        {
            const auto write_component = [&](const std::string& name,
                                             const std::vector<SternheimerRPA::Complex>& branch) {
                const std::vector<SternheimerRPA::Complex> matrix
                    = SternheimerRPA::symmetrize_chi0_imaginary_frequency(branch, num_channels);
                SternheimerRPA::write_chi0_v1_file(delta_component_v1_filename(name,
                                                                                metadata.iq,
                                                                                metadata.ifrequency),
                                                   metadata,
                                                   auxiliary_channels,
                                                   matrix);
            };
            write_component("in_sos", delta_sos_branches[static_cast<std::size_t>(ifrequency)]);
            write_component("in_pulay", delta_pulay_branches[static_cast<std::size_t>(ifrequency)]);
            write_component("out_grid", delta_out_branches[static_cast<std::size_t>(ifrequency)]);
        }
        if (write_lcao_sos)
        {
            const std::vector<SternheimerRPA::Complex> lcao_sos
                = SternheimerRPA::symmetrize_chi0_imaginary_frequency(
                    lcao_sos_branches[static_cast<std::size_t>(ifrequency)], num_channels);
            SternheimerRPA::write_chi0_v1_file(lcao_sos_v1_filename(metadata.iq, metadata.ifrequency),
                                               metadata,
                                               auxiliary_channels,
                                               lcao_sos);
        }
    }

    reduce_chi0_output_stats(all_converged,
                             solved_equations,
                             max_solver_relative_residual,
                             max_equation_residual_norm,
                             use_parallel_grid_mpi);
#ifdef __MPI
    if (use_parallel_grid_mpi && GlobalV::NPROC > 1)
    {
        MPI_Allreduce(MPI_IN_PLACE,
                      &max_full_grid_equation_residual_norm,
                      1,
                      MPI_DOUBLE,
                      MPI_MAX,
                      MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
    }
#endif

    if (GlobalV::MY_RANK != 0)
    {
        return;
    }
    std::vector<std::pair<std::string, SternheimerRPA::Chi0V1Metadata>> index_entries;
    index_entries.reserve(static_cast<std::size_t>(nfreq));
    for (int ifrequency = 0; ifrequency != nfreq; ++ifrequency)
    {
        const SternheimerRPA::Chi0V1Metadata metadata
            = make_chi0_v1_metadata(ucell,
                                    channels,
                                    response_plan.iq,
                                    ifrequency + 1,
                                    frequency_grid.omega_ha[static_cast<std::size_t>(ifrequency)],
                                    frequency_grid.weights_ha[static_cast<std::size_t>(ifrequency)]);
        index_entries.push_back(
            {chi0_v1_filename(metadata.iq,
                              metadata.ifrequency,
                              frequency_owners[static_cast<std::size_t>(ifrequency)]),
             metadata});
    }
    write_chi0_index_file("v1_sternheimer_chi0_index.dat", index_entries);

    out << "status success\n";
    out << "format v1\n";
    out << "data_files " << index_entries.size() << '\n';
    out << "index_file v1_sternheimer_chi0_index.dat\n";
    out << "sternheimer_q_index " << response_plan.iq << '\n';
    out << "sternheimer_qpoint " << response_plan.qpoint[0] << ' ' << response_plan.qpoint[1] << ' '
        << response_plan.qpoint[2] << '\n';
    out << "sternheimer_kweight_sum " << response_plan.kweight_sum << '\n';
    out << "sternheimer_kq_pairs " << response_plan.kq_pairs.size() << '\n';
    out << "source_global_k target_global_k gx gy gz target_projector_dimension target_delta_dimension "
           "delta_eigenvalue_min_Ry delta_eigenvalue_max_Ry lcao_unoccupied_min_Ry "
           "lcao_unoccupied_max_Ry lcao_occupied_raw_norm_min lcao_occupied_raw_norm_max "
           "lcao_unoccupied_raw_norm_min lcao_unoccupied_raw_norm_max lcao_occ_unocc_overlap_max "
           "delta_grid_hamiltonian_relative_difference delta_grid_hamiltonian_max_abs_difference_Ry\n";
    for (std::size_t pair_index = 0; pair_index != response_plan.kq_pairs.size(); ++pair_index)
    {
        const SternheimerKQPair& pair = response_plan.kq_pairs[pair_index];
        out << pair.source_index + 1 << ' ' << pair.target_index + 1 << ' ' << pair.reciprocal_shift[0] << ' '
            << pair.reciprocal_shift[1] << ' ' << pair.reciprocal_shift[2] << ' '
            << target_projector_dimensions[pair_index] << ' ' << target_delta_dimensions[pair_index] << ' '
            << target_delta_eigenvalue_min[pair_index] << ' ' << target_delta_eigenvalue_max[pair_index] << ' '
            << target_lcao_unoccupied_min[pair_index] << ' ' << target_lcao_unoccupied_max[pair_index] << ' '
            << target_lcao_occupied_raw_norm_min[pair_index] << ' '
            << target_lcao_occupied_raw_norm_max[pair_index] << ' '
            << target_lcao_unoccupied_raw_norm_min[pair_index] << ' '
            << target_lcao_unoccupied_raw_norm_max[pair_index] << ' '
            << target_lcao_occ_unocc_overlap_max[pair_index] << ' '
            << target_delta_grid_hamiltonian_relative_difference[pair_index] << ' '
            << target_delta_grid_hamiltonian_max_abs_difference[pair_index] << '\n';
    }
    out << "grid " << grid_data.grid.nx << ' ' << grid_data.grid.ny << ' ' << grid_data.grid.nz << " size "
        << grid_data.grid.size() << " dV " << grid_data.volume_element << '\n';
    out << "nfreq " << nfreq << '\n';
    out << "frequency_grid_source " << frequency_grid_source << '\n';
    if (use_frequency_grid_file)
    {
        out << "frequency_grid_file " << frequency_grid_file << '\n';
    }
    out << "transition_window_Ha " << transition_window.emin_ha << ' ' << transition_window.emax_ha << '\n';
    out << "pca_threshold " << pca_threshold << '\n';
    out << "ccp_rmesh_times_input " << ccp_rmesh_times << '\n';
    out << "perturbation_coulomb_kernel full_periodic_poisson\n";
    out << "periodic_kmesh " << kmesh[0] << ' ' << kmesh[1] << ' ' << kmesh[2] << '\n';
    out << "periodic_gamma_massidda_chi " << massidda_chi << '\n';
    out << "periodic_gamma_coulomb_projection "
        << (gamma_qpoint ? "diagnostic_only_physical_poisson" : "not_applicable") << '\n';
    out << "periodic_gamma_projection_relative_error " << gamma_projection_relative_error << '\n';
    out << "periodic_gamma_limit "
        << (gamma_qpoint ? "constant_mode_only_no_headwing" : "not_applicable") << '\n';
    out << "perturbation_ccp_rmesh_times_used no\n";
    out << "sternheimer_mode " << (use_delta_sternheimer ? "delta" : "standard") << '\n';
    if (use_delta_sternheimer)
    {
        out << "sternheimer_delta_a_block " << sternheimer_delta_a_block_mode_name(delta_a_block_mode) << '\n';
    }
    out << "delta_component_diagnostic " << (write_delta_components ? "yes" : "no") << '\n';
    out << "lcao_sos_diagnostic " << (write_lcao_sos ? "yes" : "no") << '\n';
    out << "abfs_channels " << num_channels << '\n';
    out << "abfs_max_channels_per_atom " << max_channels << '\n';
    out << "occupied_bands_total " << sternheimer_lcao_total_occupied_bands(occupied_kpoints) << '\n';
    out << "sternheimer_frequency_mpi " << (use_frequency_mpi ? "yes" : "no") << '\n';
    out << "sternheimer_kpoint_mpi " << (use_kpoint_mpi ? "yes" : "no") << '\n';
    out << "sternheimer_kpoint_groups " << (use_kpoint_mpi ? kpoint_groups : 1) << '\n';
    out << "sternheimer_kpoint_group_ranks "
        << (use_kpoint_mpi ? GlobalV::NPROC / kpoint_groups : GlobalV::NPROC) << '\n';
    out << "mpi_ranks " << GlobalV::NPROC << '\n';
    out << "frequency_rank_shift " << frequency_rank_shift << '\n';
    out << "solved_equations " << solved_equations << '\n';
    out << "all_converged " << (all_converged ? "yes" : "no") << '\n';
    out << "max_solver_relative_residual " << max_solver_relative_residual << '\n';
    out << "max_equation_residual_norm " << max_equation_residual_norm << '\n';
    out << "max_full_grid_equation_residual_norm " << max_full_grid_equation_residual_norm << '\n';
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
                const std::vector<double> perturbation_ry = scale_potential(channel.potential_r, kHartreeToRydberg);
                SternheimerFDHamiltonian::Vector rhs;
                SternheimerRPA::build_rhs_from_hartree_perturbation(perturbation_ry,
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

void run_sternheimer_abacus_chi0_output_impl(
    const elecstate::Potential& potential,
    const ModulePW::PW_Basis& pw_basis,
    const UnitCell& ucell,
    const elecstate::ElecState& elec_state,
    const std::string& output_dir,
    const LCAO_Orbitals* lcao_orbitals,
    const std::vector<SternheimerLCAOOccupiedKPoint>* lcao_occupied_kpoints,
    const std::array<int, 3>* lcao_kmesh)
{
    if (!PARAM.inp.out_sternheimer_librpa)
    {
        return;
    }

    const bool use_frequency_mpi = PARAM.inp.sternheimer_frequency_mpi;
    const bool use_kpoint_mpi
        = PARAM.inp.sternheimer_q_index > 0 && PARAM.globalv.kpar_lcao > 1;
    const bool use_parallel_response_mpi = use_frequency_mpi || use_kpoint_mpi;
    if (!use_parallel_response_mpi && GlobalV::MY_RANK != 0)
    {
        return;
    }

    const std::string status_path = chi0_status_path(output_dir);
    std::ofstream out;
    if (GlobalV::MY_RANK == 0)
    {
        out.open(status_path.c_str(), std::ios::out | std::ios::trunc);
        if (!out)
        {
            GlobalV::ofs_running << " Sternheimer chi0 output: failed to open " << status_path << std::endl;
#ifdef __MPI
            if (use_parallel_response_mpi && GlobalV::NPROC > 1)
            {
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
#endif
            return;
        }
        out << std::setprecision(16);
        out << "# ABACUS Sternheimer chi0 output for LibRPA\n";
    }

    try
    {
        reset_chi0_progress_file();
        const auto chi0_start_time = std::chrono::steady_clock::now();
        append_chi0_progress_event("enter",
                                   0,
                                   -1,
                                   -1,
                                   -1,
                                   0,
                                   nullptr,
                                   -1.0,
                                   elapsed_seconds_since(chi0_start_time),
                                   std::string(use_frequency_mpi ? "frequency_mpi=yes" : "frequency_mpi=no")
                                       + ","
                                       + (use_kpoint_mpi ? "kpoint_mpi=yes" : "kpoint_mpi=no"));
        if (PARAM.inp.out_librpa_reader_version != 1)
        {
            throw std::runtime_error("out_sternheimer_librpa requires out_librpa_reader_version=1.");
        }
        if (GlobalV::NPROC != 1 && !use_parallel_response_mpi)
        {
            throw std::runtime_error(
                "Sternheimer chi0 output with multiple MPI ranks requires sternheimer_frequency_mpi=true.");
        }
        if (use_frequency_mpi && use_kpoint_mpi)
        {
            throw std::runtime_error(
                "Periodic Sternheimer does not yet support nested frequency MPI and k-point MPI.");
        }
        if (use_kpoint_mpi
            && (PARAM.globalv.kpar_lcao != GlobalV::NPROC
                || PARAM.globalv.kpar_lcao > elec_state.wg.nr))
        {
            throw std::runtime_error(
                "The first Sternheimer k-point MPI implementation requires NPROC=kpar and kpar<=number of k points.");
        }
        if (use_parallel_response_mpi && GlobalV::NPROC > 1 && pw_basis.poolnproc != GlobalV::NPROC)
        {
            throw std::runtime_error(
                "Parallel Sternheimer currently requires all MPI ranks to be in the same real-space pool.");
        }
        if (elec_state.ekb.nc <= 0 || elec_state.wg.nc <= 0)
        {
            throw std::runtime_error("ABACUS DFT eigenvalues or occupations are not available.");
        }

        const bool use_lcao_zero_order = lcao_occupied_kpoints != nullptr;
        if ((lcao_orbitals == nullptr) != (lcao_occupied_kpoints == nullptr))
        {
            throw std::runtime_error("Sternheimer LCAO zero-order input is incomplete.");
        }

        if (PARAM.inp.sternheimer_q_index > 0)
        {
            if (!use_lcao_zero_order)
            {
                throw std::runtime_error("A nonzero sternheimer_q_index requires LCAO k-resolved zero-order states.");
            }
            if (lcao_kmesh == nullptr)
            {
                throw std::runtime_error("Periodic Sternheimer requires Monkhorst-Pack dimensions.");
            }
            validate_sternheimer_lcao_occupied_kpoints(*lcao_occupied_kpoints,
                                                       elec_state.wg.nr,
                                                       elec_state.wg.nr,
                                                       PARAM.inp.nspin,
                                                       PARAM.globalv.nlocal);
            run_sternheimer_periodic_lcao_chi0_output(potential,
                                                      pw_basis,
                                                      ucell,
                                                      elec_state,
                                                      *lcao_orbitals,
                                                      *lcao_occupied_kpoints,
                                                      *lcao_kmesh,
                                                      out,
                                                      use_frequency_mpi,
                                                      chi0_start_time);
            if (GlobalV::MY_RANK == 0)
            {
                GlobalV::ofs_running << " Sternheimer periodic chi0 status: " << status_path << std::endl;
            }
            return;
        }

        std::vector<const SternheimerLCAOOccupiedKPoint*> response_kpoints(1, nullptr);
        if (use_lcao_zero_order)
        {
            validate_sternheimer_lcao_occupied_kpoints(
                *lcao_occupied_kpoints,
                elec_state.wg.nr,
                elec_state.wg.nr,
                PARAM.inp.nspin,
                PARAM.globalv.nlocal);
            response_kpoints
                = select_sternheimer_gamma_spin_records(*lcao_occupied_kpoints, PARAM.inp.nspin);
        }

        const int occupied_count = use_lcao_zero_order
                                       ? sternheimer_lcao_total_occupied_bands(*lcao_occupied_kpoints)
                                       : occupied_band_count(elec_state, 0);
        if (occupied_count <= 0)
        {
            throw std::runtime_error("No occupied DFT bands are available for Sternheimer chi0 output.");
        }
        const int requested_bands = positive_int_from_env(kBandsEnv, occupied_count);
        const int num_bands = use_lcao_zero_order ? occupied_count
                                                  : std::min(requested_bands, elec_state.ekb.nc);
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
        const int default_frequency_rank_shift = use_frequency_mpi && GlobalV::NPROC > 1 ? 1 : 0;
        const int frequency_rank_shift = int_from_env(kFrequencyRankShiftEnv, default_frequency_rank_shift);
        const std::string frequency_grid_file = PARAM.inp.sternheimer_frequency_grid_file;
        const bool use_delta_sternheimer = PARAM.inp.sternheimer_delta;
        const SternheimerDeltaABlockMode delta_a_block_mode = delta_a_block_mode_from_env();
        if (use_lcao_zero_order)
        {
            if (!use_delta_sternheimer)
            {
                throw std::runtime_error("Sternheimer LCAO zero-order input currently requires sternheimer_delta=true.");
            }
            if ((PARAM.inp.nspin != 1 && PARAM.inp.nspin != 2)
                || elec_state.ekb.nr != elec_state.wg.nr)
            {
                throw std::runtime_error(
                    "Sternheimer LCAO zero-order input currently supports Gamma-point nspin=1 or nspin=2 calculations.");
            }
            if (std::any_of(response_kpoints.begin(), response_kpoints.end(), [](const auto* record) {
                    return record->coefficients.empty();
                }))
            {
                throw std::runtime_error("Sternheimer LCAO response record has no occupied coefficients.");
            }
        }

        SternheimerRPA::TransitionEnergyWindow transition_window;
        transition_window.emin_ha = std::numeric_limits<double>::infinity();
        transition_window.emax_ha = 0.0;
        for (const SternheimerLCAOOccupiedKPoint* response_kpoint: response_kpoints)
        {
            const int response_k_index = response_kpoint == nullptr ? 0 : response_kpoint->local_k_index;
            const std::vector<double> eigenvalues_ry
                = eigenvalues_ry_from_elec_state(elec_state, response_k_index);
            const std::vector<double> occupations
                = occupations_from_elec_state(elec_state, response_k_index);
            const SternheimerRPA::TransitionEnergyWindow spin_window
                = SternheimerRPA::transition_energy_window_from_eigenvalues_ry(eigenvalues_ry, occupations);
            transition_window.emin_ha = std::min(transition_window.emin_ha, spin_window.emin_ha);
            transition_window.emax_ha = std::max(transition_window.emax_ha, spin_window.emax_ha);
        }
        const bool use_frequency_grid_file = !frequency_grid_file.empty();
        const std::string frequency_grid_source = use_frequency_grid_file ? "file" : "greenx_minimax";
        const SternheimerRPA::FrequencyGrid frequency_grid
            = use_frequency_grid_file
                  ? SternheimerRPA::read_frequency_grid_file(frequency_grid_file, nfreq)
                  : SternheimerRPA::generate_greenx_minimax_frequency_grid(nfreq,
                                                                           transition_window.emin_ha,
                                                                           transition_window.emax_ha);
        append_chi0_progress_event("frequency_grid",
                                   0,
                                   -1,
                                   -1,
                                   -1,
                                   0,
                                   nullptr,
                                   -1.0,
                                   elapsed_seconds_since(chi0_start_time),
                                   "rank_shift=" + std::to_string(frequency_rank_shift));

        const SternheimerABACUSFDGridData grid_data
            = use_frequency_mpi ? make_sternheimer_fd_full_grid(pw_basis) : make_sternheimer_fd_grid(pw_basis);
        const std::vector<SternheimerABFGridChannel> channels
            = build_abfs_ccp_grid_channels(ucell, grid_data.grid, -1, pca_threshold, ccp_rmesh_times);
        if (channels.empty())
        {
            throw std::runtime_error("No ABFS CCP perturbation channels were generated.");
        }

        const int num_channels = static_cast<int>(channels.size());
        append_chi0_progress_event("channels_ready",
                                   0,
                                   -1,
                                   -1,
                                   -1,
                                   0,
                                   nullptr,
                                   -1.0,
                                   elapsed_seconds_since(chi0_start_time),
                                   "num_channels=" + std::to_string(num_channels));
        if (GlobalV::MY_RANK == 0)
        {
            write_abfs_channel_diagnostic("STERNHEIMER_ABFS_CHANNELS.dat", channels);
            if (sternheimer_abfs_diag_only_enabled())
            {
                out << "status abfs_diag_only\n";
                out << "grid " << grid_data.grid.nx << ' ' << grid_data.grid.ny << ' ' << grid_data.grid.nz
                    << " size " << grid_data.grid.size() << " dV " << grid_data.volume_element << '\n';
                out << "frequency_grid_source " << frequency_grid_source << '\n';
                out << "nfreq " << nfreq << '\n';
                out << "abfs_channels " << num_channels << '\n';
            }
        }
        if (sternheimer_abfs_diag_only_enabled())
        {
            return;
        }

        std::vector<SternheimerDeltaGridFunction> sampled_ao_functions;
        if (use_lcao_zero_order)
        {
            sampled_ao_functions
                = build_lcao_candidate_grid_functions(ucell, grid_data.grid, lcao_orbitals);
            if (sampled_ao_functions.empty())
            {
                throw std::runtime_error("Sternheimer LCAO zero-order input found no sampled AO functions.");
            }
        }

        const std::size_t response_count = response_kpoints.size();
        std::vector<SternheimerFDHamiltonian> hamiltonians;
        std::vector<SternheimerFDZeroOrderStates> states_by_response(response_count);
        std::vector<std::vector<SternheimerDeltaGridFunction>> lcao_projector_functions_by_response(response_count);
        std::vector<std::vector<SternheimerFDHamiltonian::Vector>> occupied_by_response(response_count);
        std::vector<std::vector<SternheimerFDHamiltonian::Vector>> occupied_projector_by_response(response_count);
        std::vector<SternheimerDeltaSubspace> delta_subspaces(response_count);
        hamiltonians.reserve(response_count);

        for (std::size_t response_index = 0; response_index != response_count; ++response_index)
        {
            const SternheimerLCAOOccupiedKPoint* response_kpoint = response_kpoints[response_index];
            const int response_k_index = response_kpoint == nullptr ? 0 : response_kpoint->local_k_index;
            const int response_spin_index = response_kpoint == nullptr ? 0 : response_kpoint->spin_index;
            hamiltonians.push_back(
                use_frequency_mpi
                    ? make_sternheimer_fd_full_hamiltonian(
                          potential, pw_basis, ucell, response_spin_index, 1.0)
                    : make_sternheimer_fd_hamiltonian(
                          potential, pw_basis, ucell, response_spin_index, 1.0));
            append_chi0_progress_event("hamiltonian_ready",
                                       0,
                                       -1,
                                       -1,
                                       -1,
                                       0,
                                       nullptr,
                                       -1.0,
                                       elapsed_seconds_since(chi0_start_time),
                                       "spin=" + std::to_string(response_spin_index + 1)
                                           + " grid_size=" + std::to_string(grid_data.grid.size()));

            SternheimerFDZeroOrderStates& states = states_by_response[response_index];
            append_chi0_progress_event("zero_order_start",
                                       0,
                                       -1,
                                       -1,
                                       -1,
                                       0,
                                       nullptr,
                                       -1.0,
                                       elapsed_seconds_since(chi0_start_time),
                                       std::string("spin=") + std::to_string(response_spin_index + 1) + " "
                                           + (use_lcao_zero_order
                                                  ? "lcao_sample"
                                                  : (GlobalV::MY_RANK == 0 ? "fd_solve" : "wait")));
            if (use_lcao_zero_order)
            {
                const int spin_occupied_count = static_cast<int>(response_kpoint->coefficients.size());
                std::vector<SternheimerDeltaGridFunction> lcao_occupied_functions;
                lcao_occupied_functions.reserve(static_cast<std::size_t>(spin_occupied_count));
                states.eigenvalues.reserve(static_cast<std::size_t>(spin_occupied_count));
                states.wavefunctions.reserve(static_cast<std::size_t>(spin_occupied_count));
                states.residual_norms.reserve(static_cast<std::size_t>(spin_occupied_count));
                for (int ib = 0; ib != spin_occupied_count; ++ib)
                {
                    const auto& coefficients
                        = response_kpoint->coefficients[static_cast<std::size_t>(ib)];
                    if (coefficients.size() != sampled_ao_functions.size())
                    {
                        throw std::runtime_error(
                            "Sternheimer LCAO coefficient basis size does not match sampled AO functions.");
                    }
                    SternheimerDeltaGridFunction occupied_function
                        = linear_combination_delta_sternheimer_grid_functions(sampled_ao_functions, coefficients);
                    const double norm = sternheimer_fd_grid_norm(occupied_function.values, grid_data.volume_element);
                    if (norm <= PARAM.inp.sternheimer_delta_norm_tol)
                    {
                        throw std::runtime_error("Sternheimer sampled LCAO occupied function has zero norm.");
                    }
                    const SternheimerFDHamiltonian::Complex inverse_norm(1.0 / norm, 0.0);
                    for (auto& value: occupied_function.values)
                    {
                        value *= inverse_norm;
                    }
                    for (auto& gradient: occupied_function.gradients)
                    {
                        for (auto& value: gradient)
                        {
                            value *= inverse_norm;
                        }
                    }
                    states.eigenvalues.push_back(response_kpoint->eigenvalues[static_cast<std::size_t>(ib)]);
                    states.wavefunctions.push_back(occupied_function.values);
                    states.residual_norms.push_back(0.0);
                    lcao_occupied_functions.push_back(std::move(occupied_function));
                }
                lcao_projector_functions_by_response[response_index]
                    = orthonormalize_delta_sternheimer_grid_functions(
                        lcao_occupied_functions,
                        grid_data.volume_element,
                        PARAM.inp.sternheimer_delta_norm_tol);
            }
            else
            {
                if (!use_frequency_mpi || GlobalV::MY_RANK == 0)
                {
                    states = solve_fd_zero_order_auto(hamiltonians[response_index],
                                                      num_bands,
                                                      grid_data.volume_element,
                                                      max_dense_size,
                                                      lanczos_max_subspace_size);
                }
                broadcast_zero_order_states(states, grid_data.grid.size(), use_frequency_mpi);
            }
            append_chi0_progress_event("zero_order_ready",
                                       0,
                                       -1,
                                       -1,
                                       -1,
                                       0,
                                       nullptr,
                                       -1.0,
                                       elapsed_seconds_since(chi0_start_time),
                                       "spin=" + std::to_string(response_spin_index + 1) + " source="
                                           + (use_lcao_zero_order ? "lcao_ks" : "fd_grid")
                                           + " nstates=" + std::to_string(states.wavefunctions.size()));
            occupied_by_response[response_index]
                = occupied_wavefunctions_from_states(states, elec_state, response_k_index);
            if (occupied_by_response[response_index].empty())
            {
                throw std::runtime_error("No occupied zero-order states are available for Sternheimer chi0 output.");
            }
            occupied_projector_by_response[response_index] = occupied_by_response[response_index];
            if (use_lcao_zero_order)
            {
                occupied_projector_by_response[response_index].clear();
                occupied_projector_by_response[response_index].reserve(
                    lcao_projector_functions_by_response[response_index].size());
                for (const SternheimerDeltaGridFunction& function:
                     lcao_projector_functions_by_response[response_index])
                {
                    occupied_projector_by_response[response_index].push_back(function.values);
                }
            }
        }

        const std::vector<std::vector<double>> potentials = collect_channel_potentials(channels);
        // ABFS Coulomb potentials are in Ha units; the FD Hamiltonian and omega are in Ry.
        // Keep Ha potentials for M=V chi0 V output, but use Ry perturbations in the linear equation.
        const std::vector<std::vector<double>> perturbations_ry = scale_potentials(potentials, kHartreeToRydberg);

        if (use_delta_sternheimer)
        {
            append_chi0_progress_event("delta_subspace_start",
                                       0,
                                       -1,
                                       -1,
                                       -1,
                                       0,
                                       nullptr,
                                       -1.0,
                                       elapsed_seconds_since(chi0_start_time),
                                       "candidate_orbitals");
            std::vector<SternheimerDeltaGridFunction> loaded_candidate_functions;
            const std::vector<SternheimerDeltaGridFunction>* candidate_functions = &sampled_ao_functions;
            if (candidate_functions->empty())
            {
                loaded_candidate_functions = build_lcao_candidate_grid_functions(ucell, grid_data.grid);
                candidate_functions = &loaded_candidate_functions;
            }
            if (candidate_functions->empty())
            {
                throw std::runtime_error("Sternheimer delta mode found no sampled LCAO candidate orbitals.");
            }

            for (std::size_t response_index = 0; response_index != response_count; ++response_index)
            {
                const SternheimerLCAOOccupiedKPoint* response_kpoint = response_kpoints[response_index];
                const int response_spin_index = response_kpoint == nullptr ? 0 : response_kpoint->spin_index;
                std::vector<SternheimerDeltaGridFunction> fd_occupied_functions;
                const std::vector<SternheimerDeltaGridFunction>* occupied_functions
                    = &lcao_projector_functions_by_response[response_index];
                if (occupied_functions->empty())
                {
                    fd_occupied_functions.reserve(occupied_by_response[response_index].size());
                    for (const SternheimerFDHamiltonian::Vector& occupied_wavefunction:
                         occupied_by_response[response_index])
                    {
                        fd_occupied_functions.push_back(make_delta_sternheimer_grid_function_with_fd_gradients(
                            occupied_wavefunction, grid_data.grid));
                    }
                    occupied_functions = &fd_occupied_functions;
                }

                SternheimerDeltaSubspaceOptions delta_options;
                delta_options.max_virtual_states = sternheimer_delta_virtual_state_limit(
                    PARAM.inp.sternheimer_delta_max_states,
                    static_cast<int>(candidate_functions->size()),
                    static_cast<int>(occupied_functions->size()));
                delta_options.norm_tolerance = PARAM.inp.sternheimer_delta_norm_tol;
                delta_subspaces[response_index]
                    = build_delta_sternheimer_subspace_by_mode(hamiltonians[response_index],
                                                               *occupied_functions,
                                                               *candidate_functions,
                                                               grid_data.volume_element,
                                                               delta_options,
                                                               delta_a_block_mode);
                if (delta_subspaces[response_index].virtual_states.empty())
                {
                    throw std::runtime_error("Sternheimer delta mode produced no fixed virtual states.");
                }
                append_chi0_progress_event(
                    "delta_subspace_ready",
                    0,
                    -1,
                    -1,
                    -1,
                    0,
                    nullptr,
                    -1.0,
                    elapsed_seconds_since(chi0_start_time),
                    "spin=" + std::to_string(response_spin_index + 1)
                        + " nvirtual="
                        + std::to_string(delta_subspaces[response_index].virtual_states.size()));
            }
        }

        SternheimerRPA::SolverOptions solver_options;
        solver_options.max_iter = solver_max_iter;
        solver_options.residual_tol = solver_tolerance;

        bool all_converged = true;
        int solved_equations = 0;
        double max_solver_relative_residual = 0.0;
        double max_equation_residual_norm = 0.0;
        const std::vector<SternheimerRPA::AuxiliaryChannel> auxiliary_channels = make_chi0_auxiliary_channels(channels);

        for (int ifrequency = 0; ifrequency != nfreq; ++ifrequency)
        {
            const int owner_rank = use_frequency_mpi
                                       ? SternheimerRPA::frequency_owner_rank(ifrequency,
                                                                              GlobalV::NPROC,
                                                                              frequency_rank_shift)
                                       : 0;
            if (use_frequency_mpi && owner_rank != GlobalV::MY_RANK)
            {
                continue;
            }

            std::vector<SternheimerRPA::Complex> chi0_branch(
                static_cast<std::size_t>(num_channels) * static_cast<std::size_t>(num_channels),
                SternheimerRPA::Complex(0.0, 0.0));

            const double omega_ha = frequency_grid.omega_ha[static_cast<std::size_t>(ifrequency)];
            const double omega_ry = 2.0 * omega_ha;
            const auto frequency_start_time = std::chrono::steady_clock::now();
            append_chi0_progress_event("frequency_start",
                                       ifrequency + 1,
                                       owner_rank,
                                       -1,
                                       -1,
                                       solved_equations,
                                       nullptr,
                                       -1.0,
                                       elapsed_seconds_since(chi0_start_time),
                                       "omega_Ha=" + std::to_string(omega_ha));

            for (std::size_t response_index = 0; response_index != response_count; ++response_index)
            {
                const SternheimerLCAOOccupiedKPoint* response_kpoint = response_kpoints[response_index];
                const int response_k_index = response_kpoint == nullptr ? 0 : response_kpoint->local_k_index;
                const int response_spin_index = response_kpoint == nullptr ? 0 : response_kpoint->spin_index;
                const SternheimerFDHamiltonian& hamiltonian = hamiltonians[response_index];
                const SternheimerFDZeroOrderStates& states = states_by_response[response_index];
                const auto& occupied = occupied_by_response[response_index];
                const auto& occupied_projector = occupied_projector_by_response[response_index];
                const SternheimerDeltaSubspace& delta_subspace = delta_subspaces[response_index];

                for (int ib = 0; ib != static_cast<int>(states.wavefunctions.size()); ++ib)
                {
                    const double occupation
                        = use_lcao_zero_order
                              ? sternheimer_lcao_weighted_occupation(*response_kpoint, ib)
                              : elec_state.wg(response_k_index, ib);
                    if (occupation <= 1.0e-8)
                    {
                        continue;
                    }

                    for (int ichannel = 0; ichannel != num_channels; ++ichannel)
                    {
                        const std::size_t channel_index = static_cast<std::size_t>(ichannel);
                        SternheimerFDHamiltonian::Vector rhs;
                        SternheimerRPA::build_rhs_from_hartree_perturbation(perturbations_ry[channel_index],
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
                                    perturbations_ry[channel_index],
                                    states.wavefunctions[ib],
                                    grid_data.volume_element);
                            const SternheimerDeltaLinearResponse response
                                = solve_delta_sternheimer_linear_response(hamiltonian,
                                                                          occupied_projector,
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
                        append_chi0_progress_event("equation",
                                                   ifrequency + 1,
                                                   owner_rank,
                                                   ib,
                                                   ichannel,
                                                   solved_equations,
                                                   &solver_result,
                                                   equation_residual_norm,
                                                   elapsed_seconds_since(chi0_start_time),
                                                   "spin=" + std::to_string(response_spin_index + 1) + " "
                                                       + (use_delta_sternheimer ? "delta" : "standard"));
                    }
                }
            }

            const std::vector<SternheimerRPA::Complex> chi0
                = SternheimerRPA::symmetrize_chi0_imaginary_frequency(chi0_branch, num_channels);
            const SternheimerRPA::Chi0V1Metadata metadata
                = make_chi0_v1_metadata(ucell,
                                        channels,
                                        1,
                                        ifrequency + 1,
                                        omega_ha,
                                        frequency_grid.weights_ha[static_cast<std::size_t>(ifrequency)]);
            const std::string data_file = chi0_v1_filename(metadata.iq, metadata.ifrequency);
            SternheimerRPA::write_chi0_v1_file(data_file, metadata, auxiliary_channels, chi0);
            GlobalV::ofs_running << " Sternheimer chi0 v1 output: " << data_file << std::endl;
            append_chi0_progress_event("frequency_finish",
                                       ifrequency + 1,
                                       owner_rank,
                                       -1,
                                       -1,
                                       solved_equations,
                                       nullptr,
                                       -1.0,
                                       elapsed_seconds_since(chi0_start_time),
                                       "elapsed_freq_s=" + std::to_string(
                                                               elapsed_seconds_since(frequency_start_time)));
        }

        reduce_chi0_output_stats(all_converged,
                                 solved_equations,
                                 max_solver_relative_residual,
                                 max_equation_residual_norm,
                                 use_frequency_mpi);

#ifdef __MPI
        if (use_parallel_response_mpi && GlobalV::NPROC > 1)
        {
            MPI_Barrier(MPI_COMM_WORLD);
        }
#endif

        std::vector<std::pair<std::string, SternheimerRPA::Chi0V1Metadata>> index_entries;
        if (GlobalV::MY_RANK == 0)
        {
            index_entries.reserve(frequency_grid.omega_ha.size());
            for (int ifrequency = 0; ifrequency != nfreq; ++ifrequency)
            {
                const SternheimerRPA::Chi0V1Metadata metadata
                    = make_chi0_v1_metadata(ucell,
                                            channels,
                                            1,
                                            ifrequency + 1,
                                            frequency_grid.omega_ha[static_cast<std::size_t>(ifrequency)],
                                            frequency_grid.weights_ha[static_cast<std::size_t>(ifrequency)]);
                const int owner_rank = use_frequency_mpi
                                           ? SternheimerRPA::frequency_owner_rank(ifrequency,
                                                                                  GlobalV::NPROC,
                                                                                  frequency_rank_shift)
                                           : 0;
                index_entries.push_back({chi0_v1_filename(metadata.iq, metadata.ifrequency, owner_rank), metadata});
            }
            write_chi0_index_file("v1_sternheimer_chi0_index.dat", index_entries);
        }

        const int grid_size = grid_data.grid.nx * grid_data.grid.ny * grid_data.grid.nz;
        if (GlobalV::MY_RANK != 0)
        {
            return;
        }

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
        out << "sternheimer_frequency_mpi " << (use_frequency_mpi ? "yes" : "no") << '\n';
        out << "mpi_ranks " << GlobalV::NPROC << '\n';
        out << "frequency_rank_shift " << frequency_rank_shift << '\n';
        out << "progress_file_pattern STERNHEIMER_CHI0_PROGRESS_rank*.dat\n";
        out << "ifrequency omega_Ha weight_Ha omega_Ry data_file\n";
        for (const auto& entry: index_entries)
        {
            const SternheimerRPA::Chi0V1Metadata& metadata = entry.second;
            out << metadata.ifrequency << ' ' << metadata.omega << ' ' << metadata.weight << ' '
                << 2.0 * metadata.omega << ' ' << entry.first << '\n';
        }
        out << "pca_threshold " << pca_threshold << '\n';
        out << "ccp_rmesh_times " << ccp_rmesh_times << '\n';
        out << "sternheimer_zero_order_source " << (use_lcao_zero_order ? "lcao_ks" : "fd_grid") << '\n';
        out << "sternheimer_response_records " << response_count << '\n';
        out << "sternheimer_response_spin_channels";
        for (const SternheimerLCAOOccupiedKPoint* response_kpoint: response_kpoints)
        {
            out << ' ' << (response_kpoint == nullptr ? 1 : response_kpoint->spin_index + 1);
        }
        out << '\n';
        if (response_count == 1)
        {
            const SternheimerLCAOOccupiedKPoint* response_kpoint = response_kpoints.front();
            out << "sternheimer_response_spin_channel "
                << (response_kpoint == nullptr ? 1 : response_kpoint->spin_index + 1) << '\n';
            out << "sternheimer_response_local_k_index "
                << (response_kpoint == nullptr ? 1 : response_kpoint->local_k_index + 1) << '\n';
            if (response_kpoint != nullptr)
            {
                out << "sternheimer_response_global_k_index " << response_kpoint->global_k_index + 1 << '\n';
                out << "sternheimer_response_kpoint " << response_kpoint->kpoint[0] << ' '
                    << response_kpoint->kpoint[1] << ' ' << response_kpoint->kpoint[2] << '\n';
                out << "sternheimer_response_kweight " << response_kpoint->kweight << '\n';
            }
        }
        std::size_t occupied_bands_total = 0;
        std::size_t occupied_projector_total = 0;
        for (std::size_t response_index = 0; response_index != response_count; ++response_index)
        {
            occupied_bands_total += occupied_by_response[response_index].size();
            occupied_projector_total += occupied_projector_by_response[response_index].size();
        }
        out << "occupied_bands " << occupied_bands_total << '\n';
        out << "occupied_projector_dimension " << occupied_projector_total << '\n';
        out << "abfs_channels " << num_channels << '\n';
        out << "sternheimer_delta " << (use_delta_sternheimer ? "yes" : "no") << '\n';
        if (use_delta_sternheimer)
        {
            out << "sternheimer_delta_backend " << sternheimer_delta_a_block_mode_name(delta_a_block_mode) << '\n';
            out << "sternheimer_delta_max_states " << PARAM.inp.sternheimer_delta_max_states << '\n';
            out << "sternheimer_delta_norm_tol " << PARAM.inp.sternheimer_delta_norm_tol << '\n';
            std::size_t virtual_states_total = 0;
            int accepted_candidates_total = 0;
            int discarded_candidates_total = 0;
            for (std::size_t response_index = 0; response_index != response_count; ++response_index)
            {
                virtual_states_total += delta_subspaces[response_index].virtual_states.size();
                accepted_candidates_total += delta_subspaces[response_index].accepted_candidates;
                discarded_candidates_total += delta_subspaces[response_index].discarded_candidates;
                out << "sternheimer_delta_spin " << response_index + 1 << " virtual_states "
                    << delta_subspaces[response_index].virtual_states.size() << " accepted_candidates "
                    << delta_subspaces[response_index].accepted_candidates << " discarded_candidates "
                    << delta_subspaces[response_index].discarded_candidates << '\n';
            }
            out << "sternheimer_delta_virtual_states " << virtual_states_total << '\n';
            out << "sternheimer_delta_accepted_candidates " << accepted_candidates_total << '\n';
            out << "sternheimer_delta_discarded_candidates " << discarded_candidates_total << '\n';
        }
        out << "solved_equations " << solved_equations << '\n';
        out << "all_converged " << (all_converged ? "yes" : "no") << '\n';
        out << "max_solver_relative_residual " << max_solver_relative_residual << '\n';
        out << "max_equation_residual_norm " << max_equation_residual_norm << '\n';
        GlobalV::ofs_running << " Sternheimer chi0 status: " << status_path << std::endl;
    }
    catch (const std::exception& error)
    {
        {
            std::ofstream rank_failure(chi0_failure_filename().c_str(), std::ios::out | std::ios::trunc);
            rank_failure << "rank " << GlobalV::MY_RANK << '\n';
            rank_failure << "reason " << error.what() << '\n';
        }
        if (GlobalV::MY_RANK == 0 && out)
        {
            out << "status failed\n";
            out << "reason " << error.what() << '\n';
        }
        GlobalV::ofs_running << " Sternheimer chi0 output failed: " << error.what() << std::endl;
        GlobalV::ofs_running << " Sternheimer chi0 status: " << status_path << std::endl;
        GlobalV::ofs_running.flush();
#ifdef __MPI
        if (use_parallel_response_mpi && GlobalV::NPROC > 1)
        {
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
#endif
    }
}

void run_sternheimer_abacus_chi0_output(const elecstate::Potential& potential,
                                        const ModulePW::PW_Basis& pw_basis,
                                        const UnitCell& ucell,
                                        const elecstate::ElecState& elec_state,
                                        const std::string& output_dir)
{
    run_sternheimer_abacus_chi0_output_impl(
        potential, pw_basis, ucell, elec_state, output_dir, nullptr, nullptr, nullptr);
}

void run_sternheimer_abacus_lcao_chi0_output(
    const elecstate::Potential& potential,
    const ModulePW::PW_Basis& pw_basis,
    const UnitCell& ucell,
    const elecstate::ElecState& elec_state,
    const LCAO_Orbitals& orbitals,
    const std::vector<SternheimerLCAOOccupiedKPoint>& occupied_kpoints,
    const std::array<int, 3>& kmesh,
    const std::string& output_dir)
{
    run_sternheimer_abacus_chi0_output_impl(potential,
                                            pw_basis,
                                            ucell,
                                            elec_state,
                                            output_dir,
                                            &orbitals,
                                            &occupied_kpoints,
                                            &kmesh);
}

} // namespace ModuleRI
