#include "source_lcao/module_ri/sternheimer_abacus_st_smoke.h"
#include "source_lcao/module_ri/singular_value.h"

#include "source_basis/module_ao/ORB_read.h"
#include "source_basis/module_ao/element_basis_index-ORB.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_basis/module_pw/pw_basis_k.h"
#include "source_cell/unitcell.h"
#include "source_estate/elecstate.h"
#include "source_estate/module_pot/potential_new.h"
#include "source_hamilt/module_xc/exx_info.h"
#include "source_io/module_bessel/numerical_basis.h"
#include "source_io/module_parameter/parameter.h"
#include "source_lcao/module_gint/gint_atom.h"
#include "source_lcao/module_ri/Matrix_Orbs11.h"
#include "source_lcao/module_ri/conv_coulomb_pot_k.h"
#include "source_lcao/module_ri/exx_abfs-construct_orbs.h"
#include "source_lcao/module_ri/exx_abfs-io.h"
#include "source_lcao/module_ri/rpa_abfs_preorthogonalization.h"
#include "source_lcao/module_ri/sternheimer_abfs_perturbation.h"
#include "source_lcao/module_ri/sternheimer_channel_parallel.h"
#include "source_lcao/module_ri/sternheimer_channel_resources.h"
#include "source_lcao/module_ri/sternheimer_chi0_mpi.h"
#include "source_lcao/module_ri/sternheimer_coulomb_whitening.h"
#include "source_lcao/module_ri/sternheimer_delta.h"
#include "source_lcao/module_ri/sternheimer_abacus_fd_nonlocal.h"
#include "source_lcao/module_ri/sternheimer_fd_solver.h"
#include "source_lcao/module_ri/sternheimer_grid_diagnostics.h"
#include "source_lcao/module_ri/sternheimer_galerkin_audit.h"
#include "source_lcao/module_ri/sternheimer_weak_grid.h"
#include "source_lcao/module_ri/sternheimer_weak_augmented.h"
#include "source_lcao/module_ri/sternheimer_weak_matrix_audit.h"
#include "source_lcao/module_ri/sternheimer_periodic_solver.h"
#include "source_lcao/module_ri/sternheimer_rpa.h"
#include "source_lcao/module_ri/sternheimer_response_grid.h"
#include "source_lcao/module_ri/sternheimer_runtime_options.h"
#include "source_lcao/module_ri/sternheimer_siab_mpi.h"
#include "source_lcao/module_ri/sternheimer_siab_overlap.h"
#include "source_lcao/module_ri/sternheimer_siab_provenance.h"
#include "source_lcao/module_ri/sternheimer_siab_writer.h"
#include "source_lcao/module_ri/sternheimer_supercell_perturbation.h"
#include "source_lcao/module_ri/sternheimer_supercell_sector.h"
#include "source_lcao/module_ri/sternheimer_wavefunction_diagnostic.h"
#include "source_pw/module_pwdft/structure_factor.h"
#include "source_base/module_external/blas_connector.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <memory>
#include <map>
#include <memory>
#include <limits>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

#ifdef COMMIT_INFO
#include "commit.h"
#endif

#ifdef __MPI
#include <mpi.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

namespace ModuleRI
{
namespace
{

namespace siab = ::module_ri::sternheimer_siab;
namespace sternheimer_chi0 = ::module_ri::sternheimer_chi0;

std::string sternheimer_abfs_perturbation_source(const std::vector<std::string>& explicit_abfs_files)
{
    return explicit_abfs_files.empty() ? "product_pca" : "explicit_abfs";
}

bool sternheimer_builds_product_pca_auxiliary_basis(const std::vector<std::string>& explicit_abfs_files)
{
    return explicit_abfs_files.empty();
}

constexpr const char* kSmokeEnv = "ABACUS_STERNHEIMER_FD_ST_SMOKE";
constexpr const char* kOutputEnv = "ABACUS_STERNHEIMER_FD_ST_OUT";
constexpr const char* kBandsEnv = "ABACUS_STERNHEIMER_FD_ST_BANDS";
constexpr const char* kChannelsEnv = "ABACUS_STERNHEIMER_FD_ST_CHANNELS";
constexpr const char* kChannelThreadsEnv = "ABACUS_STERNHEIMER_CHANNEL_THREADS";
constexpr const char* kChannelMaxWorkersEnv = "ABACUS_STERNHEIMER_CHANNEL_MAX_WORKERS";
constexpr const char* kMaxDenseEnv = "ABACUS_STERNHEIMER_FD_ST_MAX_DENSE";
constexpr const char* kLanczosSubspaceEnv = "ABACUS_STERNHEIMER_FD_ST_LANCZOS_SUBSPACE";
constexpr const char* kOmegaEnv = "ABACUS_STERNHEIMER_FD_ST_OMEGA";
constexpr const char* kSolverToleranceEnv = "ABACUS_STERNHEIMER_FD_ST_SOLVER_TOL";
constexpr const char* kSolverMaxIterEnv = "ABACUS_STERNHEIMER_FD_ST_MAX_ITER";
constexpr const char* kSpectralPreconditionerEnv = "ABACUS_STERNHEIMER_FD_ST_SPECTRAL_PRECONDITIONER";
constexpr const char* kSpectralPreconditionerRegularizationEnv
    = "ABACUS_STERNHEIMER_FD_ST_SPECTRAL_PRECONDITIONER_REGULARIZATION";
constexpr const char* kPCAThresholdEnv = "ABACUS_STERNHEIMER_FD_ST_PCA_THRESHOLD";
constexpr const char* kCCPRmeshTimesEnv = "ABACUS_STERNHEIMER_FD_ST_CCP_RMESH_TIMES";
constexpr const char* kOrbitalDirEnv = "ABACUS_STERNHEIMER_FD_ST_ORBITAL_DIR";
constexpr const char* kOrbitalFilesEnv = "ABACUS_STERNHEIMER_FD_ST_ORBITAL_FILES";
constexpr const char* kFrequencyRankShiftEnv = "ABACUS_STERNHEIMER_FD_ST_FREQ_RANK_SHIFT";
constexpr const char* kDeltaComponentDiagnosticEnv = "ABACUS_STERNHEIMER_DELTA_COMPONENT_DIAG";
constexpr const char* kLCAOSOSDiagnosticEnv = "ABACUS_STERNHEIMER_LCAO_SOS_DIAG";
constexpr const char* kKResolvedDiagnosticEnv = "ABACUS_STERNHEIMER_KRESOLVED_DIAG";
constexpr const char* kWavefunctionDiagnosticEnv = "ABACUS_STERNHEIMER_WAVEFUNCTION_DIAGNOSTIC";
constexpr const char* kSupercellTranslationSumEnv = "ABACUS_STERNHEIMER_SUPERCELL_TRANSLATION_SUM";
constexpr const char* kSupercellFullResponseEnv = "ABACUS_STERNHEIMER_SUPERCELL_FULL_RESPONSE";
constexpr const char* kSupercellKPointGroupsEnv = "ABACUS_STERNHEIMER_SUPERCELL_KPOINT_GROUPS";
constexpr const char* kTruncatedV1DiagnosticEnv = "ABACUS_STERNHEIMER_TRUNCATED_V1_DIAGNOSTIC";
constexpr const char* kDeltaABlockModeEnv = "ABACUS_STERNHEIMER_DELTA_A_BLOCK";
constexpr double kHartreeToRydberg = 2.0;

void hash_u64(siab::Sha256& digest, const std::uint64_t value)
{
    std::array<unsigned char, 8> bytes;
    for (std::size_t index = 0; index != bytes.size(); ++index)
    {
        bytes[index] = static_cast<unsigned char>((value >> (56U - 8U * index)) & 0xffU);
    }
    digest.update(bytes.data(), bytes.size());
}

void hash_int(siab::Sha256& digest, const int value)
{
    hash_u64(digest, static_cast<std::uint64_t>(static_cast<std::int64_t>(value)));
}

void hash_double(siab::Sha256& digest, const double value)
{
    static_assert(sizeof(double) == sizeof(std::uint64_t), "Sternheimer SIAB provenance requires binary64 doubles.");
    std::uint64_t bits = 0;
    std::memcpy(&bits, &value, sizeof(bits));
    hash_u64(digest, bits);
}

std::string hash_coulomb_whitening_transform(const SternheimerCoulombWhitening& whitening)
{
    siab::Sha256 digest;
    const std::string marker = "ABACUS_STERNHEIMER_COULOMB_WHITENING_V1";
    digest.update(reinterpret_cast<const unsigned char*>(marker.data()), marker.size());
    hash_int(digest, whitening.raw_dimension);
    hash_int(digest, whitening.retained_rank);
    hash_double(digest, whitening.relative_threshold);
    for (const double value: whitening.transform)
    {
        if (!std::isfinite(value))
        {
            throw std::runtime_error("Sternheimer Coulomb whitening transform contains a non-finite value.");
        }
        hash_double(digest, value);
    }
    return digest.finish();
}

std::string compiled_commit_metadata()
{
#ifdef COMMIT_INFO
    return COMMIT;
#else
    return std::string();
#endif
}

std::vector<double> cell_vectors_bohr(const UnitCell& ucell)
{
    return {ucell.lat0 * ucell.latvec.e11,
            ucell.lat0 * ucell.latvec.e12,
            ucell.lat0 * ucell.latvec.e13,
            ucell.lat0 * ucell.latvec.e21,
            ucell.lat0 * ucell.latvec.e22,
            ucell.lat0 * ucell.latvec.e23,
            ucell.lat0 * ucell.latvec.e31,
            ucell.lat0 * ucell.latvec.e32,
            ucell.lat0 * ucell.latvec.e33};
}

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

int sternheimer_channel_openmp_threads()
{
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
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

std::uint64_t slurm_memory_limit_bytes()
{
    const char* raw = std::getenv("SLURM_MEM_PER_NODE");
    if (raw == nullptr || raw[0] == '\0')
    {
        return 0;
    }
    char* end = nullptr;
    const unsigned long long memory_mib = std::strtoull(raw, &end, 10);
    if (end == raw || *end != '\0' || memory_mib == 0
        || memory_mib > std::numeric_limits<std::uint64_t>::max() / (1024U * 1024U))
    {
        return 0;
    }
    return static_cast<std::uint64_t>(memory_mib) * 1024U * 1024U;
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
    return join_output_path(output_dir,
                            PARAM.inp.out_sternheimer_siab ? "STERNHEIMER_SIAB_STATUS.dat"
                                                          : "STERNHEIMER_CHI0.dat");
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
    out << (PARAM.inp.out_sternheimer_siab ? "STERNHEIMER_SIAB_PROGRESS_rank"
                                          : "STERNHEIMER_CHI0_PROGRESS_rank")
        << rank << ".dat";
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

std::string precise_double_string(const double value)
{
    std::ostringstream out;
    out << std::setprecision(16) << value;
    return out.str();
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
    static std::mutex progress_mutex;
    const std::lock_guard<std::mutex> lock(progress_mutex);
    std::ofstream out(chi0_progress_filename().c_str(), std::ios::app);
    out << std::setprecision(16);
    out << event << ' ' << GlobalV::MY_RANK << ' ' << ifrequency << ' ' << owner_rank << ' ' << band << ' ' << channel
        << ' ' << solved_equations << ' ';
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
                eigenvalues_ry_from_elec_state(elec_state, record.zero_order_k_index),
                occupations_from_elec_state(elec_state, record.zero_order_k_index));
        combined.emin_ha = std::min(combined.emin_ha, window.emin_ha);
        combined.emax_ha = std::max(combined.emax_ha, window.emax_ha);
    }
    return combined;
}

SternheimerRPA::TransitionEnergyWindow transition_window_from_recovered_kpoints(
    const std::vector<SternheimerLCAOOccupiedKPoint>& records)
{
    SternheimerRPA::TransitionEnergyWindow combined;
    combined.emin_ha = std::numeric_limits<double>::max();
    combined.emax_ha = 0.0;
    for (const SternheimerLCAOOccupiedKPoint& record: records)
    {
        std::vector<double> eigenvalues = record.eigenvalues;
        eigenvalues.insert(eigenvalues.end(),
                           record.unoccupied_eigenvalues.begin(),
                           record.unoccupied_eigenvalues.end());
        std::vector<double> occupations = record.occupations;
        occupations.resize(eigenvalues.size(), 0.0);
        const auto window = SternheimerRPA::transition_energy_window_from_eigenvalues_ry(
            eigenvalues, occupations);
        combined.emin_ha = std::min(combined.emin_ha, window.emin_ha);
        combined.emax_ha = std::max(combined.emax_ha, window.emax_ha);
    }
    return combined;
}

struct SternheimerFixedQPermutation
{
    int spatial_isym = -1;
    bool time_reversal = false;
    std::vector<int> mapped_index_by_full_k;
    std::vector<std::array<int, 3>> fold_G_by_full_k;
};

std::vector<SternheimerFixedQPermutation> build_sternheimer_fixed_q_little_group_permutations(
    const UnitCell& ucell,
    const std::vector<SternheimerLCAOOccupiedKPoint>& records,
    const SternheimerReducedKPoint& qpoint,
    const std::vector<int>& allowed_spatial_operations)
{
    if (ucell.symm.nrotk <= 0 || records.empty() || allowed_spatial_operations.empty())
    {
        throw std::runtime_error("Cannot build Sternheimer fixed-q permutations without symmetry operations.");
    }
    std::vector<bool> spatial_operation_allowed(static_cast<std::size_t>(ucell.symm.nrotk), false);
    for (const int isym: allowed_spatial_operations)
    {
        if (isym < 0 || isym >= ucell.symm.nrotk)
        {
            throw std::runtime_error("Sternheimer FD stencil symmetry has an invalid spatial index.");
        }
        spatial_operation_allowed[static_cast<std::size_t>(isym)] = true;
    }
    const double tolerance = std::max(1.0e-10, ucell.symm.epsilon);
    auto canonical = [tolerance](ModuleBase::Vector3<double> point) {
        point.x = std::fmod(point.x + 100.5 - 0.5 * tolerance, 1.0) - 0.5 + 0.5 * tolerance;
        point.y = std::fmod(point.y + 100.5 - 0.5 * tolerance, 1.0) - 0.5 + 0.5 * tolerance;
        point.z = std::fmod(point.z + 100.5 - 0.5 * tolerance, 1.0) - 0.5 + 0.5 * tolerance;
        if (std::abs(point.x) < tolerance)
        {
            point.x = 0.0;
        }
        if (std::abs(point.y) < tolerance)
        {
            point.y = 0.0;
        }
        if (std::abs(point.z) < tolerance)
        {
            point.z = 0.0;
        }
        return point;
    };
    auto same_point = [tolerance](const ModuleBase::Vector3<double>& lhs,
                                  const ModuleBase::Vector3<double>& rhs) {
        return std::abs(lhs.x - rhs.x) < tolerance && std::abs(lhs.y - rhs.y) < tolerance
               && std::abs(lhs.z - rhs.z) < tolerance;
    };

    std::vector<ModuleBase::Vector3<double>> full_kpoints(records.size());
    for (const auto& record: records)
    {
        if (record.global_k_index < 0 || record.global_k_index >= static_cast<int>(records.size()))
        {
            throw std::runtime_error("Sternheimer fixed-q permutation found an invalid full-k index.");
        }
        full_kpoints[static_cast<std::size_t>(record.global_k_index)]
            = canonical({record.kpoint[0], record.kpoint[1], record.kpoint[2]});
    }
    const ModuleBase::Vector3<double> q = canonical({qpoint[0], qpoint[1], qpoint[2]});
    std::vector<SternheimerFixedQPermutation> permutations;
    for (const bool time_reversal: {false, true})
    {
        for (int isym = 0; isym != ucell.symm.nrotk; ++isym)
        {
            if (!spatial_operation_allowed[static_cast<std::size_t>(isym)])
            {
                continue;
            }
            ModuleBase::Vector3<double> transformed_q = q * ucell.symm.kgmatrix[isym];
            if (time_reversal)
            {
                transformed_q = transformed_q * -1.0;
            }
            if (!same_point(canonical(transformed_q), q))
            {
                continue;
            }

            SternheimerFixedQPermutation operation;
            operation.spatial_isym = isym;
            operation.time_reversal = time_reversal;
            operation.mapped_index_by_full_k.assign(records.size(), -1);
            operation.fold_G_by_full_k.assign(records.size(), {0, 0, 0});
            for (std::size_t ik = 0; ik != full_kpoints.size(); ++ik)
            {
                ModuleBase::Vector3<double> transformed
                    = full_kpoints[ik] * ucell.symm.kgmatrix[isym];
                if (time_reversal)
                {
                    transformed = transformed * -1.0;
                }
                const ModuleBase::Vector3<double> unfolded = transformed;
                transformed = canonical(unfolded);
                for (std::size_t target = 0; target != full_kpoints.size(); ++target)
                {
                    if (same_point(transformed, full_kpoints[target]))
                    {
                        operation.mapped_index_by_full_k[ik] = static_cast<int>(target);
                        const ModuleBase::Vector3<double> difference
                            = unfolded - full_kpoints[target];
                        operation.fold_G_by_full_k[ik]
                            = {static_cast<int>(std::llround(difference.x)),
                               static_cast<int>(std::llround(difference.y)),
                               static_cast<int>(std::llround(difference.z))};
                        if (std::abs(difference.x - operation.fold_G_by_full_k[ik][0]) >= tolerance
                            || std::abs(difference.y - operation.fold_G_by_full_k[ik][1]) >= tolerance
                            || std::abs(difference.z - operation.fold_G_by_full_k[ik][2]) >= tolerance)
                        {
                            throw std::runtime_error(
                                "A Sternheimer fixed-q symmetry route has a noninteger reciprocal fold.");
                        }
                        break;
                    }
                }
                if (operation.mapped_index_by_full_k[ik] < 0)
                {
                    throw std::runtime_error(
                        "A Sternheimer fixed-q little-group operation leaves the full k grid.");
                }
            }
            permutations.push_back(std::move(operation));
        }
    }
    if (permutations.empty())
    {
        throw std::runtime_error("Sternheimer fixed-q little group is empty.");
    }
    return permutations;
}

std::vector<SternheimerQStarPermutation> build_sternheimer_discrete_qstar_permutations(
    const UnitCell& ucell,
    const std::vector<SternheimerLCAOOccupiedKPoint>& records,
    const std::vector<int>& allowed_spatial_operations)
{
    if (ucell.symm.nrotk <= 0 || records.empty() || allowed_spatial_operations.empty())
    {
        throw std::runtime_error("Cannot build discrete Sternheimer q-star permutations.");
    }
    std::vector<bool> spatial_operation_allowed(static_cast<std::size_t>(ucell.symm.nrotk), false);
    for (const int isym: allowed_spatial_operations)
    {
        if (isym < 0 || isym >= ucell.symm.nrotk)
        {
            throw std::runtime_error("Sternheimer q-star symmetry has an invalid spatial index.");
        }
        spatial_operation_allowed[static_cast<std::size_t>(isym)] = true;
    }

    const double tolerance = std::max(1.0e-10, ucell.symm.epsilon);
    auto canonical = [tolerance](ModuleBase::Vector3<double> point) {
        point.x = std::fmod(point.x + 100.5 - 0.5 * tolerance, 1.0) - 0.5 + 0.5 * tolerance;
        point.y = std::fmod(point.y + 100.5 - 0.5 * tolerance, 1.0) - 0.5 + 0.5 * tolerance;
        point.z = std::fmod(point.z + 100.5 - 0.5 * tolerance, 1.0) - 0.5 + 0.5 * tolerance;
        if (std::abs(point.x) < tolerance)
        {
            point.x = 0.0;
        }
        if (std::abs(point.y) < tolerance)
        {
            point.y = 0.0;
        }
        if (std::abs(point.z) < tolerance)
        {
            point.z = 0.0;
        }
        return point;
    };
    auto same_point = [tolerance](const ModuleBase::Vector3<double>& lhs,
                                  const ModuleBase::Vector3<double>& rhs) {
        return std::abs(lhs.x - rhs.x) < tolerance && std::abs(lhs.y - rhs.y) < tolerance
               && std::abs(lhs.z - rhs.z) < tolerance;
    };

    std::vector<ModuleBase::Vector3<double>> full_qpoints(records.size());
    for (const auto& record: records)
    {
        if (record.global_k_index < 0 || record.global_k_index >= static_cast<int>(records.size()))
        {
            throw std::runtime_error("Sternheimer q-star permutation found an invalid full-q index.");
        }
        full_qpoints[static_cast<std::size_t>(record.global_k_index)]
            = canonical({record.kpoint[0], record.kpoint[1], record.kpoint[2]});
    }

    std::vector<SternheimerQStarPermutation> permutations;
    for (const bool time_reversal: {false, true})
    {
        for (int isym = 0; isym != ucell.symm.nrotk; ++isym)
        {
            if (!spatial_operation_allowed[static_cast<std::size_t>(isym)])
            {
                continue;
            }
            SternheimerQStarPermutation operation;
            operation.spatial_isym = isym;
            operation.time_reversal = time_reversal;
            operation.mapped_index_by_full_q.assign(records.size(), -1);
            operation.fold_G_by_full_q.assign(records.size(), {0, 0, 0});
            for (std::size_t iq = 0; iq != full_qpoints.size(); ++iq)
            {
                ModuleBase::Vector3<double> transformed
                    = full_qpoints[iq] * ucell.symm.kgmatrix[isym];
                if (time_reversal)
                {
                    transformed = transformed * -1.0;
                }
                const ModuleBase::Vector3<double> unfolded = transformed;
                transformed = canonical(unfolded);
                for (std::size_t target = 0; target != full_qpoints.size(); ++target)
                {
                    if (!same_point(transformed, full_qpoints[target]))
                    {
                        continue;
                    }
                    operation.mapped_index_by_full_q[iq] = static_cast<int>(target);
                    const ModuleBase::Vector3<double> difference = unfolded - full_qpoints[target];
                    operation.fold_G_by_full_q[iq]
                        = {static_cast<int>(std::llround(difference.x)),
                           static_cast<int>(std::llround(difference.y)),
                           static_cast<int>(std::llround(difference.z))};
                    const auto& fold = operation.fold_G_by_full_q[iq];
                    if (std::abs(difference.x - fold[0]) >= tolerance
                        || std::abs(difference.y - fold[1]) >= tolerance
                        || std::abs(difference.z - fold[2]) >= tolerance)
                    {
                        throw std::runtime_error(
                            "A Sternheimer q-star route has a noninteger reciprocal fold.");
                    }
                    break;
                }
                if (operation.mapped_index_by_full_q[iq] < 0)
                {
                    throw std::runtime_error(
                        "A discrete Sternheimer q-star operation leaves the full q grid.");
                }
            }
            permutations.push_back(std::move(operation));
        }
    }
    if (permutations.empty())
    {
        throw std::runtime_error("The discrete Sternheimer q-star group is empty.");
    }
    return permutations;
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
        throw std::runtime_error("No numerical orbital files are available. Set NUMERICAL_ORBITAL in STRU or "
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

siab::Provenance make_siab_production_provenance(const UnitCell& ucell,
                                                 const std::string& auxiliary_basis_sha256,
                                                 const SternheimerRPA::FrequencyGrid& frequency_grid,
                                                 const double pca_threshold,
                                                 const SternheimerCoulombWhitening& whitening)
{
    siab::Provenance provenance;
    provenance.abacus_commit = siab::require_source_commit(compiled_commit_metadata());
    provenance.auxiliary_basis_sha256 = auxiliary_basis_sha256;
    provenance.cell_bohr = cell_vectors_bohr(ucell);
    provenance.ecut_ry = PARAM.inp.ecutwfc;
    provenance.kernel = "full_coulomb";
    const std::vector<std::string> orbital_files
        = siab::resolve_required_input_files(orbital_dir_from_env_or_input(),
                                             orbital_files_from_env_or_cell(ucell),
                                             "initial orbital");
    const std::vector<std::string> pseudopotential_files
        = siab::resolve_required_input_files(PARAM.inp.pseudo_dir, ucell.pseudo_fn, "pseudopotential");
    provenance.orbital_sha256 = siab::sha256_file_manifest(orbital_files);
    provenance.pseudopotential_sha256 = siab::sha256_file_manifest(pseudopotential_files);
    provenance.spin_convention = "spin_resolved_occupation_in_reference_rows";
    provenance.executable_sha256 = siab::sha256_file(siab::resolve_executable_path());
    provenance.exx_pca_thr = pca_threshold;
    provenance.sternheimer_nfreq = static_cast<int>(frequency_grid.omega_ha.size());
    provenance.frequency_ha = frequency_grid.omega_ha;
    provenance.frequency_weights_ha = frequency_grid.weights_ha;
    provenance.mpi_ranks = GlobalV::NPROC;
    provenance.omp_threads = PARAM.globalv.nthread_per_proc;
    provenance.auxiliary_whitening = "global_full_coulomb_v1";
    provenance.raw_auxiliary_dimension = whitening.raw_dimension;
    provenance.whitened_auxiliary_rank = whitening.retained_rank;
    provenance.discarded_auxiliary_rank = whitening.discarded_rank;
    provenance.coulomb_relative_threshold = whitening.relative_threshold;
    provenance.coulomb_eigenvalues = whitening.eigenvalues;
    provenance.coulomb_max_orthonormality_error = whitening.max_orthonormality_error;
    provenance.coulomb_transform_sha256 = hash_coulomb_whitening_transform(whitening);
    return provenance;
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

struct SternheimerABFSInput
{
    std::vector<SternheimerABFGridChannel> channels;
    std::vector<double> full_coulomb_metric;
    std::vector<std::vector<SternheimerRadialPerturbation>> radials_by_type;
    std::vector<int> atom_types;
    std::vector<ModuleBase::Vector3<double>> atom_positions;
};

std::vector<double> build_molecular_coulomb_metric(
    const UnitCell& ucell,
    const LCAO_Orbitals& orb,
    const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>& abfs_ccp,
    const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>& abfs,
    const std::vector<SternheimerABFGridChannel>& channels)
{
    const ModuleBase::Element_Basis_Index::Range range_ccp
        = ModuleBase::Element_Basis_Index::construct_range(abfs_ccp);
    const ModuleBase::Element_Basis_Index::Range range_abfs
        = ModuleBase::Element_Basis_Index::construct_range(abfs);
    const ModuleBase::Element_Basis_Index::IndexLNM index_ccp
        = ModuleBase::Element_Basis_Index::construct_index(range_ccp);
    const ModuleBase::Element_Basis_Index::IndexLNM index_abfs
        = ModuleBase::Element_Basis_Index::construct_index(range_abfs);
    if (index_ccp.size() != index_abfs.size())
    {
        throw std::runtime_error("Sternheimer Coulomb metric ABFS type counts do not match.");
    }
    for (std::size_t type = 0; type != index_ccp.size(); ++type)
    {
        if (index_ccp[type].count_size != index_abfs[type].count_size)
        {
            throw std::runtime_error("Sternheimer Coulomb metric ABFS channel counts do not match.");
        }
    }

    std::vector<std::size_t> atom_offsets;
    std::vector<int> atom_types;
    std::vector<int> atom_local_indices;
    atom_offsets.reserve(static_cast<std::size_t>(ucell.nat));
    atom_types.reserve(static_cast<std::size_t>(ucell.nat));
    atom_local_indices.reserve(static_cast<std::size_t>(ucell.nat));
    std::size_t total_dimension = 0;
    for (int type = 0; type != ucell.ntype; ++type)
    {
        for (int atom = 0; atom != ucell.atoms[type].na; ++atom)
        {
            atom_offsets.push_back(total_dimension);
            atom_types.push_back(type);
            atom_local_indices.push_back(atom);
            total_dimension += index_abfs[static_cast<std::size_t>(type)].count_size;
        }
    }
    if (total_dimension != channels.size())
    {
        throw std::runtime_error("Sternheimer Coulomb metric dimension does not match sampled ABFS channels.");
    }
    for (const SternheimerABFGridChannel& channel: channels)
    {
        const std::size_t atom = static_cast<std::size_t>(channel.atom_index);
        if (atom >= atom_offsets.size()
            || atom_offsets[atom] + static_cast<std::size_t>(channel.atom_local_index)
                   != static_cast<std::size_t>(channel.channel_index))
        {
            throw std::runtime_error("Sternheimer Coulomb metric channel ordering is inconsistent with ABFS indices.");
        }
    }

    Matrix_Orbs11 metric_builder;
    metric_builder.init(abfs_ccp, abfs, ucell, orb, GlobalC::exx_info.info_ri.kmesh_times);
    metric_builder.init_radial_table();
    std::vector<double> metric(total_dimension * total_dimension, 0.0);
    for (std::size_t atom_a = 0; atom_a != atom_offsets.size(); ++atom_a)
    {
        const std::size_t type_a = static_cast<std::size_t>(atom_types[atom_a]);
        const auto& tau_a = ucell.atoms[type_a].tau[static_cast<std::size_t>(atom_local_indices[atom_a])];
        const std::size_t count_a = index_ccp[type_a].count_size;
        for (std::size_t atom_b = 0; atom_b != atom_offsets.size(); ++atom_b)
        {
            const std::size_t type_b = static_cast<std::size_t>(atom_types[atom_b]);
            const auto& tau_b = ucell.atoms[type_b].tau[static_cast<std::size_t>(atom_local_indices[atom_b])];
            const std::size_t count_b = index_abfs[type_b].count_size;
            const RI::Tensor<double> block
                = metric_builder.cal_overlap_matrix<double>(type_a,
                                                             type_b,
                                                             tau_a,
                                                             tau_b,
                                                             index_ccp,
                                                             index_abfs,
                                                             Matrix_Orbs11::Matrix_Order::AB);
            for (std::size_t local_a = 0; local_a != count_a; ++local_a)
            {
                for (std::size_t local_b = 0; local_b != count_b; ++local_b)
                {
                    metric[(atom_offsets[atom_a] + local_a) * total_dimension
                           + atom_offsets[atom_b] + local_b]
                        = block(local_a, local_b);
                }
            }
        }
    }
    return metric;
}

SternheimerABFBuildData build_abfs_ccp_data(const UnitCell& ucell,
                                             const SternheimerFDHamiltonian::Grid& grid,
                                             const int max_channels,
                                             const double pca_threshold,
                                             const double ccp_rmesh_times,
                                             const bool build_coulomb_metric)
{
    LCAO_Orbitals orb;
    read_sternheimer_orbitals(ucell, orb);
    std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> abfs;
    if (sternheimer_builds_product_pca_auxiliary_basis(GlobalC::exx_info.info_ri.files_abfs))
    {
        auto lcaos = Exx_Abfs::Construct_Orbs::change_orbs(orb, GlobalC::exx_info.info_ri.kmesh_times);
        Exx_Abfs::Construct_Orbs::filter_empty_orbs(lcaos);
        abfs = Exx_Abfs::Construct_Orbs::abfs_same_atom(ucell,
                                                        orb,
                                                        lcaos,
                                                        GlobalC::exx_info.info_ri.kmesh_times,
                                                        pca_threshold);
    }
    else
    {
        abfs = Exx_Abfs::IO::construct_abfs(orb,
                                            GlobalC::exx_info.info_ri.files_abfs,
                                            GlobalC::exx_info.info_ri.kmesh_times);
    }
    Exx_Abfs::Construct_Orbs::filter_empty_orbs(abfs);
    const RpaAbfsPreorthReport preorth_report
        = finalize_rpa_abfs_from_input(abfs, PARAM.inp, PARAM.inp.cal_force);
    if (GlobalV::MY_RANK == 0)
    {
        GlobalV::ofs_running << format_rpa_abfs_preorth_report(preorth_report);
    }
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

SternheimerABFSInput build_abfs_density_input(const UnitCell& ucell,
                                              const double pca_threshold,
                                              const SternheimerOrbitalSet* reusable_rpa_abfs = nullptr)
{
    if (reusable_rpa_abfs != nullptr)
    {
        return make_sternheimer_abfs_input(ucell, *reusable_rpa_abfs);
    }
    const SternheimerOrbitalSet abfs = build_sternheimer_abfs(ucell, pca_threshold);
    return make_sternheimer_abfs_input(ucell, abfs);
}

struct SternheimerABFBuildData
{
    std::vector<SternheimerABFGridChannel> channels;
    std::vector<double> full_coulomb_metric;
    std::vector<std::vector<SternheimerRadialPerturbation>> radials_by_type;
    std::vector<int> atom_types;
    std::vector<ModuleBase::Vector3<double>> atom_positions;
    RpaAbfsPreorthReport preorth_report;
};

std::vector<double> build_molecular_coulomb_metric(
    const UnitCell& ucell,
    const LCAO_Orbitals& orb,
    const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>& abfs_ccp,
    const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>& abfs,
    const std::vector<SternheimerABFGridChannel>& channels)
{
    const ModuleBase::Element_Basis_Index::Range range_ccp
        = ModuleBase::Element_Basis_Index::construct_range(abfs_ccp);
    const ModuleBase::Element_Basis_Index::Range range_abfs
        = ModuleBase::Element_Basis_Index::construct_range(abfs);
    const ModuleBase::Element_Basis_Index::IndexLNM index_ccp
        = ModuleBase::Element_Basis_Index::construct_index(range_ccp);
    const ModuleBase::Element_Basis_Index::IndexLNM index_abfs
        = ModuleBase::Element_Basis_Index::construct_index(range_abfs);
    if (index_ccp.size() != index_abfs.size())
    {
        throw std::runtime_error("Sternheimer Coulomb metric ABFS type counts do not match.");
    }
    for (std::size_t type = 0; type != index_ccp.size(); ++type)
    {
        if (index_ccp[type].count_size != index_abfs[type].count_size)
        {
            throw std::runtime_error("Sternheimer Coulomb metric ABFS channel counts do not match.");
        }
    }

    std::vector<std::size_t> atom_offsets;
    std::vector<int> atom_types;
    std::vector<int> atom_local_indices;
    atom_offsets.reserve(static_cast<std::size_t>(ucell.nat));
    atom_types.reserve(static_cast<std::size_t>(ucell.nat));
    atom_local_indices.reserve(static_cast<std::size_t>(ucell.nat));
    std::size_t total_dimension = 0;
    for (int type = 0; type != ucell.ntype; ++type)
    {
        for (int atom = 0; atom != ucell.atoms[type].na; ++atom)
        {
            atom_offsets.push_back(total_dimension);
            atom_types.push_back(type);
            atom_local_indices.push_back(atom);
            total_dimension += index_abfs[static_cast<std::size_t>(type)].count_size;
        }
    }
    if (total_dimension != channels.size())
    {
        throw std::runtime_error("Sternheimer Coulomb metric dimension does not match sampled ABFS channels.");
    }
    for (const SternheimerABFGridChannel& channel: channels)
    {
        const std::size_t atom = static_cast<std::size_t>(channel.atom_index);
        if (atom >= atom_offsets.size()
            || atom_offsets[atom] + static_cast<std::size_t>(channel.atom_local_index)
                   != static_cast<std::size_t>(channel.channel_index))
        {
            throw std::runtime_error("Sternheimer Coulomb metric channel ordering is inconsistent with ABFS indices.");
        }
    }

    Matrix_Orbs11 metric_builder;
    metric_builder.init(abfs_ccp, abfs, ucell, orb, GlobalC::exx_info.info_ri.kmesh_times);
    metric_builder.init_radial_table();
    std::vector<double> metric(total_dimension * total_dimension, 0.0);
    for (std::size_t atom_a = 0; atom_a != atom_offsets.size(); ++atom_a)
    {
        const std::size_t type_a = static_cast<std::size_t>(atom_types[atom_a]);
        const auto& tau_a = ucell.atoms[type_a].tau[static_cast<std::size_t>(atom_local_indices[atom_a])];
        const std::size_t count_a = index_ccp[type_a].count_size;
        for (std::size_t atom_b = 0; atom_b != atom_offsets.size(); ++atom_b)
        {
            const std::size_t type_b = static_cast<std::size_t>(atom_types[atom_b]);
            const auto& tau_b = ucell.atoms[type_b].tau[static_cast<std::size_t>(atom_local_indices[atom_b])];
            const std::size_t count_b = index_abfs[type_b].count_size;
            const RI::Tensor<double> block
                = metric_builder.cal_overlap_matrix<double>(type_a,
                                                             type_b,
                                                             tau_a,
                                                             tau_b,
                                                             index_ccp,
                                                             index_abfs,
                                                             Matrix_Orbs11::Matrix_Order::AB);
            for (std::size_t local_a = 0; local_a != count_a; ++local_a)
            {
                for (std::size_t local_b = 0; local_b != count_b; ++local_b)
                {
                    metric[(atom_offsets[atom_a] + local_a) * total_dimension
                           + atom_offsets[atom_b] + local_b]
                        = block(local_a, local_b);
                }
            }
        }
    }
    return metric;
}

SternheimerABFBuildData build_abfs_ccp_data(const UnitCell& ucell,
                                             const SternheimerFDHamiltonian::Grid& grid,
                                             const int max_channels,
                                             const double pca_threshold,
                                             const double ccp_rmesh_times,
                                             const bool build_coulomb_metric)
{
    LCAO_Orbitals orb;
    read_sternheimer_orbitals(ucell, orb);
    std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> abfs;
    if (sternheimer_builds_product_pca_auxiliary_basis(GlobalC::exx_info.info_ri.files_abfs))
    {
        auto lcaos = Exx_Abfs::Construct_Orbs::change_orbs(orb, GlobalC::exx_info.info_ri.kmesh_times);
        Exx_Abfs::Construct_Orbs::filter_empty_orbs(lcaos);
        abfs = Exx_Abfs::Construct_Orbs::abfs_same_atom(ucell,
                                                        orb,
                                                        lcaos,
                                                        GlobalC::exx_info.info_ri.kmesh_times,
                                                        pca_threshold);
    }
    else
    {
        abfs = Exx_Abfs::IO::construct_abfs(orb,
                                            GlobalC::exx_info.info_ri.files_abfs,
                                            GlobalC::exx_info.info_ri.kmesh_times);
    }
    Exx_Abfs::Construct_Orbs::filter_empty_orbs(abfs);
    const RpaAbfsPreorthReport preorth_report
        = finalize_rpa_abfs_from_input(abfs, PARAM.inp, PARAM.inp.cal_force);
    if (GlobalV::MY_RANK == 0)
    {
        GlobalV::ofs_running << format_rpa_abfs_preorth_report(preorth_report);
    }
    const auto abfs_ccp = Conv_Coulomb_Pot_K::cal_orbs_ccp(abfs, make_fock_hartree_coulomb_param(), ccp_rmesh_times);
    auto radials_by_type = make_sternheimer_radial_perturbations_from_orbitals(abfs_ccp);

    SternheimerFDHamiltonian::Grid potential_grid = grid;
    if (PARAM.inp.sternheimer_molecular_coulomb == "isolated_ri")
    {
        if (build_coulomb_metric)
        {
            throw std::invalid_argument("isolated_ri cannot be combined with SIAB metric construction.");
        }
        // Only the perturbation is isolated; the FD Hamiltonian keeps its existing boundary condition.
        potential_grid.periodic = false;
        std::ofstream audit;
        if (GlobalV::MY_RANK == 0)
        {
            audit.open("STERNHEIMER_COULOMB_TAIL.dat");
            if (!audit)
            {
                throw std::runtime_error("Cannot write the molecular Coulomb-tail audit.");
            }
            audit << "# type l radial_index radius_bohr endpoint_Ha charge_moment_coefficient "
                     "endpoint_extrapolated_coefficient\n"
                  << std::setprecision(17);
        }
        for (auto& type_radials: radials_by_type)
        {
            for (auto& radial: type_radials)
            {
                const auto& charges = abfs.at(radial.type_index).at(radial.angular_momentum);
                const auto charge = std::find_if(charges.begin(), charges.end(), [&](const Numerical_Orbital_Lm& orb) {
                    return orb.getChi() == radial.radial_index;
                });
                if (charge == charges.end())
                {
                    throw std::runtime_error("Cannot match the Coulomb potential to its finalized auxiliary charge.");
                }
                set_sternheimer_coulomb_tail_from_charge(radial, *charge);
                if (GlobalV::MY_RANK == 0)
                {
                    audit << radial.type_index << ' ' << radial.angular_momentum << ' ' << radial.radial_index << ' '
                          << radial.radial_grid.back() << ' ' << radial.radial_values.back() << ' '
                          << radial.coulomb_tail_coefficient << ' '
                          << radial.radial_values.back()
                                 * std::pow(radial.radial_grid.back(), radial.angular_momentum + 1)
                          << '\n';
                }
            }
        }
    }

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

    SternheimerABFBuildData result;
    result.radials_by_type = std::move(radials_by_type);
    result.atom_types = std::move(atom_types);
    result.atom_positions = std::move(atom_positions);
    result.preorth_report = preorth_report;
    result.channels = build_coulomb_metric ? describe_sternheimer_abf_grid_channels(result.radials_by_type,
                                                                                    result.atom_types,
                                                                                    result.atom_positions,
                                                                                    max_channels)
                                           : sample_sternheimer_abf_grid_channels(result.radials_by_type,
                                                                                  result.atom_types,
                                                                                  result.atom_positions,
                                                                                  potential_grid,
                                                                                  max_channels);
    if (build_coulomb_metric)
    {
        if (max_channels > 0)
        {
            throw std::invalid_argument("A full Sternheimer Coulomb metric cannot be built for truncated channels.");
        }
        result.full_coulomb_metric = build_molecular_coulomb_metric(ucell, orb, abfs_ccp, abfs, result.channels);
    }
    return result;
}

struct SIABPrimitiveExportData
{
    std::vector<siab::PrimitiveBlock> blocks;
    std::vector<std::complex<double>> reciprocal_matrix;
    std::vector<std::complex<double>> overlap_s;
    std::unique_ptr<ModulePW::PW_Basis_K> serial_pw_basis;
    int primitive_count = 0;
    int reciprocal_count = 0;
};

SIABPrimitiveExportData build_siab_primitive_export_data(const ModulePW::PW_Basis& response_pw_basis,
                                                         const Structure_Factor& structure_factor,
                                                         const UnitCell& ucell)
{
    siab::require_single_primitive_rcut(PARAM.inp.bessel_nao_rcuts);
    SIABPrimitiveExportData result;
    result.serial_pw_basis.reset(new ModulePW::PW_Basis_K("cpu", "double"));
#ifdef __MPI
    result.serial_pw_basis->initmpi(1, 0, MPI_COMM_SELF);
#endif
    result.serial_pw_basis->initgrids(ucell.lat0,
                                      ucell.latvec,
                                      response_pw_basis.nx,
                                      response_pw_basis.ny,
                                      response_pw_basis.nz);
    const ModuleBase::Vector3<double> gamma(0.0, 0.0, 0.0);
    result.serial_pw_basis->initparameters(false, PARAM.inp.ecutwfc, 1, &gamma);
    result.serial_pw_basis->fft_bundle.initfftmode(PARAM.inp.fft_mode);
    result.serial_pw_basis->setuptransform();
    result.serial_pw_basis->collect_local_pw(PARAM.inp.erf_ecut, PARAM.inp.erf_height, PARAM.inp.erf_sigma);
    if (result.serial_pw_basis->nxyz != response_pw_basis.nxyz
        || result.serial_pw_basis->nrxx != response_pw_basis.nxyz || result.serial_pw_basis->nks != 1)
    {
        throw std::runtime_error("Sternheimer SIAB serial primitive FFT basis does not cover the complete response grid.");
    }
    const Numerical_Basis::SIABPrimitiveParameters parameters
        = Numerical_Basis::siab_parameters_from_input(0, PARAM.inp.sternheimer_siab_lmax);
    Numerical_Basis numerical_basis;
    const auto reciprocal_blocks = numerical_basis.siab_primitive_reciprocal_values(
        0, result.serial_pw_basis.get(), structure_factor, ucell, parameters);
    if (reciprocal_blocks.empty())
    {
        throw std::runtime_error("Sternheimer SIAB primitive construction returned no blocks.");
    }
    for (const auto& block: reciprocal_blocks)
    {
        siab::PrimitiveBlock output_block;
        output_block.element = block.element;
        output_block.atom_index = block.atom_index;
        output_block.l = block.l;
        output_block.m = block.m;
        output_block.n_primitive = block.n_primitive;
        output_block.offset = block.offset;
        result.blocks.push_back(std::move(output_block));
        for (const auto& primitive: block.values)
        {
            if (result.reciprocal_count == 0)
            {
                result.reciprocal_count = static_cast<int>(primitive.size());
            }
            if (primitive.size() != static_cast<std::size_t>(result.reciprocal_count))
            {
                throw std::runtime_error("Sternheimer SIAB reciprocal primitive sizes are inconsistent.");
            }
            result.reciprocal_matrix.insert(result.reciprocal_matrix.end(), primitive.begin(), primitive.end());
            ++result.primitive_count;
        }
    }
    if (result.primitive_count <= 0 || result.reciprocal_count != result.serial_pw_basis->npwk[0])
    {
        throw std::runtime_error("Sternheimer SIAB reciprocal primitive matrix is empty or incomplete.");
    }
    if (GlobalV::MY_RANK == 0)
    {
        result.overlap_s.assign(static_cast<std::size_t>(result.primitive_count)
                                    * static_cast<std::size_t>(result.primitive_count), ModuleBase::ZERO);
        BlasConnector::gemm('N',
                            'C',
                            result.primitive_count,
                            result.primitive_count,
                            result.reciprocal_count,
                            std::complex<double>(1.0, 0.0),
                            result.reciprocal_matrix.data(),
                            result.reciprocal_count,
                            result.reciprocal_matrix.data(),
                            result.reciprocal_count,
                            std::complex<double>(0.0, 0.0),
                            result.overlap_s.data(),
                            result.primitive_count);
    }
    return result;
}

std::vector<std::complex<double>> project_siab_response_to_primitives(
    const std::vector<std::complex<double>>& complete_response,
    const UnitCell& ucell,
    const SIABPrimitiveExportData& primitives)
{
    if (!primitives.serial_pw_basis
        || complete_response.size() != static_cast<std::size_t>(primitives.serial_pw_basis->nrxx))
    {
        throw std::invalid_argument("Sternheimer SIAB response does not match its serial primitive FFT basis.");
    }
    std::vector<std::complex<double>> response_coefficients(
        static_cast<std::size_t>(primitives.serial_pw_basis->npwk[0]), ModuleBase::ZERO);
#ifdef _OPENMP
#pragma omp critical(sternheimer_siab_real2recip)
#endif
    {
        primitives.serial_pw_basis->real2recip(complete_response.data(), response_coefficients.data(), 0);
    }
    const double coefficient_scale = std::sqrt(ucell.omega);
    for (std::complex<double>& value: response_coefficients)
    {
        value *= coefficient_scale;
    }
    return siab::overlap_q_reciprocal_contiguous(response_coefficients,
                                                 primitives.reciprocal_matrix,
                                                 primitives.primitive_count,
                                                 primitives.reciprocal_count);
}

std::vector<SternheimerABFGridChannel> build_abfs_ccp_grid_channels(const UnitCell& ucell,
                                                                    const SternheimerFDHamiltonian::Grid& grid,
                                                                    const int max_channels,
                                                                    const double pca_threshold,
                                                                    const double ccp_rmesh_times)
{
    return build_abfs_ccp_data(ucell, grid, max_channels, pca_threshold, ccp_rmesh_times, false).channels;
}

struct SIABPrimitiveExportData
{
    std::vector<siab::PrimitiveBlock> blocks;
    std::vector<std::complex<double>> reciprocal_matrix;
    std::vector<std::complex<double>> overlap_s;
    std::unique_ptr<ModulePW::PW_Basis_K> serial_pw_basis;
    int primitive_count = 0;
    int reciprocal_count = 0;
};

SIABPrimitiveExportData build_siab_primitive_export_data(const ModulePW::PW_Basis& response_pw_basis,
                                                         const Structure_Factor& structure_factor,
                                                         const UnitCell& ucell)
{
    siab::require_single_primitive_rcut(PARAM.inp.bessel_nao_rcuts);
    SIABPrimitiveExportData result;
    result.serial_pw_basis.reset(new ModulePW::PW_Basis_K("cpu", "double"));
#ifdef __MPI
    if (!enabled || GlobalV::NPROC <= 1)
    {
        return;
    }

    int channel_count = GlobalV::MY_RANK == 0 ? static_cast<int>(channels.size()) : 0;
    MPI_Bcast(&channel_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (channel_count < 0 || grid_size <= 0)
    {
        throw std::runtime_error("Invalid periodic ABFS channel broadcast dimensions.");
    }
    if (channel_count == 0)
    {
        channels.clear();
        return;
    }
    if (GlobalV::MY_RANK != 0)
    {
        channels.resize(static_cast<std::size_t>(channel_count));
    }

    for (auto& channel: channels)
    {
        std::array<int, 7> metadata{};
        if (GlobalV::MY_RANK == 0)
        {
            metadata = {channel.channel_index,
                        channel.atom_index,
                        channel.atom_local_index,
                        channel.type_index,
                        channel.angular_momentum,
                        channel.radial_index,
                        channel.magnetic_index};
        }
        MPI_Bcast(metadata.data(), static_cast<int>(metadata.size()), MPI_INT, 0, MPI_COMM_WORLD);
        if (GlobalV::MY_RANK != 0)
        {
            channel.channel_index = metadata[0];
            channel.atom_index = metadata[1];
            channel.atom_local_index = metadata[2];
            channel.type_index = metadata[3];
            channel.angular_momentum = metadata[4];
            channel.radial_index = metadata[5];
            channel.magnetic_index = metadata[6];
        }

        int label_size = GlobalV::MY_RANK == 0 ? static_cast<int>(channel.label.size()) : 0;
        MPI_Bcast(&label_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (label_size < 0)
        {
            throw std::runtime_error("Invalid periodic ABFS channel label size broadcast.");
        }
        if (GlobalV::MY_RANK != 0)
        {
            channel.label.assign(static_cast<std::size_t>(label_size), '\0');
        }
        if (label_size > 0)
        {
            MPI_Bcast(&channel.label[0], label_size, MPI_CHAR, 0, MPI_COMM_WORLD);
        }

        MPI_Bcast(&channel.max_abs, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        int value_count = GlobalV::MY_RANK == 0 ? static_cast<int>(channel.potential_r.size()) : 0;
        MPI_Bcast(&value_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (value_count != grid_size)
        {
            throw std::runtime_error("Periodic ABFS channel size does not match the full response grid.");
        }
        if (GlobalV::MY_RANK != 0)
        {
            channel.potential_r.assign(static_cast<std::size_t>(value_count),
                                       SternheimerRPA::Complex(0.0, 0.0));
        }
        MPI_Bcast(channel.potential_r.data(),
                  value_count,
                  MPI_DOUBLE_COMPLEX,
                  0,
                  MPI_COMM_WORLD);
    }
#else
    (void)channels;
    (void)grid_size;
    (void)enabled;
#endif
}

void broadcast_periodic_abfs(SternheimerPeriodicABFGridData& data,
                             const int grid_size,
                             const bool enabled)
{
    broadcast_periodic_abf_channels(data.densities, grid_size, enabled);
    broadcast_periodic_abf_channels(data.potentials, grid_size, enabled);
}

SternheimerPeriodicABFGridData build_abfs_full_coulomb_bloch_grid_channels(
    const UnitCell& ucell,
    const SternheimerFDHamiltonian::Grid& grid,
    const SternheimerReducedKPoint& qpoint,
    const double gamma_inverse_k2,
    const int max_channels_per_atom,
    const double pca_threshold,
    const SternheimerOrbitalSet* reusable_rpa_abfs,
    const std::chrono::steady_clock::time_point& chi0_start_time)
{
    const SternheimerABFSInput input = build_abfs_density_input(ucell, pca_threshold, reusable_rpa_abfs);
    SternheimerPeriodicABFGridData result;
    result.potentials
        = limit_sternheimer_abf_channels_per_atom(sample_sternheimer_abf_bloch_grid_channels(input.radials_by_type,
                                                                                             input.atom_types,
                                                                                             input.atom_positions,
                                                                                             grid,
                                                                                             qpoint,
                                                                                             -1),
                                                  max_channels_per_atom);
    append_chi0_progress_event("abfs_grid_ready",
                               0,
                               -1,
                               -1,
                               -1,
                               0,
                               nullptr,
                               -1.0,
                               elapsed_seconds_since(chi0_start_time),
                               "channels=" + std::to_string(result.potentials.size()));
    if (sternheimer_grid_coulomb_diagnostic_enabled(static_cast<int>(result.potentials.size())))
    {
        result.densities = result.potentials;
    }
    solve_sternheimer_abf_periodic_full_coulomb_in_place(result.potentials, grid, qpoint, gamma_inverse_k2);
    append_chi0_progress_event("abfs_coulomb_ready",
                               0,
                               -1,
                               -1,
                               -1,
                               0,
                               nullptr,
                               -1.0,
                               elapsed_seconds_since(chi0_start_time),
                               "channels=" + std::to_string(result.potentials.size()));
    return result;
}

SternheimerPeriodicABFGridData build_supercell_translation_full_response_abfs_from_input(
    const SternheimerABFSInput& input,
    const SternheimerFDHamiltonian::Grid& grid,
    const SternheimerSupercellTranslationSum& translation_sum,
    const int max_channels_per_atom)
{
    std::vector<SternheimerABFBlochGridChannel> supercell_densities
        = sample_sternheimer_abf_bloch_grid_channels(input.radials_by_type,
                                                      input.atom_types,
                                                      input.atom_positions,
                                                      grid,
                                                      {0.0, 0.0, 0.0},
                                                      -1);
    if (max_channels_per_atom > 0)
    {
        supercell_densities
            = limit_sternheimer_abf_channels_per_atom(std::move(supercell_densities), max_channels_per_atom);
    }

    SternheimerPeriodicABFGridData result;
    result.potentials = combine_all_sternheimer_supercell_translation_channels(supercell_densities, translation_sum);
    supercell_densities.clear();
    supercell_densities.shrink_to_fit();
    if (sternheimer_grid_coulomb_diagnostic_enabled(static_cast<int>(result.potentials.size())))
    {
        result.densities = result.potentials;
    }
    solve_sternheimer_abf_periodic_full_coulomb_in_place(result.potentials, grid, {0.0, 0.0, 0.0}, 0.0);
    return result;
}

SternheimerPeriodicABFGridData build_supercell_translation_full_response_abfs(
    const UnitCell& ucell,
    const SternheimerFDHamiltonian::Grid& grid,
    const SternheimerSupercellTranslationSum& translation_sum,
    const int max_channels_per_atom,
    const double pca_threshold,
    const SternheimerOrbitalSet* reusable_rpa_abfs)
{
    return build_supercell_translation_full_response_abfs_from_input(
        build_abfs_density_input(ucell, pca_threshold, reusable_rpa_abfs),
        grid,
        translation_sum,
        max_channels_per_atom);
}

siab::Provenance make_siab_production_provenance(const UnitCell& ucell,
                                                 const std::string& auxiliary_basis_sha256,
                                                 const SternheimerRPA::FrequencyGrid& frequency_grid,
                                                 const double pca_threshold,
                                                 const SternheimerCoulombWhitening& whitening)
{
    siab::Provenance provenance;
    provenance.abacus_commit = siab::require_source_commit(compiled_commit_metadata());
    provenance.auxiliary_basis_sha256 = auxiliary_basis_sha256;
    provenance.cell_bohr = cell_vectors_bohr(ucell);
    provenance.ecut_ry = PARAM.inp.ecutwfc;
    provenance.kernel = "full_coulomb";
    const std::vector<std::string> orbital_files
        = siab::resolve_required_input_files(orbital_dir_from_env_or_input(),
                                             orbital_files_from_env_or_cell(ucell),
                                             "initial orbital");
    const std::vector<std::string> pseudopotential_files
        = siab::resolve_required_input_files(PARAM.inp.pseudo_dir, ucell.pseudo_fn, "pseudopotential");
    provenance.orbital_sha256 = siab::sha256_file_manifest(orbital_files);
    provenance.pseudopotential_sha256 = siab::sha256_file_manifest(pseudopotential_files);
    provenance.spin_convention = "spin_resolved_occupation_in_reference_rows";
    provenance.executable_sha256 = siab::sha256_file(siab::resolve_executable_path());
    provenance.exx_pca_thr = pca_threshold;
    provenance.sternheimer_nfreq = static_cast<int>(frequency_grid.omega_ha.size());
    provenance.frequency_ha = frequency_grid.omega_ha;
    provenance.frequency_weights_ha = frequency_grid.weights_ha;
    provenance.mpi_ranks = GlobalV::NPROC;
    provenance.omp_threads = PARAM.globalv.nthread_per_proc;
    provenance.auxiliary_whitening = "global_full_coulomb_v1";
    provenance.raw_auxiliary_dimension = whitening.raw_dimension;
    provenance.whitened_auxiliary_rank = whitening.retained_rank;
    provenance.discarded_auxiliary_rank = whitening.discarded_rank;
    provenance.coulomb_relative_threshold = whitening.relative_threshold;
    provenance.coulomb_eigenvalues = whitening.eigenvalues;
    provenance.coulomb_max_orthonormality_error = whitening.max_orthonormality_error;
    provenance.coulomb_transform_sha256 = hash_coulomb_whitening_transform(whitening);
    return provenance;
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
            sampling_atom.l_nchi[static_cast<std::size_t>(angular_momentum)] = atom_orbitals.getNchi(angular_momentum);
            sampling_atom.nw
                += (2 * angular_momentum + 1) * sampling_atom.l_nchi[static_cast<std::size_t>(angular_momentum)];
        }
        sampling_atom.set_index();

        for (int ia = 0; ia != unitcell_atom.na; ++ia, ++atom_index)
        {
            const int orbital_count = sampling_atom.nw;
            const std::size_t candidate_begin = candidates.size();
            candidates.resize(candidate_begin + static_cast<std::size_t>(orbital_count));
            for (int iw = 0; iw != orbital_count; ++iw)
            {
                SternheimerDeltaGridFunction& candidate = candidates[candidate_begin + static_cast<std::size_t>(iw)];
                candidate.values.assign(static_cast<std::size_t>(grid_size),
                                        SternheimerFDHamiltonian::Complex(0.0, 0.0));
                for (SternheimerFDHamiltonian::Vector& gradient: candidate.gradients)
                {
                    gradient.assign(static_cast<std::size_t>(grid_size), SternheimerFDHamiltonian::Complex(0.0, 0.0));
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
                                                              {atom_position.x, atom_position.y, atom_position.z},
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
                            = ModuleGint::Vec3d(ix * grid.hx - atom_position.x - image[0] * lx,
                                                iy * grid.hy - atom_position.y - image[1] * ly,
                                                iz * grid.hz - atom_position.z - image[2] * lz);
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

std::vector<SternheimerDeltaGridFunction> build_lcao_grid_functions_from_coefficients(
    const UnitCell& ucell,
    const SternheimerFDHamiltonian::Grid& grid,
    const LCAO_Orbitals& orb,
    const SternheimerReducedKPoint& kpoint,
    const std::vector<std::vector<SternheimerFDHamiltonian::Complex>>& coefficients)
{
    if (coefficients.empty())
    {
        return {};
    }
    const int grid_size = grid.size();
    constexpr int sample_chunk_size = 8192;
    std::vector<SternheimerDeltaGridFunction> functions(coefficients.size());
    for (SternheimerDeltaGridFunction& function: functions)
    {
        function.values.assign(static_cast<std::size_t>(grid_size),
                               SternheimerFDHamiltonian::Complex(0.0, 0.0));
        for (SternheimerFDHamiltonian::Vector& gradient: function.gradients)
        {
            gradient.assign(static_cast<std::size_t>(grid_size),
                            SternheimerFDHamiltonian::Complex(0.0, 0.0));
        }
    }

    std::size_t ao_begin = 0;
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
            for (const auto& state_coefficients: coefficients)
            {
                if (ao_begin > state_coefficients.size()
                    || static_cast<std::size_t>(orbital_count) > state_coefficients.size() - ao_begin)
                {
                    throw std::runtime_error(
                        "Sternheimer direct LCAO sampling coefficient basis is smaller than the sampled AO basis.");
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
                    accumulate_delta_sternheimer_lcao_state_samples(values,
                                                                    gradients,
                                                                    chunk_size,
                                                                    orbital_count,
                                                                    static_cast<std::size_t>(first),
                                                                    ao_begin,
                                                                    coefficients,
                                                                    kpoint,
                                                                    image,
                                                                    functions);
                }
            }
            ao_begin += static_cast<std::size_t>(orbital_count);
        }
    }
    for (const auto& state_coefficients: coefficients)
    {
        if (state_coefficients.size() != ao_begin)
        {
            throw std::runtime_error(
                "Sternheimer direct LCAO sampling coefficient basis exceeds the sampled AO basis.");
        }
    }
    return functions;
}

struct SternheimerSampledLCAOKPoint
{
    SternheimerFDZeroOrderStates states;
    SternheimerFDZeroOrderStates unoccupied_states;
    std::vector<SternheimerDeltaGridFunction> occupied_functions;
    std::vector<SternheimerDeltaGridFunction> unoccupied_functions;
    std::vector<SternheimerDeltaGridFunction> occupied_projector_functions;
    double occupied_raw_norm_min = std::numeric_limits<double>::infinity();
    double occupied_raw_norm_max = 0.0;
    double unoccupied_raw_norm_min = std::numeric_limits<double>::infinity();
    double unoccupied_raw_norm_max = 0.0;
    std::vector<double> occupied_raw_norms;
    std::vector<double> unoccupied_raw_norms;
};

SternheimerSampledLCAOKPoint sample_sternheimer_lcao_kpoint(
    const UnitCell& ucell,
    const SternheimerFDHamiltonian::Grid& grid,
    const LCAO_Orbitals& orbitals,
    const SternheimerLCAOOccupiedKPoint& record,
    const int requested_occupied_states,
    const bool include_unoccupied,
    const int requested_unoccupied_states,
    const double volume_element,
    const double norm_tolerance)
{
    SternheimerSampledLCAOKPoint sampled;
    const std::size_t occupied_count
        = sternheimer_sampled_occupied_count(record.coefficients.size(), requested_occupied_states);
    const std::size_t sampled_unoccupied_count
        = sternheimer_sampled_unoccupied_count(include_unoccupied,
                                                record.unoccupied_coefficients.size(),
                                                requested_unoccupied_states);
    std::vector<std::vector<SternheimerFDHamiltonian::Complex>> all_coefficients;
    all_coefficients.reserve(occupied_count + sampled_unoccupied_count);
    all_coefficients.insert(all_coefficients.end(),
                            record.coefficients.begin(),
                            record.coefficients.begin() + static_cast<std::ptrdiff_t>(occupied_count));
    if (sampled_unoccupied_count > 0)
    {
        all_coefficients.insert(all_coefficients.end(),
                                record.unoccupied_coefficients.begin(),
                                record.unoccupied_coefficients.begin()
                                    + static_cast<std::ptrdiff_t>(sampled_unoccupied_count));
    }
    std::vector<SternheimerDeltaGridFunction> all_functions
        = build_lcao_grid_functions_from_coefficients(
            ucell, grid, orbitals, sternheimer_lcao_grid_kpoint(record), all_coefficients);

    if (all_functions.size() != occupied_count + sampled_unoccupied_count)
    {
        throw std::runtime_error("Sternheimer periodic direct LCAO sampling lost selected KS states.");
    }
    sampled.states.eigenvalues.reserve(occupied_count);
    sampled.states.wavefunctions.reserve(occupied_count);
    sampled.states.residual_norms.reserve(occupied_count);
    sampled.occupied_functions.reserve(occupied_count);
    for (std::size_t ib = 0; ib != occupied_count; ++ib)
    {
        SternheimerDeltaGridFunction occupied_function = std::move(all_functions[ib]);
        const double norm = sternheimer_fd_grid_norm(occupied_function.values, volume_element);
        if (norm <= norm_tolerance)
        {
            throw std::runtime_error("Sternheimer periodic sampled occupied function has zero norm.");
        }
        sampled.occupied_raw_norm_min = std::min(sampled.occupied_raw_norm_min, norm);
        sampled.occupied_raw_norm_max = std::max(sampled.occupied_raw_norm_max, norm);
        sampled.occupied_raw_norms.push_back(norm);
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
    sampled.unoccupied_states.eigenvalues.reserve(sampled_unoccupied_count);
    sampled.unoccupied_functions.reserve(sampled_unoccupied_count);
    for (std::size_t ib = 0; ib != sampled_unoccupied_count; ++ib)
    {
        SternheimerDeltaGridFunction unoccupied_function
            = std::move(all_functions[occupied_count + ib]);
        const double norm = sternheimer_fd_grid_norm(unoccupied_function.values, volume_element);
        if (norm <= norm_tolerance)
        {
            throw std::runtime_error("Sternheimer periodic sampled unoccupied function has zero norm.");
        }
        sampled.unoccupied_raw_norm_min = std::min(sampled.unoccupied_raw_norm_min, norm);
        sampled.unoccupied_raw_norm_max = std::max(sampled.unoccupied_raw_norm_max, norm);
        sampled.unoccupied_raw_norms.push_back(norm);
        const SternheimerFDHamiltonian::Complex inverse_norm(1.0 / norm, 0.0);
        for (auto& value: unoccupied_function.values)
        {
            value *= inverse_norm;
        }
        for (auto& gradient: unoccupied_function.gradients)
        {
            for (auto& value: gradient)
            {
                value *= inverse_norm;
            }
        }
        sampled.unoccupied_states.eigenvalues.push_back(record.unoccupied_eigenvalues[ib]);
        sampled.unoccupied_functions.push_back(std::move(unoccupied_function));
    }
    sampled.occupied_projector_functions = orthonormalize_delta_sternheimer_grid_functions(
        sampled.occupied_functions, volume_element, norm_tolerance);
    if (sampled.occupied_projector_functions.size() != occupied_count)
    {
        throw std::runtime_error("Sternheimer periodic occupied projector lost a linearly dependent KS state.");
    }
    return sampled;
}

void run_sternheimer_spectrum_audit(
    const UnitCell& ucell,
    const LCAO_Orbitals& orbitals,
    const SternheimerABACUSFDGridData& grid_data,
    const std::vector<double>& local_potential,
    const SternheimerLCAOOccupiedKPoint& source_record,
    const SternheimerLCAOOccupiedKPoint& target_record,
    const int source_index,
    const int target_index)
{
    auto target = sample_sternheimer_lcao_kpoint(
        ucell, grid_data.grid, orbitals, target_record, 0, true, 0,
        grid_data.volume_element, PARAM.inp.sternheimer_delta_norm_tol);
    if (target.unoccupied_functions.empty() || source_record.eigenvalues.empty())
    {
        throw std::runtime_error("Spectrum audit requires the complete occupied and virtual KS records.");
    }
    auto target_grid = grid_data;
    target_grid.grid.kpoint = sternheimer_lcao_grid_kpoint(target_record);
    auto projectors = make_sternheimer_fd_nonlocal_projector_from_unitcell(
        ucell, target_grid.grid, target_grid.volume_element);
    const auto hamiltonian = make_sternheimer_fd_hamiltonian_from_local_potential(
        target_grid, local_potential, 1.0, std::move(projectors), PARAM.inp.sternheimer_fd_order);
    std::ofstream audit("STERNHEIMER_SPECTRUM_AUDIT.dat");
    if (!audit)
    {
        throw std::runtime_error("Cannot open spectrum audit output.");
    }
    audit << std::setprecision(17)
          << "status diagnostic_only\n"
          << "physical_result no\n"
          << "response_equations 0\n"
          << "grid " << grid_data.grid.nx << ' ' << grid_data.grid.ny << ' '
          << grid_data.grid.nz << '\n'
          << "source_target " << source_index + 1 << ' ' << target_index + 1 << '\n'
          << "occupied_raw_norm_range " << target.occupied_raw_norm_min << ' '
          << target.occupied_raw_norm_max << '\n'
          << "virtual_raw_norm_range " << target.unoccupied_raw_norm_min << ' '
          << target.unoccupied_raw_norm_max << '\n';
    std::vector<double> eigenvalues = target_record.eigenvalues;
    eigenvalues.insert(eigenvalues.end(), target_record.unoccupied_eigenvalues.begin(),
                       target_record.unoccupied_eigenvalues.end());
    const std::size_t occupied_count = target.occupied_functions.size();
    std::vector<double> raw_norms = target.occupied_raw_norms;
    raw_norms.insert(raw_norms.end(), target.unoccupied_raw_norms.begin(), target.unoccupied_raw_norms.end());
    auto functions = std::move(target.occupied_functions);
    functions.insert(functions.end(), std::make_move_iterator(target.unoccupied_functions.begin()),
                     std::make_move_iterator(target.unoccupied_functions.end()));
    target.unoccupied_functions.clear();
    {
        const auto matrices = assemble_delta_sternheimer_grid_matrices_fast(
            hamiltonian, functions, grid_data.volume_element);
        audit << "# state pbe_Ry overlap kinetic_Ry local_Ry nonlocal_Ry expectation_Ry raw_norm\n";
        for (std::size_t state = 0; state < functions.size(); ++state)
        {
            const std::size_t index = state * (functions.size() + 1);
            const double norm2 = matrices.overlap[index].real();
            audit << "state " << state + 1 << ' ' << eigenvalues.at(state) << ' '
                  << norm2 << ' ' << matrices.kinetic[index].real() / norm2 << ' '
                  << matrices.local_potential[index].real() / norm2 << ' '
                  << matrices.nonlocal[index].real() / norm2 << ' '
                  << matrices.hamiltonian[index].real() / norm2 << ' ' << raw_norms.at(state) << '\n';
        }
        std::ofstream matrix_out("STERNHEIMER_SPECTRUM_MATRICES.dat");
        if (!matrix_out)
        {
            throw std::runtime_error("Cannot open spectrum matrices output.");
        }
        matrix_out << std::setprecision(17) << "# row col S T Vloc Vnl H (real imag), column major\n";
        for (std::size_t column = 0; column < functions.size(); ++column)
        {
            for (std::size_t row = 0; row < functions.size(); ++row)
            {
                const auto index = row + column * functions.size();
                matrix_out << row + 1 << ' ' << column + 1;
                for (const auto* matrix : {&matrices.overlap, &matrices.kinetic,
                                          &matrices.local_potential, &matrices.nonlocal,
                                          &matrices.hamiltonian})
                {
                    matrix_out << ' ' << (*matrix)[index].real() << ' ' << (*matrix)[index].imag();
                }
                matrix_out << '\n';
            }
        }
        matrix_out.close();
        if (!matrix_out)
        {
            throw std::runtime_error("Spectrum matrix output did not complete.");
        }
    }
    audit.flush();
    if (!audit)
    {
        throw std::runtime_error("Spectrum state output did not complete.");
    }
    functions.erase(functions.begin(), functions.begin() + occupied_count);
    SternheimerDeltaSubspaceOptions options;
    options.max_virtual_states = 0;
    options.norm_tolerance = PARAM.inp.sternheimer_delta_norm_tol;
    options.retain_grid_functions = false;
    options.evaluate_full_grid_difference = false;
    const auto subspace = build_reference_delta_sternheimer_subspace(
        hamiltonian, target.occupied_projector_functions, std::move(functions),
        grid_data.volume_element, options);
    const double source_vbm = *std::max_element(source_record.eigenvalues.begin(),
                                                source_record.eigenvalues.end());
    const double pbe_cbm = *std::min_element(target_record.unoccupied_eigenvalues.begin(),
                                            target_record.unoccupied_eigenvalues.end());
    audit << "source_vbm_Ry " << source_vbm << '\n'
          << "target_cbm_Ry " << pbe_cbm << '\n'
          << "pbe_transition_min_Ry " << pbe_cbm - source_vbm << '\n';
    double delta_min = std::numeric_limits<double>::infinity();
    for (std::size_t state = 0; state < subspace.virtual_states.size(); ++state)
    {
        const double energy = subspace.virtual_states[state].eigenvalue;
        delta_min = std::min(delta_min, energy);
        audit << "delta_state " << state + 1 << ' ' << energy << '\n';
    }
    audit << "delta_transition_min_Ry " << delta_min - source_vbm << '\n'
          << "audit_complete yes\n";
    audit.close();
    if (!audit)
    {
        throw std::runtime_error("Spectrum audit output did not complete.");
    }
}

void run_sternheimer_galerkin_spectrum_audit(
    const UnitCell& ucell,
    const LCAO_Orbitals& orbitals,
    SternheimerABACUSFDGridData fine_grid,
    SternheimerABACUSFDGridData response_grid,
    const std::vector<double>& fine_potential,
    const SternheimerLCAOOccupiedKPoint& source_record,
    const SternheimerLCAOOccupiedKPoint& target_record,
    const int source_index,
    const int target_index)
{
    fine_grid.grid.kpoint = sternheimer_lcao_grid_kpoint(target_record);
    response_grid.grid.kpoint = fine_grid.grid.kpoint;
    auto target = sample_sternheimer_lcao_kpoint(
        ucell, fine_grid.grid, orbitals, target_record, 0, true, 0,
        fine_grid.volume_element, PARAM.inp.sternheimer_delta_norm_tol);
    if (target.unoccupied_functions.empty() || source_record.eigenvalues.empty())
    {
        throw std::runtime_error("Galerkin spectrum audit requires complete occupied and virtual KS records.");
    }
    const auto occupied_count = target.occupied_functions.size();
    auto raw_norms = target.occupied_raw_norms;
    raw_norms.insert(raw_norms.end(), target.unoccupied_raw_norms.begin(), target.unoccupied_raw_norms.end());
    auto functions = std::move(target.occupied_functions);
    functions.insert(functions.end(), std::make_move_iterator(target.unoccupied_functions.begin()),
                     std::make_move_iterator(target.unoccupied_functions.end()));
    target = SternheimerSampledLCAOKPoint();
    std::vector<const SternheimerFDHamiltonian::Vector*> fine_states;
    for (auto& function : functions)
    {
        // This diagnostic compares FullGrid FD8 forms, not the old weak-A block.
        for (auto& gradient : function.gradients)
        {
            SternheimerFDHamiltonian::Vector().swap(gradient);
        }
        fine_states.push_back(&function.values);
    }
    auto nonlocal = make_sternheimer_fd_nonlocal_projector_from_unitcell(
        ucell, fine_grid.grid, fine_grid.volume_element);
    auto fine_hamiltonian = std::make_shared<const SternheimerFDHamiltonian>(
        make_sternheimer_fd_hamiltonian_from_local_potential(
            fine_grid, fine_potential, 1.0, std::move(nonlocal), PARAM.inp.sternheimer_fd_order));
    trim_sternheimer_process_heap();
    std::ofstream audit("STERNHEIMER_GALERKIN_AUDIT.dat");
    if (!audit)
    {
        throw std::runtime_error("Cannot open Galerkin spectrum audit.");
    }
    audit << std::setprecision(17)
          << "status diagnostic_only\nphysical_result no\nresponse_equations 0\n"
          << "reference FullGrid_FD8_not_weak_A\nfull_pc_positivity_certified no\n"
          << "fine_grid " << fine_grid.grid.nx << ' ' << fine_grid.grid.ny << ' ' << fine_grid.grid.nz << '\n'
          << "coarse_grid " << response_grid.grid.nx << ' ' << response_grid.grid.ny << ' '
          << response_grid.grid.nz << '\n'
          << "source_target " << source_index + 1 << ' ' << target_index + 1 << '\n'
          << "nocc " << occupied_count << "\nnvirtual " << functions.size() - occupied_count << '\n'
          << "source_vbm_Ry " << *std::max_element(source_record.eigenvalues.begin(), source_record.eigenvalues.end())
          << '\n';
    std::vector<double> pbe_eigenvalues = target_record.eigenvalues;
    pbe_eigenvalues.insert(pbe_eigenvalues.end(), target_record.unoccupied_eigenvalues.begin(),
                          target_record.unoccupied_eigenvalues.end());
    audit << "# state case index pbe_Ry fine_norm2 projected_norm2 reconstruction_error raw_sampled_norm\n";
    constexpr std::size_t audit_workspace_budget = 32ULL * 1024 * 1024 * 1024;
    for (const bool identity : {true, false})
    {
        const std::string label = identity ? "fine" : "projected";
        const auto matrices = audit_sternheimer_galerkin_matrices(
            fine_hamiltonian, fine_states, identity ? fine_grid.grid : response_grid.grid,
            audit_workspace_budget);
        const auto n = matrices.dimension;
        if (n != functions.size())
        {
            throw std::runtime_error("Galerkin audit changed the full KS state count.");
        }
        std::ofstream output("STERNHEIMER_GALERKIN_" + label + "_MATRICES.dat");
        if (!output)
        {
            throw std::runtime_error("Cannot open Galerkin matrix output.");
        }
        output << std::setprecision(17) << "# row col S T Vloc Vnl H (real imag), column major\n";
        for (std::size_t column = 0; column != n; ++column)
        {
            for (std::size_t row = 0; row != n; ++row)
            {
                const auto index = row + column * n;
                output << row + 1 << ' ' << column + 1;
                for (const auto* matrix : {&matrices.overlap, &matrices.kinetic, &matrices.local_potential,
                                          &matrices.nonlocal, &matrices.hamiltonian})
                {
                    output << ' ' << matrix->at(index).real() << ' ' << matrix->at(index).imag();
                }
                output << '\n';
            }
        }
        output.close();
        if (!output)
        {
            throw std::runtime_error("Galerkin matrix output did not complete.");
        }
        audit << "case " << label << " workspace_bytes " << matrices.workspace_bytes
              << " filtered_identity_max_abs_error " << matrices.filtered_matrix_identity_max_abs_error << '\n';
        for (std::size_t state = 0; state != n; ++state)
        {
            audit << "state " << label << ' ' << state + 1 << ' ' << pbe_eigenvalues.at(state) << ' '
                  << matrices.fine_norm_squared.at(state) << ' ' << matrices.projected_norm_squared.at(state)
                  << ' ' << matrices.reconstruction_relative_error.at(state) << ' ' << raw_norms.at(state) << '\n';
        }
        audit.flush();
        if (!audit)
        {
            throw std::runtime_error("Galerkin state output did not complete.");
        }
    }
    audit << "audit_complete yes\n";
    audit.close();
    if (!audit)
    {
        throw std::runtime_error("Galerkin audit output did not complete.");
    }
}

void run_sternheimer_weak_response_audit(
    const UnitCell& ucell,
    const LCAO_Orbitals& orbitals,
    SternheimerABACUSFDGridData fine_grid,
    SternheimerABACUSFDGridData coarse_grid,
    const std::vector<double>& fine_potential,
    const SternheimerLCAOOccupiedKPoint& source_record,
    const SternheimerLCAOOccupiedKPoint& target_record,
    const SternheimerReducedKPoint& qpoint,
    const double omega_ha,
    const double pca_threshold,
    const SternheimerOrbitalSet* rpa_abfs,
    const int source_index,
    const int target_index)
{
    using Vector = SternheimerFDHamiltonian::Vector;
    using Complex = SternheimerFDHamiltonian::Complex;
    using Blocks = SternheimerWeakAugmented;
    const bool matrix_audit = env_is_true("ABACUS_STERNHEIMER_WEAK_MATRIX_AUDIT");
    const auto started = std::chrono::steady_clock::now();
    std::ofstream audit("STERNHEIMER_WEAK_RESPONSE_AUDIT.dat");
    if (!audit) throw std::runtime_error("Cannot open fine-weak response audit.");
    audit << std::setprecision(17)
          << "status diagnostic_only\nphysical_result no\nfull_q_response no\n"
          << "reference fine_analytic_weak_form\nsource_band_selection VBM_only\n"
          << (matrix_audit ? "auxiliary_channel_selection all\nsymmetry_restoration no\n"
                           : "auxiliary_channel_selection first_two_only\nsymmetry_restoration no\n")
          << "fine_grid " << fine_grid.grid.nx << ' ' << fine_grid.grid.ny << ' ' << fine_grid.grid.nz << '\n'
          << "coarse_grid " << coarse_grid.grid.nx << ' ' << coarse_grid.grid.ny << ' ' << coarse_grid.grid.nz << '\n'
          << "source_target " << source_index + 1 << ' ' << target_index + 1 << '\n'
          << "omega_Ha " << omega_ha << '\n';
    if (matrix_audit)
        audit << "matrix_audit yes\nresponse_matrix partial\n"
              << "global_matrix_checks deferred_external_merge\n"
              << "qpoint " << qpoint[0] << ' ' << qpoint[1] << ' ' << qpoint[2] << '\n';
    const auto stage = [&](const std::string& name, const double solve_s = -1.0,
                           const double contraction_s = -1.0, const double finecheck_s = -1.0) {
        const auto memory = detect_sternheimer_memory_snapshot();
        audit << "stage " << name << " elapsed_s " << elapsed_seconds_since(started)
              << " memory_current_bytes " << memory.current_bytes;
        if (solve_s >= 0.0)
            audit << " solve_s " << solve_s << " contraction_s " << contraction_s << " finecheck_s " << finecheck_s;
        audit << '\n';
        audit.flush();
        if (!audit) throw std::runtime_error("Fine-weak response audit write failed.");
    };
    const auto source_k = sternheimer_lcao_grid_kpoint(source_record);
    fine_grid.grid.kpoint = sternheimer_lcao_grid_kpoint(target_record);
    coarse_grid.grid.kpoint = fine_grid.grid.kpoint;
    for (int a = 0; a < 3; ++a)
        if (std::abs(fine_grid.grid.kpoint[a] - source_k[a] - qpoint[a]) > 1e-12)
            throw std::runtime_error("Fine-weak audit does not yet admit BZ-folded k pairs.");
    if (source_record.eigenvalues.empty() || target_record.unoccupied_eigenvalues.empty()
        || !(omega_ha > 0.0) || !std::isfinite(omega_ha))
        throw std::runtime_error("Fine-weak audit requires complete KS records and positive frequency.");
    if (target_record.eigenvalues.size() + target_record.unoccupied_eigenvalues.size()
            != static_cast<std::size_t>(PARAM.globalv.nlocal)
        || target_record.coefficients.size() != target_record.eigenvalues.size()
        || target_record.unoccupied_coefficients.size() != target_record.unoccupied_eigenvalues.size())
        throw std::runtime_error("Fine-weak audit requires nbands=nlocal and every target coefficient vector.");
    if (matrix_audit)
    {
        if (source_record.coefficients.size() != source_record.eigenvalues.size()
            || source_record.occupations.size() != source_record.eigenvalues.size()
            || target_record.eigenvalues.empty()
            || target_record.occupations.size() != target_record.eigenvalues.size()
            || fine_potential.size() != static_cast<std::size_t>(fine_grid.grid.size()))
            throw std::runtime_error("Weak matrix audit requires complete source/target KS records and full fine potential.");
        const auto finite_real = [](const double value) { return std::isfinite(value); };
        for (const auto* values : {&source_record.eigenvalues, &target_record.eigenvalues,
                                   &target_record.unoccupied_eigenvalues, &fine_potential})
            if (!std::all_of(values->begin(), values->end(), finite_real))
                throw std::runtime_error("Weak matrix audit input contains nonfinite real values.");
        for (const auto* coefficients : {&source_record.coefficients, &target_record.coefficients,
                                         &target_record.unoccupied_coefficients})
            for (const auto& state : *coefficients)
                if (state.size() != static_cast<std::size_t>(PARAM.globalv.nlocal)
                    || !std::all_of(state.begin(), state.end(), [](const Complex value) {
                           return std::isfinite(value.real()) && std::isfinite(value.imag());
                       }))
                    throw std::runtime_error("Weak matrix audit requires complete finite KS coefficient vectors.");
        for (int a = 0; a < 3; ++a)
            if (!std::isfinite(source_k[a]) || !std::isfinite(fine_grid.grid.kpoint[a]) || !std::isfinite(qpoint[a]))
                throw std::runtime_error("Weak matrix audit requires finite k/q coordinates.");
    }
    auto target = sample_sternheimer_lcao_kpoint(
        ucell, fine_grid.grid, orbitals, target_record, 0, true, 0,
        fine_grid.volume_element, PARAM.inp.sternheimer_delta_norm_tol);
    const int nocc = static_cast<int>(target.occupied_functions.size());
    const int nvirtual = static_cast<int>(target.unoccupied_functions.size());
    if (nocc != static_cast<int>(target_record.eigenvalues.size())
        || nvirtual != static_cast<int>(target_record.unoccupied_eigenvalues.size()))
        throw std::runtime_error("Fine-weak audit lost target KS states.");
    auto functions = std::move(target.occupied_functions);
    functions.insert(functions.end(), std::make_move_iterator(target.unoccupied_functions.begin()),
                     std::make_move_iterator(target.unoccupied_functions.end()));
    target = SternheimerSampledLCAOKPoint();
    stage("fine_states_sampled");
    const double orthogonality_error = orthonormalize_sternheimer_weak_states_in_place(functions, fine_grid.volume_element);
    audit << "orthonormalization_max_error " << orthogonality_error << '\n';
    if (orthogonality_error > 1e-8)
        throw std::runtime_error("Fine-weak state orthonormalization failed its numerical gate.");
    stage("fine_states_orthonormal");
    auto nonlocal = make_sternheimer_fd_nonlocal_projector_from_unitcell(
        ucell, fine_grid.grid, fine_grid.volume_element);
    auto h = std::make_shared<const SternheimerFDHamiltonian>(
        make_sternheimer_fd_hamiltonian_from_local_potential(
            fine_grid, fine_potential, 1.0, std::move(nonlocal), PARAM.inp.sternheimer_fd_order));
    auto op = std::make_shared<SternheimerWeakGridOperator>(h, coarse_grid.grid);
    auto assembled = op->assemble_blocks(functions);
    const double metric_error = assembled.maximum_metric_error;
    std::ofstream hu_output("STERNHEIMER_WEAK_HU.dat");
    hu_output << std::setprecision(17) << "# nocc nvirtual, then row col H_Ry_real H_Ry_imag\n"
              << nocc << ' ' << nvirtual << '\n';
    for (std::size_t j = 0; j < functions.size(); ++j)
        for (std::size_t i = 0; i < functions.size(); ++i)
        {
            const auto value = assembled.state_hamiltonian[i + functions.size() * j];
            hu_output << i + 1 << ' ' << j + 1 << ' ' << value.real() << ' ' << value.imag() << '\n';
        }
    hu_output.close();
    if (!hu_output) throw std::runtime_error("Fine-weak Hamiltonian audit write failed.");
    Blocks::Data data;
    data.nocc = nocc; data.nvirtual = nvirtual; data.ncoarse = coarse_grid.grid.size();
    data.hu = std::move(assembled.state_hamiltonian);
    data.l = std::move(assembled.state_coarse_overlap);
    data.k = std::move(assembled.state_coarse_hamiltonian);
    auto blocks = std::make_shared<const Blocks>(std::move(data));
    audit << "nocc " << nocc << "\nnvirtual " << nvirtual
          << "\nfine_state_metric_error " << metric_error
          << "\ncomplement_metric_min " << blocks->minimum_metric_eigenvalue() << '\n';
    // Couplings retain every analytic gradient contribution; only now release gradients.
    for (auto& function : functions)
        for (auto& gradient : function.gradients) Vector().swap(gradient);
    trim_sternheimer_process_heap();
    stage("weak_blocks_ready");
    auto source_grid = fine_grid.grid;
    source_grid.kpoint = source_k;
    auto source = sample_sternheimer_lcao_kpoint(
        ucell, source_grid, orbitals, source_record, 0, false, 0,
        fine_grid.volume_element, PARAM.inp.sternheimer_delta_norm_tol);
    const auto vbm = std::max_element(source_record.eigenvalues.begin(), source_record.eigenvalues.end());
    const auto band = std::size_t(vbm - source_record.eigenvalues.begin());
    const double epsilon = *vbm;
    Vector source_values = std::move(source.occupied_functions.at(band).values);
    source = SternheimerSampledLCAOKPoint();
    audit << "source_band " << band + 1 << "\nsource_epsilon_Ry " << epsilon << '\n';
    const auto input = build_abfs_density_input(ucell, pca_threshold, rpa_abfs);
    if (matrix_audit)
    {
        const auto metadata = describe_sternheimer_abf_grid_channels(
            input.radials_by_type, input.atom_types, input.atom_positions, -1);
        if (metadata.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
            throw std::length_error("Weak matrix audit channel count exceeds BLAS integer range.");
        const int channel_count = static_cast<int>(metadata.size());
        const auto shard = SternheimerWeakAuditShard::parse(
            std::getenv("WEAK_AUDIT_SHARD_INDEX"), std::getenv("WEAK_AUDIT_SHARD_COUNT"), channel_count);
        const auto owned_columns = shard.columns();
        audit << "shard_index " << shard.index << "\nshard_count " << shard.count
              << "\nauxiliary_channels " << channel_count << '\n'
              << "# auxiliary_channel IDs are 1-based; named metadata indices are 0-based; M is real-Ylm index\n";
        for (int j = 0; j < channel_count; ++j)
        {
            const auto& channel = metadata[j];
            if (channel.channel_index != j)
                throw std::runtime_error("Weak matrix audit requires contiguous global auxiliary channel IDs.");
            audit << "auxiliary_channel " << j + 1 << " atom_index " << channel.atom_index
                  << " atom_local_index " << channel.atom_local_index << " type_index " << channel.type_index
                  << " L " << channel.angular_momentum << " M " << channel.magnetic_index
                  << " radial_index " << channel.radial_index << " label " << std::quoted(channel.label) << '\n';
        }
        stage("auxiliary_metadata_ready");
        const auto sample_fine_vertex = [&](const int j) {
            auto single = sample_sternheimer_abf_bloch_grid_channels(
                input.radials_by_type, input.atom_types, input.atom_positions, fine_grid.grid, qpoint, 1, j);
            if (single.size() != 1 || single.front().channel_index != j)
                throw std::runtime_error("Weak matrix audit failed to sample the requested global auxiliary channel.");
            solve_sternheimer_abf_periodic_full_coulomb_in_place(single, fine_grid.grid, qpoint, 0.0);
            Vector g = std::move(single.front().potential_r);
            if (g.size() != source_values.size())
                throw std::runtime_error("Weak matrix audit fine vertex dimension mismatch.");
            for (std::size_t ir = 0; ir < g.size(); ++ir)
            {
                // In-place Ha -> Ry multiplication avoids a second fine vertex allocation.
                g[ir] *= 2.0 * source_values[ir];
                if (!std::isfinite(g[ir].real()) || !std::isfinite(g[ir].imag()))
                    throw std::runtime_error("Weak matrix audit fine vertex is nonfinite.");
            }
            return g;
        };
        const int ncoarse = blocks->ncoarse();
        const auto checked_elements = [](const std::size_t rows, const std::size_t columns) {
            if (columns != 0 && rows > Vector().max_size() / columns)
                throw std::length_error("Weak matrix audit projected storage size overflow.");
            return rows * columns;
        };
        // The all-channel cache contains only gF/gW, packed column-major for BLAS.
        Blocks::Vertices projected_vertices;
        projected_vertices.f.resize(checked_elements(nvirtual, channel_count));
        projected_vertices.w.resize(checked_elements(ncoarse, channel_count));
        audit << "projected_f_bytes " << projected_vertices.f.size() * sizeof(Complex)
              << "\nprojected_w_bytes " << projected_vertices.w.size() * sizeof(Complex)
              << "\nfine_vertex_cache_channels 1\n";
        const Vector first_row_vertex = sample_fine_vertex(0);
        for (int j = 0; j < channel_count; ++j)
        {
            {
                const Vector sampled = j != 0 ? sample_fine_vertex(j) : Vector();
                const Vector& g = j == 0 ? first_row_vertex : sampled;
                Vector gu(functions.size()), ge;
                for (std::size_t u = 0; u < functions.size(); ++u)
                    gu[u] = sternheimer_fd_grid_dot(functions[u].values, g, fine_grid.volume_element);
                op->project(g, ge);
                const auto vertex = blocks->project_vertices(gu, ge);
                std::copy(vertex.f.begin(), vertex.f.end(),
                          projected_vertices.f.begin() + std::size_t(j) * nvirtual);
                std::copy(vertex.w.begin(), vertex.w.end(),
                          projected_vertices.w.begin() + std::size_t(j) * ncoarse);
            }
            // Only the first-row fine vertex survives for the repeated selected-row check.
            if ((j + 1) % 32 == 0 || j + 1 == channel_count)
            {
                stage("vertices_projected_" + std::to_string(j + 1));
            }
        }
        trim_sternheimer_process_heap();
        stage("fine_vertices_ready");
        SternheimerRPA::SolverOptions options;
        options.max_iter = positive_int_from_env("ABACUS_STERNHEIMER_WEAK_AUDIT_MAX_ITER", 300);
        options.residual_tol = 1e-8;
        options.use_fd_spectral_preconditioner = false;
        audit << "coordinate_transform complement_inverse_sqrt\npreconditioner none\nmax_iter "
              << options.max_iter << "\nresidual_tolerance " << options.residual_tol << '\n';
        Vector partial_response(checked_elements(channel_count, owned_columns.size()), Complex(0.0));
        std::size_t equations = 0;
        for (int sign : {1, -1})
        {
            Blocks::Worker worker(blocks, [op](const Vector& x, Vector& result) { op->apply(x, result); },
                                  epsilon, sign * 2.0 * omega_ha);
            for (std::size_t owned = 0; owned < owned_columns.size(); ++owned)
            {
                const int j = owned_columns[owned];
                const auto f_begin = projected_vertices.f.begin() + std::size_t(j) * nvirtual;
                const auto w_begin = projected_vertices.w.begin() + std::size_t(j) * ncoarse;
                const Blocks::Vertices vertex{Vector(f_begin, f_begin + nvirtual), Vector(w_begin, w_begin + ncoarse)};
                const auto solve_started = std::chrono::steady_clock::now();
                const auto result = worker.solve(vertex, options);
                const double solve_s = elapsed_seconds_since(solve_started);
                audit << "equation " << ++equations << " channel " << j + 1 << " sign " << sign
                      << " converged " << (result.converged ? "yes" : "no") << " iterations " << result.schur.iterations
                      << " full_residual " << result.relative_residual << " elapsed_s " << elapsed_seconds_since(started) << '\n';
                audit.flush();
                if (!audit) throw std::runtime_error("Weak matrix audit equation output failed.");
                if (!result.converged) throw std::runtime_error("Fine-weak response audit equation did not converge.");
                const auto contraction_started = std::chrono::steady_clock::now();
                Vector block_column(channel_count, Complex(0.0));
                BlasConnector::gemv('C', nvirtual, channel_count, Complex(0.5), projected_vertices.f.data(), nvirtual,
                                    result.coefficients.data(), 1, Complex(0.0), block_column.data(), 1);
                BlasConnector::gemv('C', ncoarse, channel_count, Complex(0.5), projected_vertices.w.data(), ncoarse,
                                    result.coefficients.data() + nvirtual, 1, Complex(1.0), block_column.data(), 1);
                const double contraction_s = elapsed_seconds_since(contraction_started);
                const auto finecheck_started = std::chrono::steady_clock::now();
                const auto expansion = blocks->expand_coordinates(result.coefficients);
                Vector fine_response;
                op->lift(expansion.e, fine_response);
                for (std::size_t u = 0; u < functions.size(); ++u)
                {
#pragma omp parallel for schedule(static)
                    for (std::size_t ir = 0; ir < fine_response.size(); ++ir)
                        fine_response[ir] += expansion.u[u] * functions[u].values[ir];
                }
                const Vector diagonal_vertex = j != 0 ? sample_fine_vertex(j) : Vector();
                for (int i = 0; i < channel_count; ++i)
                {
                    const Complex block_value = block_column[i];
                    if (!std::isfinite(block_value.real()) || !std::isfinite(block_value.imag()))
                        throw std::runtime_error("Weak matrix audit projected response is nonfinite.");
                    if (i == j || i == 0)
                    {
                        // The first row is cached; only the solved column is regenerated.
                        const Complex fine_value = 0.5 * sternheimer_fd_grid_dot(
                            i == 0 ? first_row_vertex : diagonal_vertex, fine_response, fine_grid.volume_element);
                        const double error = std::abs(fine_value - block_value) / std::max(1.0, std::abs(fine_value));
                        audit << "vertex_check " << equations << ' ' << i + 1 << ' ' << error << '\n';
                        if (!std::isfinite(error) || error > 1e-9)
                            throw std::runtime_error("Fine-weak left/right vertex consistency failed.");
                    }
                    partial_response[std::size_t(i) + std::size_t(channel_count) * owned] += block_value;
                }
                stage("columns_complete_" + std::to_string(equations), solve_s, contraction_s,
                      elapsed_seconds_since(finecheck_started));
            }
        }
        if (equations != 2 * owned_columns.size())
            throw std::runtime_error("Weak matrix audit did not solve both signs for every owned column.");
        for (std::size_t owned = 0; owned < owned_columns.size(); ++owned)
            for (int i = 0; i < channel_count; ++i)
            {
                const auto value = partial_response[std::size_t(i) + std::size_t(channel_count) * owned];
                if (!std::isfinite(value.real()) || !std::isfinite(value.imag()))
                    throw std::runtime_error("Weak matrix audit signed response sum is nonfinite.");
                audit << "response " << i + 1 << ' ' << owned_columns[owned] + 1 << ' '
                      << value.real() << ' ' << value.imag() << '\n';
            }
        audit << "audit_complete yes\nresponse_equations " << equations << '\n';
        stage("complete");
        return;
    }
    auto channels = sample_sternheimer_abf_bloch_grid_channels(
        input.radials_by_type, input.atom_types, input.atom_positions, fine_grid.grid, qpoint, 2);
    if (channels.size() != 2) throw std::runtime_error("Fine-weak audit requires two real ABFS channels.");
    solve_sternheimer_abf_periodic_full_coulomb_in_place(channels, fine_grid.grid, qpoint, 0.0);
    std::array<Vector, 2> fine_vertices;
    std::array<Blocks::Vertices, 2> vertices;
    for (std::size_t j = 0; j < 2; ++j)
    {
        auto& g = fine_vertices[j];
        g.resize(source_values.size());
        // The real ABFS potentials are Ha; the weak Hamiltonian and RHS are Ry.
        for (std::size_t ir = 0; ir < g.size(); ++ir) g[ir] = 2.0 * channels[j].potential_r[ir] * source_values[ir];
        Vector gu(functions.size()), ge;
        for (std::size_t u = 0; u < functions.size(); ++u)
            gu[u] = sternheimer_fd_grid_dot(functions[u].values, g, fine_grid.volume_element);
        op->project(g, ge);
        vertices[j] = blocks->project_vertices(gu, ge);
    }
    stage("fine_vertices_ready");
    SternheimerRPA::SolverOptions options;
    options.max_iter = positive_int_from_env("ABACUS_STERNHEIMER_WEAK_AUDIT_MAX_ITER", 300);
    options.residual_tol = 1e-8;
    options.use_fd_spectral_preconditioner = false;
    audit << "coordinate_transform complement_inverse_sqrt\npreconditioner none\nmax_iter "
          << options.max_iter << "\nresidual_tolerance " << options.residual_tol << '\n';
    std::array<Complex, 4> response{};
    int equations = 0;
    for (int sign : {1, -1})
    {
        Blocks::Worker worker(blocks, [op](const Vector& x, Vector& result) { op->apply(x, result); },
                              epsilon, sign * 2.0 * omega_ha);
        for (std::size_t j = 0; j < 2; ++j)
        {
            const auto result = worker.solve(vertices[j], options);
            audit << "equation " << ++equations << " channel " << j + 1 << " sign " << sign
                  << " converged " << (result.converged ? "yes" : "no") << " iterations " << result.schur.iterations
                  << " full_residual " << result.relative_residual << " elapsed_s " << elapsed_seconds_since(started) << '\n';
            audit.flush();
            if (!result.converged) throw std::runtime_error("Fine-weak response audit equation did not converge.");
            const auto expansion = blocks->expand_coordinates(result.coefficients);
            Vector fine_response;
            op->lift(expansion.e, fine_response);
            for (std::size_t u = 0; u < functions.size(); ++u)
            {
#pragma omp parallel for schedule(static)
                for (std::size_t ir = 0; ir < fine_response.size(); ++ir)
                    fine_response[ir] += expansion.u[u] * functions[u].values[ir];
            }
            const Vector y(result.coefficients.begin(), result.coefficients.begin() + nvirtual);
            const Vector x(result.coefficients.begin() + nvirtual, result.coefficients.end());
            for (std::size_t i = 0; i < 2; ++i)
            {
                const Complex fine_value = 0.5 * sternheimer_fd_grid_dot(fine_vertices[i], fine_response, fine_grid.volume_element);
                const Complex block_value = 0.5 * (sternheimer_fd_grid_dot(vertices[i].f, y, 1.0)
                                                       + sternheimer_fd_grid_dot(vertices[i].w, x, 1.0));
                const double error = std::abs(fine_value - block_value) / std::max(1.0, std::abs(fine_value));
                audit << "vertex_check " << equations << ' ' << i + 1 << ' ' << error << '\n';
                if (!std::isfinite(error) || error > 1e-9)
                    throw std::runtime_error("Fine-weak left/right vertex consistency failed.");
                response[i + 2 * j] += fine_value;
            }
        }
    }
    for (int j = 0; j < 2; ++j)
        for (int i = 0; i < 2; ++i)
            audit << "response " << i + 1 << ' ' << j + 1 << ' ' << response[i + 2 * j].real()
                  << ' ' << response[i + 2 * j].imag() << '\n';
    double hermitian_error = 0.0;
    for (int j = 0; j < 2; ++j)
        for (int i = 0; i < 2; ++i)
            hermitian_error = std::max(hermitian_error, std::abs(response[i + 2 * j] - std::conj(response[j + 2 * i])));
    const auto off_diagonal = 0.5 * (response[1] + std::conj(response[2]));
    const double maximum_eigenvalue = 0.5 * (response[0].real() + response[3].real())
        + std::hypot(0.5 * (response[0].real() - response[3].real()), std::abs(off_diagonal));
    audit << "response_hermitian_error " << hermitian_error
          << "\nresponse_max_eigenvalue " << maximum_eigenvalue << '\n';
    audit.flush();
    if (!std::isfinite(maximum_eigenvalue) || hermitian_error > 1e-7 || maximum_eigenvalue > 1e-8)
        throw std::runtime_error("Fine-weak diagnostic response sign/Hermiticity failed.");
    audit << "audit_complete yes\nresponse_equations " << equations << '\n';
    stage("complete");
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

std::vector<std::vector<double>> transform_channel_potentials(
    SternheimerABFBuildData& abfs_data,
    const SternheimerFDHamiltonian::Grid& grid,
    const SternheimerCoulombWhitening& whitening)
{
    if (abfs_data.channels.size() != static_cast<std::size_t>(whitening.raw_dimension)
        || abfs_data.channels.empty())
    {
        throw std::invalid_argument("Sternheimer Coulomb transform does not match the raw ABFS channel count.");
    }
    return sample_sternheimer_abf_grid_channel_transform(abfs_data.radials_by_type,
                                                          abfs_data.atom_types,
                                                          abfs_data.atom_positions,
                                                          grid,
                                                          abfs_data.channels,
                                                          whitening.transform,
                                                          whitening.retained_rank);
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

void scale_potentials_in_place(std::vector<std::vector<double>>& potentials, const double factor)
{
    for (std::vector<double>& potential: potentials)
    {
        for (double& value: potential)
        {
            value *= factor;
        }
    }
}

void write_abfs_channel_diagnostic(const std::string& filename, const std::vector<SternheimerABFGridChannel>& channels)
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
        out << channel.channel_index << ' ' << channel.atom_index << ' ' << channel.atom_local_index << ' '
            << channel.type_index << ' ' << channel.angular_momentum << ' ' << channel.radial_index << ' '
            << channel.magnetic_index << ' ' << channel.label << ' ' << channel.max_abs << '\n';
    }
}

void write_coulomb_whitening_diagnostic(const std::string& filename,
                                         const SternheimerCoulombWhitening& whitening)
{
    std::ofstream out(filename.c_str(), std::ios::out | std::ios::trunc);
    if (!out)
    {
        throw std::runtime_error("Failed to open Sternheimer Coulomb whitening diagnostic file: " + filename);
    }
    out << std::setprecision(16);
    out << "# ABACUS Sternheimer global full-Coulomb whitening\n";
    out << "raw_auxiliary_dimension " << whitening.raw_dimension << '\n';
    out << "whitened_auxiliary_rank " << whitening.retained_rank << '\n';
    out << "discarded_auxiliary_rank " << whitening.discarded_rank << '\n';
    out << "relative_threshold " << whitening.relative_threshold << '\n';
    out << "largest_eigenvalue " << whitening.largest_eigenvalue << '\n';
    out << "smallest_retained_eigenvalue " << whitening.smallest_retained_eigenvalue << '\n';
    out << "max_orthonormality_error " << whitening.max_orthonormality_error << '\n';
    out << "transform_sha256 " << hash_coulomb_whitening_transform(whitening) << '\n';
    out << "# index eigenvalue retained\n";
    const int first_retained = whitening.discarded_rank;
    for (int index = 0; index != whitening.raw_dimension; ++index)
    {
        out << index << ' ' << whitening.eigenvalues[static_cast<std::size_t>(index)] << ' '
            << (index >= first_retained ? "yes" : "no") << '\n';
    }
}

SternheimerRPA::Chi0V1Metadata make_chi0_v1_metadata(const UnitCell& ucell,
                                                     const std::vector<Channel>& channels,
                                                     const int iq,
                                                     const int ifrequency,
                                                     const double omega_ha,
                                                     const double weight_ha,
                                                     const int output_atom_count = -1)
{
    SternheimerRPA::Chi0V1Metadata metadata;
    metadata.iq = iq;
    metadata.ifrequency = ifrequency;
    metadata.omega = omega_ha;
    metadata.weight = weight_ha;
    const int atom_count = output_atom_count > 0 ? output_atom_count : ucell.nat;
    metadata.atom_naux.assign(static_cast<std::size_t>(atom_count), 0);
    for (const Channel& channel: channels)
    {
        if (channel.atom_index < 0 || channel.atom_index >= atom_count)
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
                                    SternheimerFDHamiltonian::Vector(static_cast<std::size_t>(grid_size),
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

void run_sternheimer_periodic_lcao_chi0_output(const elecstate::Potential& potential,
                                               const ModulePW::PW_Basis& pw_basis,
                                               const UnitCell& ucell,
                                               const elecstate::ElecState& elec_state,
                                               const LCAO_Orbitals& orbitals,
                                               const std::vector<SternheimerLCAOOccupiedKPoint>& occupied_kpoints,
                                               const std::array<int, 3>& kmesh,
                                               std::ofstream& out,
                                               const bool use_frequency_mpi,
                                               const int kpoint_groups,
                                               const std::chrono::steady_clock::time_point& chi0_start_time,
                                               const SternheimerOrbitalSet* reusable_rpa_abfs)
{
    if (PARAM.inp.nspin != 1)
    {
        throw std::runtime_error("The first periodic Sternheimer driver supports only nspin=1 insulators.");
    }
    if (PARAM.inp.symmetry != "-1" && PARAM.inp.symmetry != "1")
    {
        throw std::runtime_error("Periodic Sternheimer output requires symmetry=-1 or symmetry=1.");
    }

    const char* supercell_translation_sum_raw = std::getenv(kSupercellTranslationSumEnv);
    const bool use_supercell_translation_sum
        = supercell_translation_sum_raw != nullptr && supercell_translation_sum_raw[0] != '\0';
    const bool full_supercell_response
        = use_supercell_translation_sum && env_is_true(kSupercellFullResponseEnv);
    if (env_is_true(kSupercellFullResponseEnv) && !use_supercell_translation_sum)
    {
        throw std::runtime_error(
            "Full supercell Sternheimer response requires a translation-sum specification.");
    }
    if (full_supercell_response && PARAM.inp.symmetry != "-1")
    {
        throw std::runtime_error(
            "The first full supercell Sternheimer response gate requires symmetry=-1.");
    }

    SternheimerSupercellTranslationSum supercell_translation_sum;
    std::vector<SternheimerLCAOOccupiedKPoint> recovered_supercell_kpoints;
    int supercell_sector_dimension = 0;
    double supercell_sector_kweight = 0.0;
    double supercell_sector_max_orthonormality_error = 0.0;
    double supercell_sector_max_full_space_residual = 0.0;
    if (use_supercell_translation_sum)
    {
        supercell_translation_sum
            = parse_sternheimer_supercell_translation_sum(supercell_translation_sum_raw);
        const int expected_atoms
            = sternheimer_supercell_primitive_cell_count(supercell_translation_sum)
              * supercell_translation_sum.atoms_per_primitive;
        constexpr double gamma_tolerance = 1.0e-10;
        if (ucell.nat != expected_atoms || occupied_kpoints.size() != 1
            || std::any_of(occupied_kpoints.front().kpoint.begin(),
                           occupied_kpoints.front().kpoint.end(),
                           [gamma_tolerance](const double coordinate) {
                               return std::abs(coordinate) > gamma_tolerance;
                           }))
        {
            throw std::runtime_error(
                "Supercell translation response requires one Gamma record and a matching translation-major atom list.");
        }
        if (full_supercell_response)
        {
            const auto& gamma_record = occupied_kpoints.front();
            std::vector<double> complete_eigenvalues = gamma_record.eigenvalues;
            complete_eigenvalues.insert(complete_eigenvalues.end(),
                                        gamma_record.unoccupied_eigenvalues.begin(),
                                        gamma_record.unoccupied_eigenvalues.end());
            std::vector<std::vector<SternheimerRPA::Complex>> complete_coefficients
                = gamma_record.coefficients;
            complete_coefficients.insert(complete_coefficients.end(),
                                         gamma_record.unoccupied_coefficients.begin(),
                                         gamma_record.unoccupied_coefficients.end());
            const auto sectors = recover_all_sternheimer_supercell_sectors(
                complete_eigenvalues,
                complete_coefficients,
                supercell_translation_sum.repeats);
            for (const auto& sector: sectors)
            {
                supercell_sector_dimension
                    = std::max(supercell_sector_dimension,
                               static_cast<int>(sector.sector.eigenvalues.size()));
                supercell_sector_max_orthonormality_error
                    = std::max(supercell_sector_max_orthonormality_error,
                               sector.sector.max_orthonormality_error);
                supercell_sector_max_full_space_residual
                    = std::max(supercell_sector_max_full_space_residual,
                               sector.sector.max_full_space_residual);
            }
            recovered_supercell_kpoints
                = build_sternheimer_supercell_full_kpoint_records(gamma_record, sectors);
            supercell_sector_kweight = recovered_supercell_kpoints.front().kweight;
        }
    }

    const auto& response_kpoints
        = full_supercell_response ? recovered_supercell_kpoints : occupied_kpoints;
    const int supercell_primitive_cell_count
        = use_supercell_translation_sum
              ? sternheimer_supercell_primitive_cell_count(supercell_translation_sum)
              : 1;
    const double response_matrix_scale
        = sternheimer_supercell_response_matrix_scale(full_supercell_response,
                                                       supercell_primitive_cell_count);
    const std::array<int, 3> response_kmesh
        = full_supercell_response ? supercell_translation_sum.repeats : kmesh;
    const int response_q_index
        = full_supercell_response
              ? sternheimer_find_kpoint_one_based(response_kpoints,
                                                   supercell_translation_sum.primitive_qpoint)
              : PARAM.inp.sternheimer_q_index;
    const SternheimerPeriodicResponsePlan response_plan
        = build_sternheimer_periodic_response_plan(
            response_kpoints,
            response_q_index,
            use_supercell_translation_sum && !full_supercell_response);
    if (PARAM.inp.sternheimer_q_index <= 0)
    {
        throw std::runtime_error("Internal error: the periodic Sternheimer path requires a positive q index.");
    }
    validate_sternheimer_periodic_kmesh(response_kmesh,
                                        static_cast<int>(response_kpoints.size()));
    const bool use_delta_sternheimer = PARAM.inp.sternheimer_delta;
    const std::array<int, 3> requested_response_dimensions = {PARAM.inp.sternheimer_response_nx,
                                                              PARAM.inp.sternheimer_response_ny,
                                                              PARAM.inp.sternheimer_response_nz};
    const bool explicit_response_dimensions_requested
        = std::any_of(requested_response_dimensions.begin(),
                      requested_response_dimensions.end(),
                      [](const int value) { return value != 0; });
    if ((PARAM.inp.sternheimer_response_ecutwfc > 0.0 || explicit_response_dimensions_requested)
        && !use_delta_sternheimer)
    {
        throw std::runtime_error("An independent Sternheimer response grid requires sternheimer_delta=true.");
    }
    SternheimerResponseGrid response_grid = make_sternheimer_response_grid(pw_basis,
                                                                            PARAM.inp.ecutwfc,
                                                                            PARAM.inp.sternheimer_response_ecutwfc,
                                                                            PARAM.inp.fft_mode,
                                                                            requested_response_dimensions,
                                                                            PARAM.inp.sternheimer_fd_order);
    const char* response_grid_source = sternheimer_response_grid_source_name(response_grid.source);
    const double response_grid_ecutwfc
        = response_grid.source == SternheimerResponseGridSource::Cutoff
              ? PARAM.inp.sternheimer_response_ecutwfc
              : (response_grid.source == SternheimerResponseGridSource::Pbe ? PARAM.inp.ecutwfc : 0.0);
    const ModulePW::PW_Basis& response_pw_basis = *response_grid.basis;
    const bool use_kpoint_mpi = kpoint_groups > 1;
    const bool use_parallel_grid_mpi = use_frequency_mpi || use_kpoint_mpi;
    const SternheimerABACUSFDGridData grid_data
        = use_parallel_grid_mpi ? make_sternheimer_fd_full_grid(response_pw_basis)
                                : make_sternheimer_fd_grid(response_pw_basis);
    const bool use_symmetry_partial_response = PARAM.inp.symmetry == "1";
    const bool write_kresolved_diagnostic
        = !use_symmetry_partial_response && env_is_true(kKResolvedDiagnosticEnv);
    const bool write_partial_kresolved = use_symmetry_partial_response || write_kresolved_diagnostic;
    std::vector<int> canonical_q_indices
        = sternheimer_canonical_q_indices_one_based(response_kpoints);
    std::vector<SternheimerFixedQKOrbit> fixed_q_orbits;
    std::vector<SternheimerFixedQKRoute> fixed_q_routes;
    std::vector<SternheimerQStarRoute> qstar_routes;
    int fixed_q_little_group_order = 1;
    int fixed_q_discrete_spatial_order = 1;
    std::vector<bool> fixed_q_representative(response_kpoints.size(), true);
    if (use_symmetry_partial_response)
    {
        std::vector<SternheimerFDReducedRotation> reduced_rotations;
        reduced_rotations.reserve(static_cast<std::size_t>(ucell.symm.nrotk));
        for (int isym = 0; isym != ucell.symm.nrotk; ++isym)
        {
            const auto& rotation = ucell.symm.gmatrix[isym];
            reduced_rotations.push_back({{{rotation.e11, rotation.e12, rotation.e13},
                                          {rotation.e21, rotation.e22, rotation.e23},
                                          {rotation.e31, rotation.e32, rotation.e33}}});
        }
        const auto discrete_spatial_operations
            = sternheimer_fd_second_order_stencil_symmetry_indices(grid_data.grid,
                                                                    reduced_rotations);
        fixed_q_discrete_spatial_order = static_cast<int>(discrete_spatial_operations.size());
        const auto qstar_permutations = build_sternheimer_discrete_qstar_permutations(
            ucell, response_kpoints, discrete_spatial_operations);
        qstar_routes = build_sternheimer_qstar_routes_from_permutations(
            static_cast<int>(response_kpoints.size()), qstar_permutations);
        canonical_q_indices.clear();
        for (const auto& route: qstar_routes)
        {
            if (route.representative_iq == route.member_iq)
            {
                canonical_q_indices.push_back(route.representative_iq);
            }
        }
        const auto selected_q_route = std::find_if(
            qstar_routes.begin(), qstar_routes.end(), [&](const auto& route) {
                return route.member_iq == response_plan.iq;
            });
        if (selected_q_route == qstar_routes.end()
            || selected_q_route->representative_iq != response_plan.iq)
        {
            throw std::runtime_error(
                "Sternheimer symmetry output requires sternheimer_q_index to select a discrete q-star representative.");
        }
        const auto permutations = build_sternheimer_fixed_q_little_group_permutations(
            ucell, response_kpoints, response_plan.qpoint, discrete_spatial_operations);
        fixed_q_little_group_order = static_cast<int>(permutations.size());
        std::vector<std::vector<int>> index_permutations;
        index_permutations.reserve(permutations.size());
        for (const auto& permutation: permutations)
        {
            index_permutations.push_back(permutation.mapped_index_by_full_k);
        }
        fixed_q_orbits = build_sternheimer_fixed_q_k_orbits_from_permutations(
            static_cast<int>(response_kpoints.size()), index_permutations);
        std::fill(fixed_q_representative.begin(), fixed_q_representative.end(), false);
        for (const auto& orbit: fixed_q_orbits)
        {
            fixed_q_representative[static_cast<std::size_t>(orbit.representative_ik_full)] = true;
            for (const int member_ik_full: orbit.members)
            {
                const auto route = std::find_if(
                    permutations.begin(),
                    permutations.end(),
                    [&](const auto& operation) {
                        return operation.mapped_index_by_full_k[
                                   static_cast<std::size_t>(member_ik_full)]
                               == orbit.representative_ik_full;
                    });
                if (route == permutations.end())
                {
                    throw std::runtime_error(
                        "A Sternheimer fixed-q orbit member has no inverse route to its representative.");
                }
                fixed_q_routes.push_back(
                    {response_plan.iq,
                     orbit.representative_ik_full,
                     member_ik_full,
                     route->spatial_isym,
                     route->time_reversal,
                     route->fold_G_by_full_k[static_cast<std::size_t>(member_ik_full)]});
            }
        }
    }
    const SternheimerKPointSchedule kpoint_schedule(fixed_q_representative,
                                                     use_kpoint_mpi ? kpoint_groups : 1);
    constexpr double q_tolerance = 1.0e-10;
    const bool gamma_qpoint
        = std::all_of(response_plan.qpoint.begin(), response_plan.qpoint.end(), [q_tolerance](const double coordinate) {
              return std::abs(coordinate) <= q_tolerance;
          });
    const double massidda_chi = gamma_qpoint
                                    ? Singular_Value::cal_massidda(ucell, response_kmesh, 2, 1.0, 5, 1.0e-4)
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

    const double solver_tolerance = positive_double_from_env(kSolverToleranceEnv, 1.0e-8);
    const int solver_max_iter = positive_int_from_env(kSolverMaxIterEnv, 300);
    const double pca_threshold = nonnegative_double_from_env(kPCAThresholdEnv, PARAM.inp.exx_pca_threshold);
    const double threshold_scale = std::max(1.0, std::abs(PARAM.inp.exx_pca_threshold));
    const bool reuse_rpa_abfs = reusable_rpa_abfs != nullptr && !reusable_rpa_abfs->empty()
                                && std::abs(pca_threshold - PARAM.inp.exx_pca_threshold) <= 1.0e-14 * threshold_scale;
    const SternheimerOrbitalSet* periodic_abfs_orbitals = reuse_rpa_abfs ? reusable_rpa_abfs : nullptr;
    const double ccp_rmesh_times = positive_double_from_env(kCCPRmeshTimesEnv, PARAM.inp.rpa_ccp_rmesh_times);
    const int max_channels = positive_int_from_env(kChannelsEnv, -1);
    const int max_bands = positive_int_from_env(kBandsEnv, -1);
    const int channel_threads = positive_int_from_env(kChannelThreadsEnv, 0);
    const int nfreq = PARAM.inp.sternheimer_nfreq;
    const bool use_channel_mpi = PARAM.inp.sternheimer_channel_mpi;
    const bool use_nested_response_mpi = use_frequency_mpi && use_kpoint_mpi;
    const SternheimerNestedMPIReplicaLayout nested_replica_layout
        = use_nested_response_mpi
              ? sternheimer_nested_mpi_replica_layout(kpoint_groups,
                                                       nfreq,
                                                       GlobalV::NPROC,
                                                       GlobalV::MY_RANK,
                                                       use_channel_mpi)
              : SternheimerNestedMPIReplicaLayout();
    const int response_replicas
        = use_nested_response_mpi
              ? nested_replica_layout.replicas_per_slot
              : (use_frequency_mpi && use_channel_mpi ? GlobalV::NPROC / nfreq : 1);
    const int local_channel_replica
        = use_nested_response_mpi
              ? nested_replica_layout.local_replica
              : (use_frequency_mpi && use_channel_mpi ? GlobalV::MY_RANK % response_replicas : 0);
    const int local_kpoint_group
        = use_kpoint_mpi
              ? (use_nested_response_mpi
                     ? nested_replica_layout.local_response_slot / nfreq
                     : GlobalV::MY_RANK)
              : 0;
    const int local_frequency_slot
        = use_nested_response_mpi ? nested_replica_layout.local_response_slot % nfreq : 0;
    const int default_frequency_rank_shift = use_frequency_mpi && GlobalV::NPROC > 1 ? 1 : 0;
    const int frequency_rank_shift = int_from_env(kFrequencyRankShiftEnv, default_frequency_rank_shift);
    const SternheimerDeltaABlockMode delta_a_block_mode = delta_a_block_mode_from_env();
    const bool write_delta_components = use_delta_sternheimer && env_is_true(kDeltaComponentDiagnosticEnv);
    const bool write_lcao_sos = env_is_true(kLCAOSOSDiagnosticEnv);
    const char* wavefunction_diagnostic_raw = std::getenv(kWavefunctionDiagnosticEnv);
    const bool write_wavefunction_diagnostic
        = wavefunction_diagnostic_raw != nullptr && wavefunction_diagnostic_raw[0] != '\0';
    SternheimerWavefunctionDiagnostic::Configuration wavefunction_diagnostic_config;
    if (write_wavefunction_diagnostic)
    {
        wavefunction_diagnostic_config
            = SternheimerWavefunctionDiagnostic::parse_configuration(wavefunction_diagnostic_raw);
    }
    const bool bands_are_truncated
        = max_bands > 0
          && std::any_of(response_kpoints.begin(),
                         response_kpoints.end(),
                         [max_bands](const SternheimerLCAOOccupiedKPoint& record) {
                             return static_cast<int>(record.coefficients.size()) > max_bands;
                         });
    const bool write_truncated_v1_diagnostic
        = bands_are_truncated && env_is_true(kTruncatedV1DiagnosticEnv);
    const bool write_periodic_v1
        = sternheimer_write_periodic_v1(use_supercell_translation_sum,
                                        bands_are_truncated,
                                        full_supercell_response,
                                        write_truncated_v1_diagnostic);
    validate_sternheimer_periodic_output_mode(write_periodic_v1, write_partial_kresolved);
    if (use_symmetry_partial_response
        && (write_delta_components || write_lcao_sos || write_wavefunction_diagnostic))
    {
        throw std::runtime_error(
            "Symmetry-reduced Sternheimer partial output does not yet support component diagnostics.");
    }
    if (write_lcao_sos)
    {
        for (const SternheimerLCAOOccupiedKPoint& record: response_kpoints)
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
        = full_supercell_response
              ? transition_window_from_recovered_kpoints(response_kpoints)
              : transition_window_from_all_kpoints(elec_state, response_kpoints);
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

    const auto restrict_response_field = [&](const std::vector<double>& fine_values) {
        return response_grid.source == SternheimerResponseGridSource::Explicit
                   ? restrict_sternheimer_real_field_rectangular(*response_grid.serial_fine_basis,
                                                                 *response_grid.serial_response_basis,
                                                                 fine_values)
                   : restrict_sternheimer_real_field(*response_grid.serial_fine_basis,
                                                     *response_grid.serial_response_basis,
                                                     fine_values);
    };
    const std::vector<double> independent_response_full_potential
        = response_grid.independent
              ? restrict_response_field(copy_sternheimer_full_local_potential(potential, pw_basis, 0))
              : std::vector<double>();
    const std::vector<double> kpoint_parallel_full_potential
        = use_kpoint_mpi
              ? (response_grid.independent
                     ? independent_response_full_potential
                     : copy_sternheimer_full_local_potential(potential, pw_basis, 0))
              : std::vector<double>();
    const std::vector<double> diagnostic_fixed_local_potential
        = write_wavefunction_diagnostic
              ? (response_grid.independent
                     ? restrict_response_field(copy_sternheimer_full_fixed_local_potential(potential, pw_basis))
                     : (use_kpoint_mpi
                            ? copy_sternheimer_full_fixed_local_potential(potential, pw_basis)
                            : copy_sternheimer_fixed_local_potential(potential, pw_basis)))
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
                               "pbe_grid=" + std::to_string(pw_basis.nx) + "x"
                                   + std::to_string(pw_basis.ny) + "x" + std::to_string(pw_basis.nz)
                                   + ",response_grid=" + std::to_string(grid_data.grid.nx) + "x"
                                   + std::to_string(grid_data.grid.ny) + "x"
                                   + std::to_string(grid_data.grid.nz) + ",grid_size="
                                   + std::to_string(grid_data.grid.size()));
    if (env_is_true("ABACUS_STERNHEIMER_WEAK_RESPONSE_AUDIT")
        || env_is_true("ABACUS_STERNHEIMER_WEAK_MATRIX_AUDIT"))
    {
        if (GlobalV::NPROC != 1 || use_supercell_translation_sum || !use_delta_sternheimer
            || PARAM.inp.sternheimer_delta_max_states != 0 || gamma_qpoint)
            throw std::runtime_error("Fine-weak response audit requires one rank, full states, and nonzero q.");
        const auto& pair = response_plan.kq_pairs.at(0);
        const auto& source = response_kpoints.at(response_plan.record_index_by_global_k.at(pair.source_index));
        const auto& target = response_kpoints.at(response_plan.record_index_by_global_k.at(pair.target_index));
        run_sternheimer_weak_response_audit(
            ucell, orbitals, make_sternheimer_fd_full_grid(pw_basis), grid_data,
            copy_sternheimer_full_local_potential(potential, pw_basis, 0), source, target,
            response_plan.qpoint, frequency_grid.omega_ha.at(0), pca_threshold, periodic_abfs_orbitals,
            pair.source_index, pair.target_index);
        out << "status diagnostic_only\nphysical_result no\nweak_response_audit_complete yes\n";
        return;
    }
    if (env_is_true("ABACUS_STERNHEIMER_SPECTRUM_AUDIT")
        || env_is_true("ABACUS_STERNHEIMER_GALERKIN_AUDIT"))
    {
        if (GlobalV::NPROC != 1 || use_supercell_translation_sum || !use_delta_sternheimer)
        {
            throw std::runtime_error("Spectrum-only audit requires one MPI rank and the direct periodic Delta route.");
        }
        const auto& pair = response_plan.kq_pairs.at(0);
        const auto& source = response_kpoints.at(response_plan.record_index_by_global_k.at(pair.source_index));
        const auto& target = response_kpoints.at(response_plan.record_index_by_global_k.at(pair.target_index));
        if (env_is_true("ABACUS_STERNHEIMER_GALERKIN_AUDIT"))
        {
            run_sternheimer_galerkin_spectrum_audit(
                ucell, orbitals, make_sternheimer_fd_full_grid(pw_basis), grid_data,
                copy_sternheimer_full_local_potential(potential, pw_basis, 0), source, target,
                pair.source_index, pair.target_index);
        }
        else
        {
            const auto audit_potential = response_grid.independent
                ? independent_response_full_potential
                : copy_sternheimer_full_local_potential(potential, pw_basis, 0);
            run_sternheimer_spectrum_audit(ucell, orbitals, grid_data, audit_potential, source, target,
                                           pair.source_index, pair.target_index);
        }
        out << "status diagnostic_only\n" << "physical_result no\n"
            << "response_equations 0\n" << "spectrum_audit_complete yes\n";
        return;
    }
    append_chi0_progress_event("abfs_source_ready",
                               0,
                               -1,
                               -1,
                               -1,
                               0,
                               nullptr,
                               -1.0,
                               elapsed_seconds_since(chi0_start_time),
                               std::string("source=") + (reuse_rpa_abfs ? "rpa_handoff" : "rebuilt"));
    SternheimerPeriodicABFGridData periodic_abfs;
    if (full_supercell_response)
    {
        if (use_kpoint_mpi)
        {
            const SternheimerABFSInput collective_abfs_input
                = build_abfs_density_input(ucell, pca_threshold, periodic_abfs_orbitals);
            if (GlobalV::MY_RANK == 0)
            {
                periodic_abfs = build_supercell_translation_full_response_abfs_from_input(
                    collective_abfs_input,
                    grid_data.grid,
                    supercell_translation_sum,
                    max_channels);
            }
        }
        else
        {
            periodic_abfs = build_supercell_translation_full_response_abfs(ucell,
                                                                           grid_data.grid,
                                                                           supercell_translation_sum,
                                                                           max_channels,
                                                                           pca_threshold,
                                                                           periodic_abfs_orbitals);
        }
        broadcast_periodic_abfs(periodic_abfs, grid_data.grid.size(), use_kpoint_mpi);
    }
    else
    {
        periodic_abfs = build_abfs_full_coulomb_bloch_grid_channels(ucell,
                                                                    grid_data.grid,
                                                                    response_plan.qpoint,
                                                                    gamma_inverse_k2,
                                                                    max_channels,
                                                                    pca_threshold,
                                                                    periodic_abfs_orbitals,
                                                                    chi0_start_time);
    }
    if (use_supercell_translation_sum && !full_supercell_response)
    {
        if (!periodic_abfs.densities.empty())
        {
            periodic_abfs.densities = {
                combine_sternheimer_supercell_translation_channel(periodic_abfs.densities, supercell_translation_sum)};
        }
        periodic_abfs.potentials = {combine_sternheimer_supercell_translation_channel(
            periodic_abfs.potentials, supercell_translation_sum)};
    }
    if (periodic_abfs.potentials.empty())
    {
        throw std::runtime_error("No periodic ABFS full-Coulomb perturbation channels were generated.");
    }
    const int num_channels = static_cast<int>(periodic_abfs.potentials.size());
    const int channel_batch_width = use_delta_sternheimer ? sternheimer_channel_batch_width() : 1;
    const SternheimerMemorySnapshot channel_memory = detect_sternheimer_memory_snapshot();
    SternheimerChannelWorkerPlan channel_worker_plan
        = plan_sternheimer_channel_workers(num_channels,
                                            sternheimer_channel_openmp_threads(),
                                            grid_data.grid.size(),
                                            channel_threads,
                                            channel_memory,
                                            channel_batch_width);
    append_chi0_progress_event(
        "channel_workers_ready",
        0,
        -1,
        -1,
        -1,
        0,
        nullptr,
        -1.0,
        elapsed_seconds_since(chi0_start_time),
        "global_channels=" + std::to_string(num_channels)
            + " local_channels=" + std::to_string(num_channels) + " "
            + format_sternheimer_channel_worker_diagnostic(channel_memory,
                                                           channel_worker_plan,
                                                           grid_data.grid.size(),
                                                           channel_threads));
    const int output_atom_count
        = full_supercell_response ? supercell_translation_sum.atoms_per_primitive : ucell.nat;
    double gamma_projection_relative_error = -1.0;
    const bool evaluate_gamma_projection
        = gamma_qpoint && !use_supercell_translation_sum && !periodic_abfs.densities.empty();
    if (evaluate_gamma_projection)
    {
        const auto target = SternheimerRPA::read_coulomb_v1_files(
            find_coulomb_v1_rank_files(response_plan.iq, GlobalV::NPROC));
        if (target.iq != response_plan.iq
            || target.atom_naux != atom_auxiliary_sizes(periodic_abfs.potentials, ucell.nat))
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
            if (gamma_qpoint && !use_supercell_translation_sum)
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
            out << "pbe_ecutwfc_Ry " << PARAM.inp.ecutwfc << '\n';
            out << "sternheimer_response_ecutwfc_requested_Ry "
                << PARAM.inp.sternheimer_response_ecutwfc << '\n';
            out << "sternheimer_response_ecutwfc_Ry " << response_grid_ecutwfc << '\n';
            out << "sternheimer_response_grid_requested " << response_grid.requested_dimensions[0] << ' '
                << response_grid.requested_dimensions[1] << ' ' << response_grid.requested_dimensions[2] << '\n';
            out << "sternheimer_response_grid_source " << response_grid_source << '\n';
            out << "pbe_grid " << pw_basis.nx << ' ' << pw_basis.ny << ' ' << pw_basis.nz << " size "
                << pw_basis.nxyz << '\n';
            out << "response_grid " << grid_data.grid.nx << ' ' << grid_data.grid.ny << ' '
                << grid_data.grid.nz << " size " << grid_data.grid.size() << '\n';
            out << "abfs_channels " << num_channels << '\n';
            out << "perturbation_coulomb_kernel full_periodic_poisson\n";
            out << "periodic_kmesh " << response_kmesh[0] << ' ' << response_kmesh[1] << ' '
                << response_kmesh[2] << '\n';
            out << "periodic_gamma_massidda_chi " << massidda_chi << '\n';
            out << "periodic_gamma_coulomb_projection "
                << (use_supercell_translation_sum ? "skipped_supercell_translation_sum"
                                                  : (evaluate_gamma_projection ? "diagnostic_only_physical_poisson"
                                                                               : "skipped_large_auxiliary_dimension"))
                << '\n';
            out << "periodic_gamma_projection_relative_error " << gamma_projection_relative_error << '\n';
            out << "periodic_gamma_limit constant_mode_only_no_headwing\n";
            out << "perturbation_ccp_rmesh_times_used no\n";
        }
    }
    if (sternheimer_abfs_diag_only_enabled())
    {
        return;
    }

    const std::vector<SternheimerRPA::AuxiliaryChannel> auxiliary_channels
        = make_chi0_auxiliary_channels(channels);

    const auto response_owner_rank = [&](const int source_kpoint_index,
                                         const int ifrequency) {
        if (use_nested_response_mpi)
        {
            const auto assignment = kpoint_schedule.assignment(
                source_kpoint_index,
                ifrequency,
                nfreq,
                kpoint_groups * nfreq,
                frequency_rank_shift);
            return assignment.owner_rank * response_replicas;
        }
        if (use_kpoint_mpi)
        {
            return kpoint_schedule.owner_group(source_kpoint_index);
        }
        if (use_frequency_mpi && use_channel_mpi)
        {
            return SternheimerRPA::frequency_mpi_assignment(ifrequency,
                                                             nfreq,
                                                             GlobalV::NPROC,
                                                             GlobalV::MY_RANK,
                                                             frequency_rank_shift,
                                                             true)
                .frequency_leader_rank;
        }
        return SternheimerRPA::frequency_owner_rank(ifrequency,
                                                    GlobalV::NPROC,
                                                    frequency_rank_shift);
    };
    const auto response_owned_by_rank = [&](const int source_kpoint_index,
                                             const int ifrequency) {
        const int leader = response_owner_rank(source_kpoint_index, ifrequency);
        return GlobalV::MY_RANK >= leader && GlobalV::MY_RANK < leader + response_replicas;
    };

#ifdef __MPI
    MPI_Comm response_group_comm = MPI_COMM_NULL;
    if (response_replicas > 1)
    {
        const int response_group_color = GlobalV::MY_RANK / response_replicas;
        MPI_Comm_split(MPI_COMM_WORLD,
                       response_group_color,
                       local_channel_replica,
                       &response_group_comm);
    }
#endif

    std::vector<int> frequency_owners(static_cast<std::size_t>(nfreq), 0);
    std::vector<std::vector<SternheimerRPA::Complex>> chi0_branches(static_cast<std::size_t>(nfreq));
    std::vector<std::vector<SternheimerRPA::Complex>> delta_sos_branches(static_cast<std::size_t>(nfreq));
    std::vector<std::vector<SternheimerRPA::Complex>> delta_pulay_branches(static_cast<std::size_t>(nfreq));
    std::vector<std::vector<SternheimerRPA::Complex>> delta_out_branches(static_cast<std::size_t>(nfreq));
    std::vector<std::vector<SternheimerRPA::Complex>> lcao_sos_branches(static_cast<std::size_t>(nfreq));
    for (int ifrequency = 0; ifrequency != nfreq; ++ifrequency)
    {
        const int owner_rank = use_frequency_mpi && !use_kpoint_mpi
                                   ? response_owner_rank(0, ifrequency)
                                   : 0;
        frequency_owners[static_cast<std::size_t>(ifrequency)] = owner_rank;
        if (use_kpoint_mpi || response_owned_by_rank(0, ifrequency))
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
    solver_options.use_fd_spectral_preconditioner
        = sternheimer_environment_flag(kSpectralPreconditionerEnv, true);
    solver_options.fd_spectral_preconditioner_regularization
        = nonnegative_double_from_env(kSpectralPreconditionerRegularizationEnv, 0.0);
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
    int local_wavefunction_diagnostic_count = 0;

    const std::vector<std::size_t> owned_pair_indices
        = kpoint_schedule.owned_pair_indices(response_plan, use_kpoint_mpi ? local_kpoint_group : 0);
    std::vector<SternheimerPartialResponseRecord> local_partial_records;
    for (const std::size_t pair_index: owned_pair_indices)
    {
        const SternheimerKQPair& pair = response_plan.kq_pairs[pair_index];
        const int source_record_index
            = response_plan.record_index_by_global_k[static_cast<std::size_t>(pair.source_index)];
        const int target_record_index
            = response_plan.record_index_by_global_k[static_cast<std::size_t>(pair.target_index)];
        const SternheimerLCAOOccupiedKPoint& full_source_record
            = response_kpoints[static_cast<std::size_t>(source_record_index)];
        const SternheimerLCAOOccupiedKPoint& full_target_record
            = response_kpoints[static_cast<std::size_t>(target_record_index)];
        SternheimerLCAOOccupiedKPoint source_sector_record;
        SternheimerLCAOOccupiedKPoint target_sector_record;
        const SternheimerLCAOOccupiedKPoint* source_record_pointer = &full_source_record;
        const SternheimerLCAOOccupiedKPoint* target_record_pointer = &full_target_record;
        if (use_supercell_translation_sum && !full_supercell_response)
        {
            const int cell_count
                = sternheimer_supercell_primitive_cell_count(supercell_translation_sum);
            if (full_source_record.coefficients.size() % static_cast<std::size_t>(cell_count) != 0
                || full_target_record.coefficients.size() % static_cast<std::size_t>(cell_count) != 0)
            {
                throw std::runtime_error(
                    "Supercell translation-sector occupied bands do not divide evenly over primitive cells.");
            }
            const int source_occupied_count
                = static_cast<int>(full_source_record.coefficients.size()) / cell_count;
            const int target_occupied_count
                = static_cast<int>(full_target_record.coefficients.size()) / cell_count;
            auto recover_record
                = [&](const SternheimerLCAOOccupiedKPoint& full_record,
                      const SternheimerReducedKPoint& primitive_kpoint,
                      const int occupied_count,
                      SternheimerLCAOOccupiedKPoint& sector_record) {
                      std::vector<double> complete_eigenvalues = full_record.eigenvalues;
                      complete_eigenvalues.insert(complete_eigenvalues.end(),
                                                  full_record.unoccupied_eigenvalues.begin(),
                                                  full_record.unoccupied_eigenvalues.end());
                      std::vector<std::vector<SternheimerRPA::Complex>> complete_coefficients
                          = full_record.coefficients;
                      complete_coefficients.insert(complete_coefficients.end(),
                                                   full_record.unoccupied_coefficients.begin(),
                                                   full_record.unoccupied_coefficients.end());
                      const SternheimerSupercellSector sector
                          = recover_sternheimer_supercell_sector(complete_eigenvalues,
                                                                 complete_coefficients,
                                                                 supercell_translation_sum.repeats,
                                                                 primitive_kpoint);
                      supercell_sector_dimension = static_cast<int>(sector.eigenvalues.size());
                      supercell_sector_max_orthonormality_error
                          = std::max(supercell_sector_max_orthonormality_error,
                                     sector.max_orthonormality_error);
                      supercell_sector_max_full_space_residual
                          = std::max(supercell_sector_max_full_space_residual,
                                     sector.max_full_space_residual);
                      if (occupied_count <= 0
                          || occupied_count >= static_cast<int>(sector.eigenvalues.size())
                          || full_record.occupations.empty())
                      {
                          throw std::runtime_error(
                              "Supercell translation-sector recovery found an invalid occupied dimension.");
                      }
                      validate_sternheimer_supercell_sector_occupations(
                          full_record.occupations,
                          static_cast<int>(full_record.coefficients.size()));
                      sector_record = full_record;
                      sector_record.kweight
                          = sternheimer_supercell_sector_kweight(full_record.kweight, cell_count);
                      supercell_sector_kweight = sector_record.kweight;
                      sector_record.eigenvalues.assign(sector.eigenvalues.begin(),
                                                       sector.eigenvalues.begin() + occupied_count);
                      sector_record.coefficients.assign(sector.coefficients.begin(),
                                                        sector.coefficients.begin() + occupied_count);
                      sector_record.occupations.assign(
                          static_cast<std::size_t>(occupied_count),
                          full_record.occupations.front());
                      sector_record.unoccupied_eigenvalues.assign(
                          sector.eigenvalues.begin() + occupied_count,
                          sector.eigenvalues.end());
                      sector_record.unoccupied_coefficients.assign(
                          sector.coefficients.begin() + occupied_count,
                          sector.coefficients.end());
                      append_chi0_progress_event(
                          "supercell_sector_ready",
                          0,
                          -1,
                          -1,
                          -1,
                          solved_equations,
                          nullptr,
                          -1.0,
                          elapsed_seconds_since(chi0_start_time),
                          "primitive_k=" + std::to_string(primitive_kpoint[0]) + ":"
                              + std::to_string(primitive_kpoint[1]) + ":"
                              + std::to_string(primitive_kpoint[2])
                              + " dimension=" + std::to_string(sector.eigenvalues.size())
                              + " orth_error="
                              + precise_double_string(sector.max_orthonormality_error)
                              + " residual="
                              + precise_double_string(sector.max_full_space_residual));
                  };
            recover_record(full_source_record,
                           SternheimerReducedKPoint{0.0, 0.0, 0.0},
                           source_occupied_count,
                           source_sector_record);
            recover_record(full_target_record,
                           supercell_translation_sum.primitive_qpoint,
                           target_occupied_count,
                           target_sector_record);
            source_record_pointer = &source_sector_record;
            target_record_pointer = &target_sector_record;
        }
        const SternheimerLCAOOccupiedKPoint& source_record = *source_record_pointer;
        const SternheimerLCAOOccupiedKPoint& target_record = *target_record_pointer;
        const SternheimerLCAOSamplingPlan sampling_plan = sternheimer_lcao_sampling_plan(
            use_delta_sternheimer,
            write_lcao_sos,
            !target_record.unoccupied_coefficients.empty());

        std::vector<std::vector<SternheimerRPA::Complex>> partial_pair_branches;
        if (write_partial_kresolved)
        {
            partial_pair_branches.resize(static_cast<std::size_t>(nfreq));
            for (int ifrequency = 0; ifrequency != nfreq; ++ifrequency)
            {
                if (response_owned_by_rank(pair.source_index, ifrequency))
                {
                    partial_pair_branches[static_cast<std::size_t>(ifrequency)].assign(
                        static_cast<std::size_t>(num_channels) * num_channels,
                        SternheimerRPA::Complex(0.0, 0.0));
                }
            }
        }

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

        const bool reuse_target_sampling
            = can_reuse_sternheimer_target_lcao_sampling(source_record_pointer == target_record_pointer, sampling_plan);
        std::unique_ptr<SternheimerSampledLCAOKPoint> sampled_source;
        if (!reuse_target_sampling)
        {
            sampled_source = std::make_unique<SternheimerSampledLCAOKPoint>(
                sample_sternheimer_lcao_kpoint(ucell,
                                               grid_data.grid,
                                               orbitals,
                                               source_record,
                                               max_bands,
                                               sampling_plan.sample_source_unoccupied,
                                               write_lcao_sos ? 0 : delta_options.max_virtual_states,
                                               grid_data.volume_element,
                                               PARAM.inp.sternheimer_delta_norm_tol));
        }
        SternheimerSampledLCAOKPoint target
            = sample_sternheimer_lcao_kpoint(ucell,
                                             grid_data.grid,
                                             orbitals,
                                             target_record,
                                             0,
                                             sampling_plan.sample_target_unoccupied,
                                             write_lcao_sos ? 0 : delta_options.max_virtual_states,
                                             grid_data.volume_element,
                                             PARAM.inp.sternheimer_delta_norm_tol);
        const SternheimerSampledLCAOKPoint& source = reuse_target_sampling ? target : *sampled_source;
        append_chi0_progress_event("lcao_sampling_ready",
                                   0,
                                   -1,
                                   -1,
                                   -1,
                                   solved_equations,
                                   nullptr,
                                   -1.0,
                                   elapsed_seconds_since(chi0_start_time),
                                   std::string("source_reused=") + (reuse_target_sampling ? "yes" : "no")
                                       + " target_unoccupied="
                                       + std::to_string(target.unoccupied_functions.size()));
        const SternheimerFDHamiltonian hamiltonian = [&]() {
            if (response_grid.independent || use_kpoint_mpi)
            {
                SternheimerABACUSFDGridData target_grid_data = grid_data;
                target_grid_data.grid.kpoint = sternheimer_lcao_grid_kpoint(target_record);
                auto nonlocal_projector = make_sternheimer_fd_nonlocal_projector_from_unitcell(
                    ucell, target_grid_data.grid, target_grid_data.volume_element);
                return make_sternheimer_fd_hamiltonian_from_local_potential(target_grid_data,
                                                                             response_grid.independent
                                                                                 ? independent_response_full_potential
                                                                                 : kpoint_parallel_full_potential,
                                                                             1.0,
                                                                             std::move(nonlocal_projector),
                                                                             PARAM.inp.sternheimer_fd_order);
            }
            return use_frequency_mpi
                       ? make_sternheimer_fd_full_hamiltonian(
                             potential,
                             pw_basis,
                             ucell,
                             0,
                             1.0,
                             sternheimer_lcao_grid_kpoint(target_record),
                             PARAM.inp.sternheimer_fd_order)
                       : make_sternheimer_fd_hamiltonian(
                             potential,
                             pw_basis,
                             ucell,
                             0,
                             1.0,
                             sternheimer_lcao_grid_kpoint(target_record),
                             PARAM.inp.sternheimer_fd_order);
        }();

        std::vector<SternheimerFDHamiltonian::Vector> target_occupied_projector;
        target_occupied_projector.reserve(target.occupied_projector_functions.size());
        for (const SternheimerDeltaGridFunction& function: target.occupied_projector_functions)
        {
            target_occupied_projector.push_back(function.values);
        }
        SternheimerDeltaSubspace delta_subspace;
        std::unique_ptr<SternheimerPeriodicResponseProjectors> pair_projectors;
        SternheimerChannelWorkerPlan pair_channel_worker_plan = channel_worker_plan;
        target_projector_dimensions[pair_index] = static_cast<int>(target_occupied_projector.size());
        if (use_delta_sternheimer)
        {
            std::vector<SternheimerDeltaGridFunction> delta_candidates;
            int candidate_occupied_count = 0;
            if (sampling_plan.build_target_ao_candidates)
            {
                delta_candidates = build_lcao_candidate_grid_functions(
                    ucell,
                    grid_data.grid,
                    &orbitals,
                    sternheimer_lcao_grid_kpoint(target_record));
                candidate_occupied_count = static_cast<int>(target_occupied_projector.size());
            }
            else if (write_lcao_sos)
            {
                delta_candidates = target.unoccupied_functions;
            }
            else
            {
                delta_candidates = std::move(target.unoccupied_functions);
            }
            SternheimerDeltaSubspaceOptions pair_delta_options = delta_options;
            pair_delta_options.max_virtual_states = sternheimer_delta_virtual_state_limit(
                delta_options.max_virtual_states,
                static_cast<int>(delta_candidates.size()),
                candidate_occupied_count);
            pair_delta_options.retain_grid_functions = false;
            pair_delta_options.evaluate_full_grid_difference = false;
            delta_subspace = build_delta_sternheimer_subspace_by_mode(hamiltonian,
                                                                      target.occupied_projector_functions,
                                                                      std::move(delta_candidates),
                                                                      grid_data.volume_element,
                                                                      pair_delta_options,
                                                                      delta_a_block_mode);
            std::ostringstream excitation_context;
            excitation_context << std::setprecision(17)
                               << "q_index=" << response_plan.iq
                               << " source_k=" << pair.source_index + 1
                               << " target_k=" << pair.target_index + 1
                               << " rank=" << GlobalV::MY_RANK
                               << " response_grid=" << grid_data.grid.nx << 'x' << grid_data.grid.ny
                               << 'x' << grid_data.grid.nz
                               << " pbe_grid=" << pw_basis.nx << 'x' << pw_basis.ny << 'x' << pw_basis.nz
                               << " a_block=" << sternheimer_delta_a_block_mode_name(delta_a_block_mode);
            for (const auto* record : {&source_record, &target_record})
            {
                const auto& k = sternheimer_lcao_grid_kpoint(*record);
                excitation_context << " k_reduced=" << k[0] << ',' << k[1] << ',' << k[2];
            }
            require_positive_periodic_delta_excitations(
                delta_subspace.virtual_states, source_record.eigenvalues, excitation_context.str());
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

            append_chi0_progress_event("delta_subspace_ready",
                                       0,
                                       -1,
                                       -1,
                                       -1,
                                       solved_equations,
                                       nullptr,
                                       -1.0,
                                       elapsed_seconds_since(chi0_start_time),
                                       "nvirtual="
                                           + std::to_string(delta_subspace.virtual_states.size()));

            std::vector<SternheimerDeltaGridFunction>().swap(target.occupied_functions);
            std::vector<SternheimerDeltaGridFunction>().swap(target.occupied_projector_functions);
            if (!write_lcao_sos)
            {
                target.unoccupied_states = SternheimerFDZeroOrderStates();
            }
            const auto memory_before_trim = detect_sternheimer_memory_snapshot();
            const auto trim_status = trim_sternheimer_process_heap();
            const auto memory_after_trim = detect_sternheimer_memory_snapshot();
            bool pack_projectors = channel_batch_width > 1;
            std::uint64_t shared_projector_bytes
                = SternheimerPeriodicResponseProjectors::estimated_storage_bytes(
                    grid_data.grid.size(), target_occupied_projector.size(), delta_subspace.virtual_states.size(),
                    pack_projectors);
            append_chi0_progress_event(
                "delta_memory_reclaimed", 0, -1, -1, -1, solved_equations, nullptr, -1.0,
                elapsed_seconds_since(chi0_start_time),
                "resource_source=" + memory_after_trim.source
                    + " accounting_mode=" + sternheimer_memory_accounting_mode_name(memory_after_trim.mode)
                    + " memory_before_trim_bytes=" + std::to_string(memory_before_trim.current_bytes)
                    + " memory_after_trim_bytes=" + std::to_string(memory_after_trim.current_bytes)
                    + " node_memory_limit_bytes=" + std::to_string(memory_after_trim.limit_bytes)
                    + " trim_status=" + std::to_string(static_cast<int>(trim_status))
                    + " shared_projector_bytes=" + std::to_string(shared_projector_bytes));
            // Admit the persistent packed projectors before allocating them, then
            // measure them as baseline rather than duplicating them in each worker.
            if (pack_projectors)
            {
                try
                {
                    const auto admission = plan_sternheimer_channel_workers(
                        num_channels, sternheimer_channel_openmp_threads(), grid_data.grid.size(), channel_threads,
                        reserve_sternheimer_shared_memory(memory_after_trim, shared_projector_bytes), channel_batch_width);
                    pack_projectors = admission.channel_batch_width > 1;
                }
                catch (const std::runtime_error&)
                {
                    // A scalar solve borrows the same complete basis without packing it.
                    // It still has to pass the unchanged admission check below.
                    pack_projectors = false;
                }
            }
            shared_projector_bytes = SternheimerPeriodicResponseProjectors::estimated_storage_bytes(
                grid_data.grid.size(), target_occupied_projector.size(), delta_subspace.virtual_states.size(),
                pack_projectors);
            const int admitted_batch_width = pack_projectors ? channel_batch_width : 1;
            plan_sternheimer_channel_workers(
                num_channels, sternheimer_channel_openmp_threads(), grid_data.grid.size(), channel_threads,
                reserve_sternheimer_shared_memory(memory_after_trim, shared_projector_bytes), admitted_batch_width);
            pair_projectors = std::make_unique<SternheimerPeriodicResponseProjectors>(
                target_occupied_projector, delta_subspace.virtual_states, grid_data.volume_element, pack_projectors);
            trim_sternheimer_process_heap();
            SternheimerMemorySnapshot pair_channel_memory = detect_sternheimer_memory_snapshot();
            append_chi0_progress_event(
                "shared_projectors_ready", 0, -1, -1, -1, solved_equations, nullptr, -1.0,
                elapsed_seconds_since(chi0_start_time),
                "memory_current_bytes=" + std::to_string(pair_channel_memory.current_bytes)
                    + " node_memory_limit_bytes=" + std::to_string(pair_channel_memory.limit_bytes)
                    + " shared_projector_bytes=" + std::to_string(shared_projector_bytes)
                    + " owned_projector_bytes="
                    + std::to_string(pair_projectors->occupied.storage_bytes() + pair_projectors->fixed.storage_bytes()));
            pair_channel_worker_plan
                = plan_sternheimer_channel_workers(num_channels,
                                                    sternheimer_channel_openmp_threads(),
                                                    grid_data.grid.size(),
                                                    channel_threads,
                                                    pair_channel_memory,
                                                    admitted_batch_width);
            if (pack_projectors && pair_channel_worker_plan.channel_batch_width == 1)
            {
                pair_projectors.reset();
                trim_sternheimer_process_heap();
                pair_projectors = std::make_unique<SternheimerPeriodicResponseProjectors>(
                    target_occupied_projector, delta_subspace.virtual_states, grid_data.volume_element, false);
                pair_channel_memory = detect_sternheimer_memory_snapshot();
                pair_channel_worker_plan = plan_sternheimer_channel_workers(
                    num_channels, sternheimer_channel_openmp_threads(), grid_data.grid.size(), channel_threads,
                    pair_channel_memory, 1);
            }
            channel_worker_plan = pair_channel_worker_plan;
            append_chi0_progress_event(
                "equation_workers_ready",
                0,
                -1,
                -1,
                -1,
                solved_equations,
                nullptr,
                -1.0,
                elapsed_seconds_since(chi0_start_time),
                format_sternheimer_channel_worker_diagnostic(pair_channel_memory,
                                                             pair_channel_worker_plan,
                                                             grid_data.grid.size(),
                                                             channel_threads));
        }
        const int pair_channel_batch_width = pair_channel_worker_plan.channel_batch_width;

        if (!target_record.unoccupied_eigenvalues.empty())
        {
            const auto lcao_unoccupied_minmax
                = std::minmax_element(target_record.unoccupied_eigenvalues.begin(),
                                      target_record.unoccupied_eigenvalues.end());
            target_lcao_unoccupied_min[pair_index] = *lcao_unoccupied_minmax.first;
            target_lcao_unoccupied_max[pair_index] = *lcao_unoccupied_minmax.second;
        }
        else
        {
            const int target_occupied_count = static_cast<int>(target.states.eigenvalues.size());
            if (target_occupied_count >= elec_state.ekb.nc)
            {
                throw std::runtime_error("Periodic Sternheimer target k point has no LCAO unoccupied states.");
            }
            double lcao_unoccupied_min = std::numeric_limits<double>::infinity();
            double lcao_unoccupied_max = -std::numeric_limits<double>::infinity();
            for (int ib = target_occupied_count; ib != elec_state.ekb.nc; ++ib)
            {
                const double eigenvalue = elec_state.ekb(target_record.zero_order_k_index, ib);
                lcao_unoccupied_min = std::min(lcao_unoccupied_min, eigenvalue);
                lcao_unoccupied_max = std::max(lcao_unoccupied_max, eigenvalue);
            }
            target_lcao_unoccupied_min[pair_index] = lcao_unoccupied_min;
            target_lcao_unoccupied_max[pair_index] = lcao_unoccupied_max;
        }
        target_lcao_occupied_raw_norm_min[pair_index] = target.occupied_raw_norm_min;
        target_lcao_occupied_raw_norm_max[pair_index] = target.occupied_raw_norm_max;

        std::vector<SternheimerDeltaVirtualState> lcao_virtual_states;
        if (write_lcao_sos)
        {
            lcao_virtual_states.reserve(target.unoccupied_functions.size());
            for (std::size_t ia = 0; ia != target.unoccupied_functions.size(); ++ia)
            {
                SternheimerDeltaVirtualState state;
                state.orbital = target.unoccupied_functions[ia].values;
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
            const int owner_rank = response_owner_rank(pair.source_index, ifrequency);
            if (!response_owned_by_rank(pair.source_index, ifrequency))
            {
                continue;
            }
            const double omega_ry = 2.0 * frequency_grid.omega_ha[static_cast<std::size_t>(ifrequency)];
            std::vector<SternheimerRPA::Complex>& chi0_branch
                = write_partial_kresolved
                      ? partial_pair_branches[static_cast<std::size_t>(ifrequency)]
                      : chi0_branches[static_cast<std::size_t>(ifrequency)];
            std::vector<SternheimerRPA::Complex>* delta_sos_branch
                = write_delta_components ? &delta_sos_branches[static_cast<std::size_t>(ifrequency)] : nullptr;
            std::vector<SternheimerRPA::Complex>* delta_pulay_branch
                = write_delta_components ? &delta_pulay_branches[static_cast<std::size_t>(ifrequency)] : nullptr;
            std::vector<SternheimerRPA::Complex>* delta_out_branch
                = write_delta_components ? &delta_out_branches[static_cast<std::size_t>(ifrequency)] : nullptr;
            std::vector<SternheimerRPA::Complex>* lcao_sos_branch
                = write_lcao_sos ? &lcao_sos_branches[static_cast<std::size_t>(ifrequency)] : nullptr;

            const int source_band_count = sternheimer_periodic_band_count(
                static_cast<int>(source.states.wavefunctions.size()), max_bands);
            for (int ib = 0; ib != source_band_count; ++ib)
            {
                const double occupation = sternheimer_lcao_weighted_occupation(source_record, ib);
                const double matrix_occupation = response_matrix_scale * occupation;
                struct PeriodicChannelEquationResult
                {
                    int channel_index = -1;
                    SternheimerRPA::SolverResult solver;
                    double equation_residual_norm = 0.0;
                    double full_grid_equation_residual_norm = 0.0;
                    bool has_wavefunction_diagnostic = false;
                    SternheimerWavefunctionDiagnostic::Record wavefunction_diagnostic;
                };
                const std::vector<SternheimerChannelBatch> channel_batches
                    = make_sternheimer_owned_channel_batches(num_channels, channel_batch_width,
                                                             pair_channel_batch_width, response_replicas,
                                                             local_channel_replica, ib);
                const int solved_equations_before_band = solved_equations;
                std::atomic<int> completed_batch_equations{0};
                std::vector<std::vector<PeriodicChannelEquationResult>> grouped_channel_results
                    = run_sternheimer_channel_tasks<std::vector<PeriodicChannelEquationResult>>(
                        static_cast<int>(channel_batches.size()),
                        [&](const int batch_task) {
                            const SternheimerChannelBatch batch
                                = channel_batches[static_cast<std::size_t>(batch_task)];
                            SternheimerFDHamiltonian::Matrix rhs_batch(static_cast<std::size_t>(batch.size));
                            std::vector<SternheimerFDHamiltonian::Vector> batch_perturbations_ry(
                                static_cast<std::size_t>(batch.size));
                            std::vector<std::vector<SternheimerFDHamiltonian::Complex>>
                                perturbation_matrix_elements(static_cast<std::size_t>(batch.size));
                            for (int offset = 0; offset != batch.size; ++offset)
                            {
                                const int ichannel = batch.begin + offset;
                                const std::size_t channel_index = static_cast<std::size_t>(ichannel);
                                auto& perturbation_ry = batch_perturbations_ry[static_cast<std::size_t>(offset)];
                                perturbation_ry = channels[channel_index].potential_r;
                                for (auto& value: perturbation_ry)
                                {
                                    value *= kHartreeToRydberg;
                                }
                                SternheimerRPA::build_rhs_from_hartree_perturbation(
                                    perturbation_ry,
                                    source.states.wavefunctions[ib],
                                    rhs_batch[static_cast<std::size_t>(offset)]);
                                if (use_delta_sternheimer)
                                {
                                    perturbation_matrix_elements[static_cast<std::size_t>(offset)]
                                        = delta_sternheimer_perturbation_matrix_elements(
                                            delta_subspace.virtual_states,
                                            perturbation_ry,
                                            source.states.wavefunctions[ib],
                                            grid_data.volume_element);
                                }
                            }
                            std::vector<SternheimerPeriodicLinearResponse> responses;
                            if (pair_channel_batch_width == 1 || !use_delta_sternheimer)
                            {
                                responses.push_back(solve_sternheimer_periodic_linear_response(
                                    use_delta_sternheimer,
                                    hamiltonian,
                                    target_occupied_projector,
                                    source.states.eigenvalues[ib],
                                    rhs_batch.front(),
                                    delta_subspace.virtual_states,
                                    perturbation_matrix_elements.front(),
                                    omega_ry,
                                    grid_data.volume_element,
                                    solver_options,
                                    pair_projectors.get()));
                            }
                            else
                            {
                                responses = solve_sternheimer_periodic_linear_response_batch(
                                    true,
                                    hamiltonian,
                                    target_occupied_projector,
                                    source.states.eigenvalues[ib],
                                    rhs_batch,
                                    delta_subspace.virtual_states,
                                    perturbation_matrix_elements,
                                    omega_ry,
                                    grid_data.volume_element,
                                    solver_options,
                                    pair_projectors.get());
                            }
                            std::vector<const SternheimerFDHamiltonian::Vector*> response_views;
                            response_views.reserve(responses.size());
                            for (const auto& response: responses)
                            {
                                response_views.push_back(&response.wavefunction);
                            }
                            SternheimerRPA::accumulate_chi0_branch_columns(channels,
                                                                           source.states.wavefunctions[ib],
                                                                           response_views,
                                                                           grid_data.volume_element,
                                                                           matrix_occupation,
                                                                           batch.begin,
                                                                           chi0_branch);
                            std::vector<PeriodicChannelEquationResult> batch_results;
                            batch_results.reserve(static_cast<std::size_t>(batch.size));
                            for (int offset = 0; offset != batch.size; ++offset)
                            {
                                const int ichannel = batch.begin + offset;
                                const std::size_t channel_index = static_cast<std::size_t>(ichannel);
                                const SternheimerPeriodicLinearResponse& response
                                    = responses[static_cast<std::size_t>(offset)];
                                if (write_delta_components)
                                {
                                    SternheimerRPA::accumulate_chi0_branch_column(
                                        channels,
                                        source.states.wavefunctions[ib],
                                        response.delta_components.in_sos_wavefunction,
                                        grid_data.volume_element,
                                        matrix_occupation,
                                        ichannel,
                                        *delta_sos_branch);
                                    SternheimerRPA::accumulate_chi0_branch_column(
                                        channels,
                                        source.states.wavefunctions[ib],
                                        response.delta_components.in_pulay_wavefunction,
                                        grid_data.volume_element,
                                        matrix_occupation,
                                        ichannel,
                                        *delta_pulay_branch);
                                    SternheimerRPA::accumulate_chi0_branch_column(
                                        channels,
                                        source.states.wavefunctions[ib],
                                        response.delta_components.out_wavefunction,
                                        grid_data.volume_element,
                                        matrix_occupation,
                                        ichannel,
                                        *delta_out_branch);
                                }
                                SternheimerFDHamiltonian::Vector lcao_response;
                                if (write_lcao_sos)
                                {
                                    const auto lcao_matrix_elements = delta_sternheimer_perturbation_matrix_elements(
                                        lcao_virtual_states,
                                        batch_perturbations_ry[static_cast<std::size_t>(offset)],
                                        source.states.wavefunctions[ib],
                                        grid_data.volume_element);
                                    lcao_response = build_delta_sternheimer_sos_wavefunction(
                                        lcao_virtual_states,
                                        lcao_matrix_elements,
                                        source.states.eigenvalues[ib],
                                        omega_ry);
                                    SternheimerRPA::accumulate_chi0_branch_column(channels,
                                                                                  source.states.wavefunctions[ib],
                                                                                  lcao_response,
                                                                                  grid_data.volume_element,
                                                                                  matrix_occupation,
                                                                                  ichannel,
                                                                                  *lcao_sos_branch);
                                }
                                PeriodicChannelEquationResult result;
                                result.channel_index = ichannel;
                                result.solver = response.solver;
                                result.equation_residual_norm = response.residual_norm;
                                result.full_grid_equation_residual_norm
                                    = response.full_grid_equation_residual_norm;
                                if (write_wavefunction_diagnostic
                                    && wavefunction_diagnostic_config.selector.matches(response_plan.iq,
                                                                                       pair.source_index,
                                                                                       ib,
                                                                                       ifrequency + 1,
                                                                                       ichannel))
                                {
                                    result.has_wavefunction_diagnostic = true;
                                    auto& diagnostic = result.wavefunction_diagnostic;
                                    diagnostic.metadata.nx = grid_data.grid.nx;
                                    diagnostic.metadata.ny = grid_data.grid.ny;
                                    diagnostic.metadata.nz = grid_data.grid.nz;
                                    diagnostic.metadata.iq = response_plan.iq;
                                    diagnostic.metadata.ik_full = pair.source_index;
                                    diagnostic.metadata.ib = ib;
                                    diagnostic.metadata.ifrequency = ifrequency + 1;
                                    diagnostic.metadata.channel = ichannel;
                                    const auto lattice = sternheimer_fd_grid_lattice_vectors(grid_data.grid);
                                    for (int i = 0; i != 3; ++i)
                                    {
                                        for (int j = 0; j != 3; ++j)
                                        {
                                            diagnostic.metadata.lattice[static_cast<std::size_t>(3 * i + j)]
                                                = lattice[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
                                        }
                                    }
                                    diagnostic.metadata.qpoint
                                        = use_supercell_translation_sum
                                              ? supercell_translation_sum.primitive_qpoint
                                              : response_plan.qpoint;
                                    diagnostic.metadata.source_kpoint = source_record.kpoint;
                                    diagnostic.metadata.target_kpoint = target_record.kpoint;
                                    diagnostic.metadata.omega_ha
                                        = frequency_grid.omega_ha[static_cast<std::size_t>(ifrequency)];
                                    diagnostic.metadata.omega_ry = omega_ry;
                                    diagnostic.metadata.volume_element = grid_data.volume_element;
                                    diagnostic.metadata.reference_eigenvalue_ry = source.states.eigenvalues[ib];
                                    diagnostic.metadata.weighted_occupation = occupation;
                                    diagnostic.metadata.rhs_norm
                                        = sternheimer_fd_grid_norm(response.projected_rhs,
                                                                   grid_data.volume_element);
                                    diagnostic.metadata.solver_relative_residual = response.solver.relative_residual;
                                    diagnostic.metadata.equation_relative_residual
                                        = diagnostic.metadata.rhs_norm > 0.0
                                              ? response.full_grid_equation_residual_norm
                                                    / diagnostic.metadata.rhs_norm
                                              : response.full_grid_equation_residual_norm;
                                    diagnostic.metadata.diagonal_branch_element
                                        = matrix_occupation
                                          * SternheimerRPA::accumulate_polarizability_grid_element(
                                              channels[channel_index].potential_r,
                                              source.states.wavefunctions[ib],
                                              response.wavefunction,
                                              grid_data.volume_element);
                                    diagnostic.vectors.push_back(
                                        {"psi0", source.states.wavefunctions[ib]});
                                    diagnostic.vectors.push_back(
                                        {"hartree_potential_ha", channels[channel_index].potential_r});
                                    diagnostic.vectors.push_back(
                                        {"perturbation_ry", batch_perturbations_ry[static_cast<std::size_t>(offset)]});
                                    diagnostic.vectors.push_back({"rhs_ry", response.projected_rhs});
                                    diagnostic.vectors.push_back({"delta_psi", response.wavefunction});
                                    if (response.has_delta_components)
                                    {
                                        diagnostic.vectors.push_back(
                                            {"delta_in_sos", response.delta_components.in_sos_wavefunction});
                                        diagnostic.vectors.push_back(
                                            {"delta_in_pulay", response.delta_components.in_pulay_wavefunction});
                                        diagnostic.vectors.push_back(
                                            {"delta_out_grid", response.delta_components.out_wavefunction});
                                        const SternheimerDeltaPulayOperatorComponents pulay_terms
                                            = decompose_delta_sternheimer_pulay_operator_terms(
                                                hamiltonian,
                                                diagnostic_fixed_local_potential,
                                                target_occupied_projector,
                                                delta_subspace.virtual_states,
                                                response.delta_components.out_wavefunction,
                                                source.states.eigenvalues[ib],
                                                omega_ry,
                                                grid_data.volume_element);
                                        diagnostic.vectors.push_back(
                                            {"delta_in_pulay_kinetic", pulay_terms.kinetic});
                                        diagnostic.vectors.push_back(
                                            {"delta_in_pulay_fixed_local", pulay_terms.fixed_local});
                                        diagnostic.vectors.push_back(
                                            {"delta_in_pulay_hxc_local", pulay_terms.hxc_local});
                                        diagnostic.vectors.push_back(
                                            {"delta_in_pulay_nonlocal", pulay_terms.nonlocal});
                                        diagnostic.vectors.push_back(
                                            {"delta_in_pulay_eigenvalue", pulay_terms.eigenvalue});
                                        diagnostic.vectors.push_back(
                                            {"delta_in_pulay_operator_total", pulay_terms.total});
                                    }
                                    if (!lcao_response.empty())
                                    {
                                        diagnostic.vectors.push_back({"lcao_sos", std::move(lcao_response)});
                                    }
                                }
                                    batch_results.push_back(std::move(result));
                            }
                            SternheimerRPA::SolverResult batch_solver;
                            batch_solver.converged = true;
                            double batch_equation_residual_norm = 0.0;
                            for (const PeriodicChannelEquationResult& result: batch_results)
                            {
                                batch_solver.converged = batch_solver.converged && result.solver.converged;
                                batch_solver.iterations = std::max(batch_solver.iterations,
                                                                   result.solver.iterations);
                                batch_solver.relative_residual = std::max(batch_solver.relative_residual,
                                                                          result.solver.relative_residual);
                                batch_solver.absolute_residual = std::max(batch_solver.absolute_residual,
                                                                          result.solver.absolute_residual);
                                batch_equation_residual_norm
                                    = std::max(batch_equation_residual_norm,
                                               result.equation_residual_norm);
                            }
                            const int batch_completed
                                = completed_batch_equations.fetch_add(batch.size) + batch.size;
                            append_chi0_progress_event(
                                "equation_batch",
                                ifrequency + 1,
                                owner_rank,
                                ib,
                                batch.begin,
                                solved_equations_before_band + batch_completed,
                                &batch_solver,
                                batch_equation_residual_norm,
                                elapsed_seconds_since(chi0_start_time),
                                "batch_size=" + std::to_string(batch.size)
                                    + ",source_k=" + std::to_string(pair.source_index + 1)
                                    + ",target_k=" + std::to_string(pair.target_index + 1));
                            return batch_results;
                        },
                        pair_channel_worker_plan.effective_workers);
                std::vector<PeriodicChannelEquationResult> channel_results;
                channel_results.reserve(static_cast<std::size_t>(num_channels / response_replicas + 1));
                for (std::vector<PeriodicChannelEquationResult>& group: grouped_channel_results)
                {
                    for (PeriodicChannelEquationResult& result: group)
                    {
                        channel_results.push_back(std::move(result));
                    }
                }

                for (PeriodicChannelEquationResult& result: channel_results)
                {
                    all_converged = all_converged && result.solver.converged;
                    ++solved_equations;
                    max_solver_relative_residual
                        = std::max(max_solver_relative_residual, result.solver.relative_residual);
                    max_equation_residual_norm
                        = std::max(max_equation_residual_norm, result.equation_residual_norm);
                    max_full_grid_equation_residual_norm
                        = std::max(max_full_grid_equation_residual_norm,
                                   result.full_grid_equation_residual_norm);
                    if (result.has_wavefunction_diagnostic)
                    {
                        SternheimerWavefunctionDiagnostic::write(
                            wavefunction_diagnostic_config.output_filename,
                            result.wavefunction_diagnostic);
                        ++local_wavefunction_diagnostic_count;
                        GlobalV::ofs_running << " Sternheimer wavefunction diagnostic: "
                                             << wavefunction_diagnostic_config.output_filename << std::endl;
                    }
                    append_chi0_progress_event("equation",
                                               ifrequency + 1,
                                               owner_rank,
                                               ib,
                                               result.channel_index,
                                               solved_equations,
                                               &result.solver,
                                               result.equation_residual_norm,
                                               elapsed_seconds_since(chi0_start_time),
                                               "source_k=" + std::to_string(pair.source_index + 1)
                                                   + ",target_k=" + std::to_string(pair.target_index + 1));
                }
            }
        }
        if (write_partial_kresolved)
        {
#ifdef __MPI
            if (response_replicas > 1)
            {
                for (auto& branch: partial_pair_branches)
                {
                    if (!branch.empty())
                    {
                        MPI_Allreduce(MPI_IN_PLACE,
                                      branch.data(),
                                      static_cast<int>(branch.size()),
                                      MPI_DOUBLE_COMPLEX,
                                      MPI_SUM,
                                      response_group_comm);
                    }
                }
            }
#endif
            for (int ifrequency = 0; ifrequency != nfreq; ++ifrequency)
            {
                const auto& branch = partial_pair_branches[static_cast<std::size_t>(ifrequency)];
                if (branch.empty() || local_channel_replica != 0)
                {
                    continue;
                }
                SternheimerPartialResponseRecord partial
                    = make_sternheimer_partial_response_record(response_plan.iq,
                                                               pair.source_index,
                                                               ifrequency + 1,
                                                               branch,
                                                               num_channels);
                const SternheimerRPA::Chi0V1Metadata metadata
                    = make_chi0_v1_metadata(ucell,
                                            channels,
                                            response_plan.iq,
                                            ifrequency + 1,
                                            frequency_grid.omega_ha[static_cast<std::size_t>(ifrequency)],
                                            frequency_grid.weights_ha[static_cast<std::size_t>(ifrequency)],
                                            output_atom_count);
                SternheimerRPA::write_chi0_v1_file(
                    partial.filename, metadata, auxiliary_channels, partial.matrix);
                GlobalV::ofs_running << " Sternheimer periodic partial response: "
                                     << partial.filename << std::endl;
                local_partial_records.push_back(std::move(partial));

                if (write_kresolved_diagnostic)
                {
                    auto& aggregate = chi0_branches[static_cast<std::size_t>(ifrequency)];
                    for (std::size_t entry = 0; entry != branch.size(); ++entry)
                    {
                        aggregate[entry] += branch[entry];
                    }
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

#ifdef __MPI
    if (response_replicas > 1 && !write_partial_kresolved)
    {
        const auto reduce_response_replica_branches
            = [&](std::vector<std::vector<SternheimerRPA::Complex>>& branches) {
                  for (auto& branch: branches)
                  {
                      if (branch.empty())
                      {
                          continue;
                      }
                      MPI_Allreduce(MPI_IN_PLACE,
                                    branch.data(),
                                    static_cast<int>(branch.size()),
                                    MPI_DOUBLE_COMPLEX,
                                    MPI_SUM,
                                    response_group_comm);
                      if (local_channel_replica != 0)
                      {
                          std::fill(branch.begin(), branch.end(), SternheimerRPA::Complex(0.0, 0.0));
                      }
                  }
              };
        reduce_response_replica_branches(chi0_branches);
        if (write_delta_components)
        {
            reduce_response_replica_branches(delta_sos_branches);
            reduce_response_replica_branches(delta_pulay_branches);
            reduce_response_replica_branches(delta_out_branches);
        }
        if (write_lcao_sos)
        {
            reduce_response_replica_branches(lcao_sos_branches);
        }
    }
    if (response_group_comm != MPI_COMM_NULL)
    {
        MPI_Comm_free(&response_group_comm);
    }
#endif

    if (use_nested_response_mpi
        && (local_frequency_slot != 0 || local_channel_replica != 0))
    {
        std::fill(target_projector_dimensions.begin(), target_projector_dimensions.end(), 0);
        std::fill(target_delta_dimensions.begin(), target_delta_dimensions.end(), 0);
        std::fill(target_delta_eigenvalue_min.begin(), target_delta_eigenvalue_min.end(), 0.0);
        std::fill(target_delta_eigenvalue_max.begin(), target_delta_eigenvalue_max.end(), 0.0);
        std::fill(target_delta_grid_hamiltonian_relative_difference.begin(),
                  target_delta_grid_hamiltonian_relative_difference.end(),
                  0.0);
        std::fill(target_delta_grid_hamiltonian_max_abs_difference.begin(),
                  target_delta_grid_hamiltonian_max_abs_difference.end(),
                  0.0);
        std::fill(target_lcao_unoccupied_min.begin(), target_lcao_unoccupied_min.end(), 0.0);
        std::fill(target_lcao_unoccupied_max.begin(), target_lcao_unoccupied_max.end(), 0.0);
        std::fill(target_lcao_occupied_raw_norm_min.begin(), target_lcao_occupied_raw_norm_min.end(), 0.0);
        std::fill(target_lcao_occupied_raw_norm_max.begin(), target_lcao_occupied_raw_norm_max.end(), 0.0);
        std::fill(target_lcao_unoccupied_raw_norm_min.begin(), target_lcao_unoccupied_raw_norm_min.end(), 0.0);
        std::fill(target_lcao_unoccupied_raw_norm_max.begin(), target_lcao_unoccupied_raw_norm_max.end(), 0.0);
        std::fill(target_lcao_occ_unocc_overlap_max.begin(), target_lcao_occ_unocc_overlap_max.end(), 0.0);
    }

    if (use_kpoint_mpi)
    {
        if (!use_symmetry_partial_response)
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

    if (!use_symmetry_partial_response && write_periodic_v1)
    {
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
                                        frequency_grid.weights_ha[static_cast<std::size_t>(ifrequency)],
                                        output_atom_count);
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
    }
    int global_wavefunction_diagnostic_count = local_wavefunction_diagnostic_count;
#ifdef __MPI
    MPI_Allreduce(&local_wavefunction_diagnostic_count,
                  &global_wavefunction_diagnostic_count,
                  1,
                  MPI_INT,
                  MPI_SUM,
                  MPI_COMM_WORLD);
#endif
    if (write_wavefunction_diagnostic && global_wavefunction_diagnostic_count != 1)
    {
        throw std::runtime_error("Sternheimer wavefunction diagnostic selector matched "
                                 + std::to_string(global_wavefunction_diagnostic_count)
                                 + " equations; expected exactly one.");
    }

    std::vector<SternheimerPartialResponseRecord> gathered_partial_records;
    if (write_partial_kresolved)
    {
        std::vector<int> local_keys;
        local_keys.reserve(3 * local_partial_records.size());
        for (const auto& record: local_partial_records)
        {
            local_keys.push_back(record.iq);
            local_keys.push_back(record.ik_full);
            local_keys.push_back(record.ifrequency);
        }
#ifdef __MPI
        const int local_record_count = static_cast<int>(local_partial_records.size());
        std::vector<int> record_counts;
        if (GlobalV::MY_RANK == 0)
        {
            record_counts.resize(static_cast<std::size_t>(GlobalV::NPROC), 0);
        }
        MPI_Gather(&local_record_count,
                   1,
                   MPI_INT,
                   record_counts.data(),
                   1,
                   MPI_INT,
                   0,
                   MPI_COMM_WORLD);

        std::vector<int> key_counts;
        std::vector<int> key_displacements;
        std::vector<int> gathered_keys;
        if (GlobalV::MY_RANK == 0)
        {
            key_counts.resize(record_counts.size(), 0);
            key_displacements.resize(record_counts.size(), 0);
            int total_key_count = 0;
            for (std::size_t rank = 0; rank != record_counts.size(); ++rank)
            {
                key_counts[rank] = 3 * record_counts[rank];
                key_displacements[rank] = total_key_count;
                total_key_count += key_counts[rank];
            }
            gathered_keys.resize(static_cast<std::size_t>(total_key_count));
        }
        MPI_Gatherv(local_keys.data(),
                    static_cast<int>(local_keys.size()),
                    MPI_INT,
                    gathered_keys.data(),
                    key_counts.data(),
                    key_displacements.data(),
                    MPI_INT,
                    0,
                    MPI_COMM_WORLD);
        if (GlobalV::MY_RANK == 0)
        {
            gathered_partial_records.reserve(gathered_keys.size() / 3);
            for (std::size_t index = 0; index != gathered_keys.size(); index += 3)
            {
                SternheimerPartialResponseRecord record;
                record.iq = gathered_keys[index];
                record.ik_full = gathered_keys[index + 1];
                record.ifrequency = gathered_keys[index + 2];
                record.filename = sternheimer_partial_response_filename(
                    record.iq, record.ik_full, record.ifrequency);
                gathered_partial_records.push_back(std::move(record));
            }
        }
#else
        gathered_partial_records = local_partial_records;
#endif
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
    if (write_partial_kresolved)
    {
        const std::size_t expected_records
            = (use_symmetry_partial_response ? fixed_q_orbits.size() : response_plan.kq_pairs.size())
              * static_cast<std::size_t>(nfreq);
        if (gathered_partial_records.size() != expected_records)
        {
            throw std::runtime_error(
                "Sternheimer k-resolved output did not produce one partial response per expected k record and frequency.");
        }
        const std::string manifest_file = use_symmetry_partial_response
                                              ? "v1_sternheimer_partial_manifest_iq_"
                                                    + std::to_string(response_plan.iq) + ".dat"
                                              : "v1_sternheimer_kresolved_manifest_iq_"
                                                    + std::to_string(response_plan.iq) + ".dat";
        std::ofstream manifest(manifest_file);
        if (!manifest)
        {
            throw std::runtime_error("Cannot open Sternheimer partial-response manifest: " + manifest_file);
        }
        manifest << format_sternheimer_partial_manifest(gathered_partial_records);

        constexpr const char* full_kpoint_manifest = "v1_sternheimer_full_kpoints.dat";
        std::ofstream full_kpoint_out(full_kpoint_manifest);
        if (!full_kpoint_out)
        {
            throw std::runtime_error("Cannot open Sternheimer full-k-point manifest: "
                                     + std::string(full_kpoint_manifest));
        }
        full_kpoint_out << format_sternheimer_full_kpoint_manifest(response_kpoints);

        if (use_symmetry_partial_response)
        {
            const std::string route_fragment
                = "v1_sternheimer_symmetry_routes_iq_"
                  + std::to_string(response_plan.iq) + ".dat";
            std::ofstream route_out(route_fragment);
            if (!route_out)
            {
                throw std::runtime_error(
                    "Cannot open Sternheimer fixed-q route fragment: " + route_fragment);
            }
            route_out << format_sternheimer_fixed_q_routes(fixed_q_routes);

            std::vector<SternheimerQStarRoute> selected_qstar_routes;
            std::copy_if(qstar_routes.begin(),
                         qstar_routes.end(),
                         std::back_inserter(selected_qstar_routes),
                         [&](const auto& route) {
                             return route.representative_iq == response_plan.iq;
                         });
            if (selected_qstar_routes.empty())
            {
                throw std::runtime_error(
                    "The selected Sternheimer discrete q-star has no route members.");
            }
            const std::string qstar_route_fragment
                = "v1_sternheimer_qstar_routes_iq_"
                  + std::to_string(response_plan.iq) + ".dat";
            std::ofstream qstar_route_out(qstar_route_fragment);
            if (!qstar_route_out)
            {
                throw std::runtime_error(
                    "Cannot open Sternheimer q-star route fragment: " + qstar_route_fragment);
            }
            qstar_route_out << format_sternheimer_qstar_routes(selected_qstar_routes);

            const int qstar_size = static_cast<int>(selected_qstar_routes.size());
            const double qweight
                = static_cast<double>(qstar_size) / static_cast<double>(response_kpoints.size());
            const std::string qpoint_fragment
                = "v1_sternheimer_qpoint_iq_" + std::to_string(response_plan.iq) + ".dat";
            std::ofstream qpoint_out(qpoint_fragment);
            if (!qpoint_out)
            {
                throw std::runtime_error("Cannot open Sternheimer q-point fragment: " + qpoint_fragment);
            }
            qpoint_out << std::setprecision(17) << response_plan.iq << ' ' << response_plan.qpoint[0] << ' '
                       << response_plan.qpoint[1] << ' ' << response_plan.qpoint[2] << ' ' << qweight << '\n';
        }
    }
    std::vector<std::pair<std::string, SternheimerRPA::Chi0V1Metadata>> index_entries;
    if (!use_symmetry_partial_response && write_periodic_v1)
    {
        index_entries.reserve(static_cast<std::size_t>(nfreq));
        for (int ifrequency = 0; ifrequency != nfreq; ++ifrequency)
        {
            const SternheimerRPA::Chi0V1Metadata metadata
                = make_chi0_v1_metadata(ucell,
                                        channels,
                                        response_plan.iq,
                                        ifrequency + 1,
                                        frequency_grid.omega_ha[static_cast<std::size_t>(ifrequency)],
                                        frequency_grid.weights_ha[static_cast<std::size_t>(ifrequency)],
                                        output_atom_count);
            index_entries.push_back(
                {chi0_v1_filename(metadata.iq,
                                  metadata.ifrequency,
                                  frequency_owners[static_cast<std::size_t>(ifrequency)]),
                 metadata});
        }
        write_chi0_index_file("v1_sternheimer_chi0_index.dat", index_entries);
    }

    out << "status success\n";
    out << "format " << (!write_periodic_v1
                              ? "diagnostic_only"
                              : (write_truncated_v1_diagnostic
                                     ? (use_symmetry_partial_response
                                            ? "v1_partial_diagnostic_truncated_bands"
                                            : "v1_diagnostic_truncated_bands")
                                     : (use_symmetry_partial_response
                                            ? "v1_partial"
                                            : (write_kresolved_diagnostic ? "v1_kresolved" : "v1"))))
        << '\n';
    if (write_truncated_v1_diagnostic)
    {
        out << "diagnostic_only_reason truncated_bands\n";
        out << "physical_result no\n";
    }
    if (!write_periodic_v1)
    {
        out << "data_files 0\n";
        out << "response_output diagnostic_only_no_v1\n";
        out << "diagnostic_only_reason "
            << (use_supercell_translation_sum ? "supercell_translation_sum" : "truncated_bands") << '\n';
    }
    else if (use_symmetry_partial_response)
    {
        out << "data_files " << gathered_partial_records.size() << '\n';
        out << "partial_manifest v1_sternheimer_partial_manifest_iq_" << response_plan.iq << ".dat\n";
        out << "symmetry_routes v1_sternheimer_symmetry_routes_iq_" << response_plan.iq << ".dat\n";
        out << "qstar_routes v1_sternheimer_qstar_routes_iq_" << response_plan.iq << ".dat\n";
        out << "qpoint_fragment v1_sternheimer_qpoint_iq_" << response_plan.iq << ".dat\n";
        out << "full_kpoints_manifest v1_sternheimer_full_kpoints.dat\n";
        out << "discrete_spatial_group_order " << fixed_q_discrete_spatial_order << '\n';
        out << "fixed_q_little_group_order " << fixed_q_little_group_order << '\n';
    }
    else if (write_kresolved_diagnostic)
    {
        out << "data_files " << index_entries.size() << '\n';
        out << "index_file v1_sternheimer_chi0_index.dat\n";
        out << "kresolved_data_files " << gathered_partial_records.size() << '\n';
        out << "kresolved_manifest v1_sternheimer_kresolved_manifest_iq_" << response_plan.iq << ".dat\n";
        out << "full_kpoints_manifest v1_sternheimer_full_kpoints.dat\n";
    }
    else
    {
        out << "data_files " << index_entries.size() << '\n';
        out << "index_file v1_sternheimer_chi0_index.dat\n";
    }
    out << "sternheimer_q_index " << response_plan.iq << '\n';
    out << "sternheimer_qpoint " << response_plan.qpoint[0] << ' ' << response_plan.qpoint[1] << ' '
        << response_plan.qpoint[2] << '\n';
    out << "sternheimer_kweight_sum " << response_plan.kweight_sum << '\n';
    out << "sternheimer_kq_pairs " << response_plan.kq_pairs.size() << '\n';
    out << "sternheimer_canonical_q_indices";
    for (const int iq: canonical_q_indices)
    {
        out << ' ' << iq;
    }
    out << '\n';
    out << "sternheimer_fixed_q_little_group_order " << fixed_q_little_group_order << '\n';
    out << "sternheimer_fixed_q_representatives "
        << (use_symmetry_partial_response ? fixed_q_orbits.size() : response_plan.kq_pairs.size()) << '\n';
    out << "source_global_k target_global_k fixed_q_representative gx gy gz target_projector_dimension target_delta_dimension "
           "delta_eigenvalue_min_Ry delta_eigenvalue_max_Ry lcao_unoccupied_min_Ry "
           "lcao_unoccupied_max_Ry lcao_occupied_raw_norm_min lcao_occupied_raw_norm_max "
           "lcao_unoccupied_raw_norm_min lcao_unoccupied_raw_norm_max lcao_occ_unocc_overlap_max "
           "delta_grid_hamiltonian_relative_difference delta_grid_hamiltonian_max_abs_difference_Ry\n";
    for (std::size_t pair_index = 0; pair_index != response_plan.kq_pairs.size(); ++pair_index)
    {
        const SternheimerKQPair& pair = response_plan.kq_pairs[pair_index];
        out << pair.source_index + 1 << ' ' << pair.target_index + 1 << ' '
            << (fixed_q_representative[static_cast<std::size_t>(pair.source_index)] ? "yes" : "no") << ' '
            << pair.reciprocal_shift[0] << ' '
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
    out << "pbe_ecutwfc_Ry " << PARAM.inp.ecutwfc << '\n';
    out << "sternheimer_response_ecutwfc_requested_Ry " << PARAM.inp.sternheimer_response_ecutwfc << '\n';
    out << "sternheimer_response_ecutwfc_Ry " << response_grid_ecutwfc << '\n';
    out << "sternheimer_response_grid_requested " << response_grid.requested_dimensions[0] << ' '
        << response_grid.requested_dimensions[1] << ' ' << response_grid.requested_dimensions[2] << '\n';
    out << "sternheimer_response_grid_source " << response_grid_source << '\n';
    out << "pbe_grid " << pw_basis.nx << ' ' << pw_basis.ny << ' ' << pw_basis.nz << " size " << pw_basis.nxyz
        << '\n';
    out << "response_grid " << grid_data.grid.nx << ' ' << grid_data.grid.ny << ' ' << grid_data.grid.nz << " size "
        << grid_data.grid.size() << '\n';
    out << "nfreq " << nfreq << '\n';
    out << "frequency_grid_source " << frequency_grid_source << '\n';
    if (use_frequency_grid_file)
    {
        out << "frequency_grid_file " << frequency_grid_file << '\n';
    }
    out << "transition_window_Ha " << transition_window.emin_ha << ' ' << transition_window.emax_ha << '\n';
    out << "pca_threshold " << pca_threshold << '\n';
    out << "sternheimer_abfs_source " << (reuse_rpa_abfs ? "rpa_handoff" : "rebuilt") << '\n';
    out << "ccp_rmesh_times_input " << ccp_rmesh_times << '\n';
    out << "perturbation_coulomb_kernel full_periodic_poisson\n";
    out << "periodic_kmesh " << response_kmesh[0] << ' ' << response_kmesh[1] << ' '
        << response_kmesh[2] << '\n';
    out << "periodic_gamma_massidda_chi " << massidda_chi << '\n';
    out << "periodic_gamma_coulomb_projection "
        << (use_supercell_translation_sum
                ? "skipped_supercell_translation_sum"
                : (evaluate_gamma_projection ? "diagnostic_only_physical_poisson"
                                             : (gamma_qpoint ? "skipped_large_auxiliary_dimension" : "not_applicable")))
        << '\n';
    out << "periodic_gamma_projection_relative_error " << gamma_projection_relative_error << '\n';
    out << "periodic_gamma_limit "
        << (gamma_qpoint ? "constant_mode_only_no_headwing" : "not_applicable") << '\n';
    out << "perturbation_ccp_rmesh_times_used no\n";
    out << "supercell_translation_perturbation " << (use_supercell_translation_sum ? "yes" : "no") << '\n';
    out << "supercell_translation_full_response " << (full_supercell_response ? "yes" : "no") << '\n';
    out << "supercell_abfs_root_build_broadcast "
        << (full_supercell_response && use_kpoint_mpi ? "yes" : "no") << '\n';
    out << "supercell_response_matrix_scale " << response_matrix_scale << '\n';
    if (use_supercell_translation_sum)
    {
        out << "supercell_translation_repeats " << supercell_translation_sum.repeats[0] << ' '
            << supercell_translation_sum.repeats[1] << ' ' << supercell_translation_sum.repeats[2] << '\n';
        out << "supercell_translation_primitive_q " << supercell_translation_sum.primitive_qpoint[0] << ' '
            << supercell_translation_sum.primitive_qpoint[1] << ' '
            << supercell_translation_sum.primitive_qpoint[2] << '\n';
        out << "supercell_translation_atoms_per_primitive "
            << supercell_translation_sum.atoms_per_primitive << '\n';
        out << "supercell_translation_basis_atom " << supercell_translation_sum.basis_atom << '\n';
        out << "supercell_translation_channel_within_atom "
            << supercell_translation_sum.channel_within_atom << '\n';
        out << "supercell_sector_dimension " << supercell_sector_dimension << '\n';
        out << "supercell_sector_kweight " << supercell_sector_kweight << '\n';
        out << "supercell_sector_max_orthonormality_error "
            << supercell_sector_max_orthonormality_error << '\n';
        out << "supercell_sector_max_full_space_residual "
            << supercell_sector_max_full_space_residual << '\n';
    }
    out << "sternheimer_channel_threads " << channel_worker_plan.effective_workers << '\n';
    out << "sternheimer_channel_batch_width " << channel_worker_plan.channel_batch_width << '\n';
    out << "sternheimer_mode " << (use_delta_sternheimer ? "delta" : "standard") << '\n';
    if (use_delta_sternheimer)
    {
        out << "sternheimer_delta_a_block " << sternheimer_delta_a_block_mode_name(delta_a_block_mode) << '\n';
    }
    out << "delta_component_diagnostic " << (write_delta_components ? "yes" : "no") << '\n';
    out << "lcao_sos_diagnostic " << (write_lcao_sos ? "yes" : "no") << '\n';
    out << "wavefunction_diagnostic " << (write_wavefunction_diagnostic ? "yes" : "no") << '\n';
    if (write_wavefunction_diagnostic)
    {
        out << "wavefunction_diagnostic_file " << wavefunction_diagnostic_config.output_filename << '\n';
    }
    out << "abfs_channels " << num_channels << '\n';
    out << "abfs_max_channels_per_atom " << max_channels << '\n';
    out << "occupied_bands_total " << sternheimer_lcao_total_occupied_bands(response_kpoints) << '\n';
    out << "sternheimer_bands_per_k_limit " << max_bands << '\n';
    out << "sternheimer_bands_truncated " << (bands_are_truncated ? "yes" : "no") << '\n';
    out << "sternheimer_truncated_v1_diagnostic "
        << (write_truncated_v1_diagnostic ? "yes" : "no") << '\n';
    out << "sternheimer_frequency_mpi " << (use_frequency_mpi ? "yes" : "no") << '\n';
    out << "sternheimer_channel_mpi " << (PARAM.inp.sternheimer_channel_mpi ? "yes" : "no") << '\n';
    out << "sternheimer_response_slots "
        << (use_nested_response_mpi ? kpoint_groups * nfreq : nfreq) << '\n';
    out << "sternheimer_channel_replicas_per_response_slot " << response_replicas << '\n';
    out << "sternheimer_mpi_layout " << PARAM.inp.sternheimer_mpi_layout << '\n';
    out << "equation_owner_formula "
        << (PARAM.inp.sternheimer_mpi_layout == "global_equation"
                ? "occupied_frequency_channel_modulo"
                : "frequency_group_assignment")
        << '\n';
    out << "sternheimer_fd_order " << PARAM.inp.sternheimer_fd_order << '\n';
    out << "sternheimer_preconditioner "
        << (solver_options.use_fd_spectral_preconditioner ? "fd_spectral" : "none") << '\n';
    out << "sternheimer_preconditioner_regularization_Ry "
        << solver_options.fd_spectral_preconditioner_regularization << '\n';
    out << "sternheimer_kpoint_mpi " << (use_kpoint_mpi ? "yes" : "no") << '\n';
    out << "sternheimer_nested_k_frequency_mpi " << (use_nested_response_mpi ? "yes" : "no") << '\n';
    out << "sternheimer_kpoint_groups " << (use_kpoint_mpi ? kpoint_groups : 1) << '\n';
    out << "sternheimer_frequency_groups " << (use_frequency_mpi ? (use_nested_response_mpi ? nfreq
                                                                                              : GlobalV::NPROC)
                                                                  : 1)
        << '\n';
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
        const double solver_tolerance
            = positive_double_from_env(kSolverToleranceEnv, default_sternheimer_solver_tolerance());
        const int solver_max_iter = positive_int_from_env(kSolverMaxIterEnv, 300);
        const double pca_threshold = nonnegative_double_from_env(kPCAThresholdEnv, PARAM.inp.exx_pca_threshold);
        const double ccp_rmesh_times = positive_double_from_env(kCCPRmeshTimesEnv, PARAM.inp.rpa_ccp_rmesh_times);

        SternheimerABACUSSTSmokeResult result;
        result.grid_data = make_sternheimer_fd_grid(pw_basis);
        result.omega = omega;
        result.pca_threshold = pca_threshold;
        result.ccp_rmesh_times = ccp_rmesh_times;
        result.perturbation_source = "abfs_ccp";

        const SternheimerFDHamiltonian hamiltonian
            = make_sternheimer_fd_hamiltonian(
                potential, pw_basis, ucell, 0, 1.0, {0.0, 0.0, 0.0}, PARAM.inp.sternheimer_fd_order);
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

        auto dot = [&result](const SternheimerFDHamiltonian::Vector& lhs, const SternheimerFDHamiltonian::Vector& rhs) {
            return sternheimer_fd_grid_dot(lhs, rhs, result.grid_data.volume_element);
        };

        SternheimerRPA::SolverOptions solver_options;
        solver_options.max_iter = solver_max_iter;
        solver_options.residual_tol = solver_tolerance;
        solver_options.use_fd_spectral_preconditioner
            = sternheimer_environment_flag(kSpectralPreconditionerEnv, true);
        solver_options.fd_spectral_preconditioner_regularization
            = nonnegative_double_from_env(kSpectralPreconditionerRegularizationEnv, 0.0);

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
                SternheimerRPA::build_rhs_from_hartree_perturbation(perturbation_ry, states.wavefunctions[ib], rhs);
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
                channel_result.polarizability
                    = SternheimerRPA::accumulate_polarizability_grid_element(channel.potential_r,
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

void run_sternheimer_abacus_chi0_output_impl(const elecstate::Potential& potential,
                                             const ModulePW::PW_Basis& pw_basis,
                                             const UnitCell& ucell,
                                             const elecstate::ElecState& elec_state,
                                             const std::string& output_dir,
                                             const LCAO_Orbitals* lcao_orbitals,
                                             const std::vector<SternheimerLCAOOccupiedKPoint>* lcao_occupied_kpoints,
                                             const std::array<int, 3>* lcao_kmesh,
                                             const ModulePW::PW_Basis_K* siab_pw_wfc,
                                             const Structure_Factor* siab_structure_factor,
                                             const SternheimerOrbitalSet* reusable_rpa_abfs)
{
    const bool write_librpa = PARAM.inp.out_sternheimer_librpa;
    const bool write_siab = PARAM.inp.out_sternheimer_siab;
    const bool write_grid_diagnostics = PARAM.inp.sternheimer_grid_diagnostics;
    if (!write_librpa && !write_siab && !write_grid_diagnostics)
    {
        return;
    }

    const bool use_frequency_mpi = PARAM.inp.sternheimer_frequency_mpi;
    const bool use_channel_mpi = PARAM.inp.sternheimer_channel_mpi;
    const std::string mpi_layout = PARAM.inp.sternheimer_mpi_layout;
    const bool use_global_equation_mpi = mpi_layout == "global_equation";
    const int nfreq = PARAM.inp.sternheimer_nfreq;
    const char* supercell_translation_sum_raw = std::getenv(kSupercellTranslationSumEnv);
    const bool full_supercell_response
        = supercell_translation_sum_raw != nullptr
          && supercell_translation_sum_raw[0] != '\0'
          && env_is_true(kSupercellFullResponseEnv);
    int response_kpoint_count = elec_state.wg.nr;
    int requested_supercell_kpoint_groups = 1;
    if (full_supercell_response)
    {
        const auto translation_sum
            = parse_sternheimer_supercell_translation_sum(supercell_translation_sum_raw);
        response_kpoint_count = sternheimer_supercell_primitive_cell_count(translation_sum);
        requested_supercell_kpoint_groups
            = positive_int_from_env(kSupercellKPointGroupsEnv, 1);
    }
    const int response_kpoint_groups
        = PARAM.inp.sternheimer_q_index > 0
              ? sternheimer_response_kpoint_group_count(full_supercell_response,
                                                         requested_supercell_kpoint_groups,
                                                         PARAM.globalv.kpar_lcao,
                                                         response_kpoint_count)
              : 1;
    const bool use_kpoint_mpi
        = PARAM.inp.sternheimer_q_index > 0 && response_kpoint_groups > 1;
    const bool use_nested_response_mpi = use_frequency_mpi && use_kpoint_mpi;
    const bool use_parallel_response_mpi = use_frequency_mpi || use_kpoint_mpi;
    const bool use_distributed_mpi = use_parallel_response_mpi || use_channel_mpi;
    initialize_sternheimer_memory_detection();
    if (!use_distributed_mpi && GlobalV::MY_RANK != 0)
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
            if (use_distributed_mpi && GlobalV::NPROC > 1)
            {
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
#endif
            return;
        }
        out << std::setprecision(16);
        out << (write_siab ? "# ABACUS Coulomb-whitened Sternheimer target output for SIAB\n"
                           : "# ABACUS Sternheimer chi0 output for LibRPA\n");
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
                                       + (use_kpoint_mpi ? "kpoint_mpi=yes" : "kpoint_mpi=no")
                                       + ","
                                       + (use_nested_response_mpi ? "nested_mpi=yes" : "nested_mpi=no"));
        if (write_librpa && PARAM.inp.out_librpa_reader_version != 1)
        {
            throw std::runtime_error("out_sternheimer_librpa requires out_librpa_reader_version=1.");
        }
        SternheimerRPA::validate_mpi_layout(mpi_layout,
                                            use_frequency_mpi,
                                            use_channel_mpi,
                                            write_siab,
                                            write_librpa,
                                            nfreq,
                                            GlobalV::NPROC);
        if (GlobalV::NPROC != 1 && !use_parallel_response_mpi)
        {
            throw std::runtime_error(
                "Sternheimer chi0 output with multiple MPI ranks requires sternheimer_frequency_mpi=true.");
        }
        const int nested_response_slots
            = response_kpoint_groups * PARAM.inp.sternheimer_nfreq;
        if (use_nested_response_mpi
            && ((!use_channel_mpi && GlobalV::NPROC != nested_response_slots)
                || (use_channel_mpi
                    && (GlobalV::NPROC < nested_response_slots
                        || GlobalV::NPROC % nested_response_slots != 0))))
        {
            throw std::runtime_error(
                "Nested Sternheimer MPI requires one rank or an integer channel-MPI multiple per k-frequency slot.");
        }
        if (use_kpoint_mpi && !use_frequency_mpi
            && response_kpoint_groups != GlobalV::NPROC)
        {
            throw std::runtime_error(
                "Sternheimer k-point MPI without frequency MPI requires NPROC=k-point-groups.");
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
        const std::string analytic_ao_path = PARAM.inp.sternheimer_ao_potential_file;
        const bool use_analytic_ao = validate_sternheimer_molecular_coulomb(PARAM.inp);
        if (std::getenv("ABACUS_STERNHEIMER_ANALYTIC_AO_POTENTIALS") != nullptr
            || sternheimer_environment_flag("ABACUS_STERNHEIMER_MOLECULAR_COULOMB_TAIL", false))
        {
            throw std::invalid_argument(
                "Use the explicit sternheimer_molecular_coulomb and "
                "sternheimer_ao_potential_file inputs instead of experimental environment flags.");
        }
        if (use_analytic_ao
            && (!use_lcao_zero_order
                || (!GlobalC::exx_info.info_ri.files_abfs.empty() && PARAM.inp.exx_pca_threshold < 1.0)))
        {
            throw std::invalid_argument("isolated_ri requires LCAO states and either product-PCA ABFS or "
                                        "external-only ABFS (exx_pca_threshold >= 1), not a mixed basis.");
        }
        if ((lcao_orbitals == nullptr) != (lcao_occupied_kpoints == nullptr))
        {
            throw std::runtime_error("Sternheimer LCAO zero-order input is incomplete.");
        }
        if (write_siab
            && (!use_lcao_zero_order || siab_pw_wfc == nullptr || siab_structure_factor == nullptr))
        {
            throw std::runtime_error(
                "out_sternheimer_siab requires LCAO Delta-ST plus the PW FFT basis and structure factor.");
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
            validate_sternheimer_full_lcao_occupied_kpoints(*lcao_occupied_kpoints,
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
                                                      response_kpoint_groups,
                                                      chi0_start_time,
                                                      reusable_rpa_abfs);
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
                PARAM.globalv.nlocal,
                -1,
                false);
            response_kpoints
                = select_sternheimer_gamma_spin_records(*lcao_occupied_kpoints, PARAM.inp.nspin);
        }

        std::vector<int> occupied_band_counts;
        occupied_band_counts.reserve(response_kpoints.size());
        for (const SternheimerLCAOOccupiedKPoint* response_kpoint: response_kpoints)
        {
            const int response_k_index = response_kpoint == nullptr ? 0 : response_kpoint->local_k_index;
            const int occupied_count = use_lcao_zero_order
                                           ? static_cast<int>(response_kpoint->coefficients.size())
                                           : occupied_band_count(elec_state, response_k_index);
            if (occupied_count <= 0)
            {
                throw std::runtime_error("No occupied DFT bands are available for a Sternheimer response record.");
            }
            occupied_band_counts.push_back(occupied_count);
        }
        const int requested_bands = positive_int_from_env(kBandsEnv, occupied_band_counts.front());
        const int fd_num_bands = std::min(requested_bands, elec_state.ekb.nc);
        if (!use_lcao_zero_order && fd_num_bands < occupied_band_counts.front())
        {
            throw std::runtime_error("Sternheimer chi0 output requires all occupied bands; use the smoke test for "
                                     "band-limited debugging.");
        }

        const int max_dense_size = positive_int_from_env(kMaxDenseEnv, 4096);
        const int lanczos_max_subspace_size = positive_int_from_env(kLanczosSubspaceEnv, 320);
        const double solver_tolerance
            = positive_double_from_env(kSolverToleranceEnv, default_sternheimer_solver_tolerance());
        const int solver_max_iter = positive_int_from_env(kSolverMaxIterEnv, 300);
        const double pca_threshold = nonnegative_double_from_env(kPCAThresholdEnv, PARAM.inp.exx_pca_threshold);
        const double ccp_rmesh_times = positive_double_from_env(kCCPRmeshTimesEnv, PARAM.inp.rpa_ccp_rmesh_times);
        const int channel_worker_user_cap = int_from_env(kChannelMaxWorkersEnv, 0);
        if (channel_worker_user_cap < 0)
        {
            throw std::invalid_argument(std::string("Invalid non-negative integer in ") + kChannelMaxWorkersEnv + ".");
        }
        const int default_frequency_rank_shift = use_frequency_mpi && GlobalV::NPROC > 1 ? 1 : 0;
        const int frequency_rank_shift = int_from_env(kFrequencyRankShiftEnv, default_frequency_rank_shift);
        const auto frequency_assignment = [&](const int ifrequency) {
            if (use_global_equation_mpi)
            {
                SternheimerRPA::FrequencyMPIAssignment assignment;
                assignment.owns_frequency = true;
                assignment.frequency_leader_rank = 0;
                assignment.frequency_group_size = GlobalV::NPROC;
                assignment.frequency_group_local_rank = GlobalV::MY_RANK;
                return assignment;
            }
            if (!use_frequency_mpi)
            {
                return SternheimerRPA::frequency_mpi_assignment(ifrequency, nfreq, 1, 0, 0, false);
            }
            return SternheimerRPA::frequency_mpi_assignment(ifrequency,
                                                             nfreq,
                                                             GlobalV::NPROC,
                                                             GlobalV::MY_RANK,
                                                             frequency_rank_shift,
                                                             use_channel_mpi);
        };
        const int frequency_group_size = frequency_assignment(0).frequency_group_size;
#ifdef __MPI
        MPI_Comm chi0_frequency_group_communicator = MPI_COMM_NULL;
        bool owns_chi0_frequency_group_communicator = false;
        if (write_librpa && use_channel_mpi && frequency_group_size > 1)
        {
            if (use_global_equation_mpi)
            {
                chi0_frequency_group_communicator = MPI_COMM_WORLD;
            }
            else
            {
                const int frequency_group_color = GlobalV::MY_RANK / frequency_group_size;
                if (MPI_Comm_split(MPI_COMM_WORLD,
                                   frequency_group_color,
                                   GlobalV::MY_RANK,
                                   &chi0_frequency_group_communicator)
                    != MPI_SUCCESS)
                {
                    throw std::runtime_error("Failed to create the Sternheimer chi0 frequency-group communicator.");
                }
                owns_chi0_frequency_group_communicator = true;
            }
        }
#endif
        const std::string frequency_grid_file = PARAM.inp.sternheimer_frequency_grid_file;
        const bool use_frequency_grid_file = !frequency_grid_file.empty();
        const bool use_delta_sternheimer = PARAM.inp.sternheimer_delta;
        const int channel_batch_width = use_delta_sternheimer ? sternheimer_channel_batch_width() : 1;
        if (use_analytic_ao && pca_threshold != PARAM.inp.exx_pca_threshold)
        {
            throw std::invalid_argument("isolated_ri cannot override the producer auxiliary-basis PCA threshold.");
        }
        const SternheimerDeltaABlockMode delta_a_block_mode = delta_a_block_mode_from_env();
        if (use_lcao_zero_order)
        {
            if (!use_delta_sternheimer)
            {
                throw std::runtime_error(
                    "Sternheimer LCAO zero-order input currently requires sternheimer_delta=true.");
            }
            if (PARAM.inp.nspin != 1 && PARAM.inp.nspin != 2)
            {
                throw std::runtime_error("Sternheimer LCAO zero-order input currently supports Gamma-point nspin=1 or "
                                         "nspin=2 calculations.");
            }
            for (std::size_t ispin = 0; ispin != lcao_occupied_channels->size(); ++ispin)
            {
                if ((*lcao_occupied_channels)[ispin].coefficients.size()
                    != static_cast<std::size_t>(occupied_band_counts[ispin]))
                {
                    throw std::runtime_error("Sternheimer LCAO occupied coefficient count does not match occupations.");
                }
                if (lcao_virtual_source == SternheimerLCAOVirtualSource::KSBands)
                {
                    const auto& channel = (*lcao_occupied_channels)[ispin];
                    const std::size_t available_bands
                        = channel.coefficients.size() + channel.unoccupied_coefficients.size();
                    if (available_bands != static_cast<std::size_t>(PARAM.globalv.nlocal))
                    {
                        throw std::runtime_error(
                            "Sternheimer ks_bands virtual source requires nbands=nlocal so the LCAO virtual space is "
                            "complete.");
                    }
                    if (channel.unoccupied_coefficients.empty())
                    {
                        throw std::runtime_error("Sternheimer ks_bands virtual source found no unoccupied KS bands.");
                    }
                }
            }
        }

        std::vector<SternheimerRPA::TransitionEnergyWindow> spin_transition_windows;
        spin_transition_windows.reserve(response_spin_indices.size());
        bool transition_window_available = true;
        for (const int spin_index: response_spin_indices)
        {
            const std::vector<double> eigenvalues_ry = eigenvalues_ry_from_elec_state(elec_state, spin_index);
            const std::vector<double> occupations = occupations_from_elec_state(elec_state, spin_index);
            SternheimerRPA::TransitionEnergyWindow spin_window;
            if (SternheimerRPA::try_transition_energy_window_from_eigenvalues_ry(eigenvalues_ry,
                                                                                 occupations,
                                                                                 spin_window))
            {
                spin_transition_windows.push_back(spin_window);
            }
            else
            {
                transition_window_available = false;
            }
        }
        if (!use_frequency_grid_file && !transition_window_available)
        {
            throw std::runtime_error("Sternheimer GreenX minimax grid requires a positive occupied-to-empty "
                                     "transition in every response spin channel.");
        }
        SternheimerRPA::TransitionEnergyWindow transition_window;
        if (transition_window_available)
        {
            transition_window = SternheimerRPA::merge_transition_energy_windows(spin_transition_windows);
        }
        const std::string frequency_grid_source = use_frequency_grid_file ? "file" : "greenx_minimax";
        const SternheimerRPA::FrequencyGrid frequency_grid
            = use_frequency_grid_file
                  ? SternheimerRPA::read_frequency_grid_file(frequency_grid_file, nfreq)
                  : SternheimerRPA::generate_greenx_minimax_frequency_grid(nfreq,
                                                                           transition_window.emin_ha,
                                                                           transition_window.emax_ha);
        const bool transition_window_available
            = std::isfinite(transition_window.emin_ha) && std::isfinite(transition_window.emax_ha)
              && transition_window.emax_ha > transition_window.emin_ha;
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

        const std::array<int, 3> requested_response_dimensions = {PARAM.inp.sternheimer_response_nx,
                                                                  PARAM.inp.sternheimer_response_ny,
                                                                  PARAM.inp.sternheimer_response_nz};
        const bool response_grid_requested
            = PARAM.inp.sternheimer_response_ecutwfc > 0.0
              || std::any_of(requested_response_dimensions.begin(),
                             requested_response_dimensions.end(),
                             [](const int value) { return value != 0; });
        if (response_grid_requested && (!use_delta_sternheimer || !use_lcao_zero_order || write_siab))
        {
            throw std::invalid_argument(
                "An independent molecular Sternheimer response grid requires LCAO Delta-ST without SIAB output.");
        }
        SternheimerResponseGrid response_grid = make_sternheimer_response_grid(pw_basis,
                                                                                PARAM.inp.ecutwfc,
                                                                                PARAM.inp.sternheimer_response_ecutwfc,
                                                                                PARAM.inp.fft_mode,
                                                                                requested_response_dimensions,
                                                                                PARAM.inp.sternheimer_fd_order);
        const ModulePW::PW_Basis& response_pw_basis = *response_grid.basis;
        const SternheimerABACUSFDGridData grid_data
            = use_frequency_mpi ? make_sternheimer_fd_full_grid(response_pw_basis)
                                : make_sternheimer_fd_grid(response_pw_basis);
        if (GlobalV::MY_RANK == 0)
        {
            out << "sternheimer_response_ecutwfc_requested_Ry " << PARAM.inp.sternheimer_response_ecutwfc << '\n';
            out << "sternheimer_response_ecutwfc_Ry "
                << (response_grid.source == SternheimerResponseGridSource::Cutoff
                        ? PARAM.inp.sternheimer_response_ecutwfc
                        : (response_grid.source == SternheimerResponseGridSource::Pbe ? PARAM.inp.ecutwfc : 0.0))
                << '\n';
            out << "sternheimer_response_grid_requested " << requested_response_dimensions[0] << ' '
                << requested_response_dimensions[1] << ' ' << requested_response_dimensions[2] << '\n';
            out << "sternheimer_response_grid_source " << sternheimer_response_grid_source_name(response_grid.source)
                << '\n';
            out << "pbe_grid " << pw_basis.nx << ' ' << pw_basis.ny << ' ' << pw_basis.nz << " size "
                << pw_basis.nxyz << '\n';
            out << "response_grid " << grid_data.grid.nx << ' ' << grid_data.grid.ny << ' ' << grid_data.grid.nz
                << " size " << grid_data.grid.size() << " dV " << grid_data.volume_element << '\n';
        }
        SternheimerABFBuildData abfs_data
            = build_abfs_ccp_data(ucell, grid_data.grid, -1, pca_threshold, ccp_rmesh_times, write_siab);
        std::vector<SternheimerABFGridChannel>& channels = abfs_data.channels;
        const std::string abfs_source
            = sternheimer_abfs_perturbation_source(GlobalC::exx_info.info_ri.files_abfs);
        if (channels.empty())
        {
            throw std::runtime_error("No ABFS CCP perturbation channels were generated.");
        }

        const int raw_num_channels = static_cast<int>(channels.size());
        if (write_siab && GlobalC::exx_info.info_ri.files_abfs.empty())
        {
            throw std::runtime_error(
                "out_sternheimer_siab requires explicit ABFS_ORBITAL files for immutable target provenance.");
        }
        const std::string auxiliary_basis_sha256
            = write_siab && GlobalV::MY_RANK == 0
                  ? siab::sha256_unique_file_manifest(GlobalC::exx_info.info_ri.files_abfs)
                  : std::string();
        SternheimerCoulombWhitening coulomb_whitening;
        if (write_siab)
        {
            coulomb_whitening = make_sternheimer_coulomb_whitening(
                abfs_data.full_coulomb_metric,
                raw_num_channels,
                PARAM.inp.sternheimer_siab_coulomb_threshold);
        }
        const int num_channels = write_siab ? coulomb_whitening.retained_rank : raw_num_channels;
        append_chi0_progress_event("channels_ready",
                                   0,
                                   -1,
                                   -1,
                                   -1,
                                   0,
                                   nullptr,
                                   -1.0,
                                   elapsed_seconds_since(chi0_start_time),
                                   "raw_channels=" + std::to_string(raw_num_channels)
                                       + " response_channels=" + std::to_string(num_channels));

        SIABPrimitiveExportData siab_primitives;
        SternheimerSIABMemoryEstimate siab_memory_estimate;
        std::uint64_t siab_slurm_memory_bytes = 0;
        if (write_siab)
        {
            siab_primitives = build_siab_primitive_export_data(pw_basis,
                                                               *siab_structure_factor,
                                                               ucell);
            append_chi0_progress_event("siab_primitives_ready",
                                       0,
                                       -1,
                                       -1,
                                       -1,
                                       0,
                                       nullptr,
                                       -1.0,
                                       elapsed_seconds_since(chi0_start_time),
                                       "nprimitive=" + std::to_string(siab_primitives.primitive_count));
            int occupied_total = 0;
            for (const int count: occupied_band_counts)
            {
                occupied_total += count;
            }
            siab_memory_estimate
                = estimate_sternheimer_siab_dense_memory(grid_data.grid.size(),
                                                         raw_num_channels,
                                                         num_channels,
                                                         siab_primitives.primitive_count,
                                                         siab_primitives.reciprocal_count,
                                                         occupied_total,
                                                         nfreq);
            siab_slurm_memory_bytes = slurm_memory_limit_bytes();
            constexpr std::uint64_t solver_reserve_bytes = 16ULL * 1024ULL * 1024ULL * 1024ULL;
            if (siab_slurm_memory_bytes > 0
                && (siab_memory_estimate.total_bytes > siab_slurm_memory_bytes
                    || solver_reserve_bytes > siab_slurm_memory_bytes - siab_memory_estimate.total_bytes))
            {
                throw std::runtime_error(
                    "Sternheimer SIAB estimated dense allocation plus the 16-GiB solver reserve exceeds "
                    "SLURM_MEM_PER_NODE.");
            }
            if (GlobalV::MY_RANK == 0)
            {
                std::ofstream memory("STERNHEIMER_SIAB_MEMORY.dat");
                if (!memory)
                {
                    throw std::runtime_error("Failed to open Sternheimer SIAB memory diagnostic file.");
                }
                memory << "# ABACUS Sternheimer SIAB conservative dense-memory estimate\n";
                memory << "representation streamed_raw_hartree_and_serial_reciprocal_primitives_v1\n";
                memory << "coulomb_metric_bytes " << siab_memory_estimate.coulomb_metric_bytes << '\n';
                memory << "transformed_potential_bytes " << siab_memory_estimate.transformed_potential_bytes << '\n';
                memory << "channel_transform_workspace_bytes "
                       << siab_memory_estimate.channel_transform_workspace_bytes << '\n';
                memory << "reciprocal_primitive_bytes " << siab_memory_estimate.reciprocal_primitive_bytes << '\n';
                memory << "primitive_overlap_bytes " << siab_memory_estimate.primitive_overlap_bytes << '\n';
                memory << "gathered_reference_row_bytes " << siab_memory_estimate.gathered_reference_row_bytes << '\n';
                memory << "estimated_dense_total_bytes " << siab_memory_estimate.total_bytes << '\n';
                memory << "solver_reserve_bytes " << solver_reserve_bytes << '\n';
                memory << "slurm_memory_per_node_bytes " << siab_slurm_memory_bytes << '\n';
            }
            append_chi0_progress_event("siab_memory_ready",
                                       0,
                                       -1,
                                       -1,
                                       -1,
                                       0,
                                       nullptr,
                                       -1.0,
                                       elapsed_seconds_since(chi0_start_time),
                                       "dense_bytes=" + std::to_string(siab_memory_estimate.total_bytes));
        }
        if (GlobalV::MY_RANK == 0)
        {
            write_abfs_channel_diagnostic("STERNHEIMER_ABFS_CHANNELS.dat", channels);
            if (write_siab)
            {
                write_coulomb_whitening_diagnostic("STERNHEIMER_SIAB_COULOMB_WHITENING.dat",
                                                    coulomb_whitening);
            }
            if (sternheimer_abfs_diag_only_enabled())
            {
                out << "status abfs_diag_only\n";
                out << "grid " << grid_data.grid.nx << ' ' << grid_data.grid.ny << ' ' << grid_data.grid.nz << " size "
                    << grid_data.grid.size() << " dV " << grid_data.volume_element << '\n';
                out << "frequency_grid_source " << frequency_grid_source << '\n';
                out << "nfreq " << nfreq << '\n';
                out << "abfs_channels " << raw_num_channels << '\n';
                out << "response_channels " << num_channels << '\n';
                out << "abfs_source " << abfs_source << '\n';
            }
        }
        if (sternheimer_abfs_diag_only_enabled())
        {
            return;
        }

        std::vector<std::vector<double>> potentials;
        std::vector<std::vector<double>> perturbations_ry;
        if (write_librpa)
        {
            potentials = collect_channel_potentials(channels);
            // Raw ABFS Coulomb potentials are retained in Ha for M=V chi0 V output.
            perturbations_ry = scale_potentials(potentials, kHartreeToRydberg);
        }
        else
        {
            perturbations_ry = transform_channel_potentials(abfs_data, grid_data.grid, coulomb_whitening);
            scale_potentials_in_place(perturbations_ry, kHartreeToRydberg);
            if (GlobalV::MY_RANK == 0)
            {
                write_abfs_channel_diagnostic("STERNHEIMER_ABFS_CHANNELS.dat", channels);
            }
            for (SternheimerABFGridChannel& channel: channels)
            {
                channel.potential_r.clear();
                channel.potential_r.shrink_to_fit();
            }
        }

        std::vector<SternheimerDeltaGridFunction> sampled_ao_functions;
        if (use_lcao_zero_order)
        {
            sampled_ao_functions = build_lcao_candidate_grid_functions(ucell, grid_data.grid, lcao_orbitals);
            if (sampled_ao_functions.empty())
            {
                throw std::runtime_error("Sternheimer LCAO zero-order input found no sampled AO functions.");
            }
        }

        SternheimerAOPotentials analytic_ao_potentials;
        if (use_analytic_ao)
        {
            analytic_ao_potentials = read_sternheimer_ao_potentials(analytic_ao_path,
                                                                    static_cast<int>(sampled_ao_functions.size()),
                                                                    num_channels);
        }
        const std::size_t response_count = response_kpoints.size();
        std::vector<SternheimerPerturbationTensor> analytic_ao_tensors(response_count);
        double analytic_ao_max_reconstruction_error = 0.0;
        std::vector<SternheimerFDHamiltonian> hamiltonians;
        std::vector<SternheimerFDZeroOrderStates> states_by_response(response_count);
        std::vector<std::vector<SternheimerDeltaGridFunction>> lcao_projector_functions_by_response(response_count);
        std::vector<std::vector<SternheimerFDHamiltonian::Vector>> occupied_by_response(response_count);
        std::vector<std::vector<SternheimerFDHamiltonian::Vector>> occupied_projector_by_response(response_count);
        std::vector<SternheimerDeltaSubspace> delta_subspaces(response_count);
        std::vector<SternheimerDeltaFixedSubspace> delta_fixed_subspaces(response_count);
        hamiltonians.reserve(response_count);

        for (std::size_t response_index = 0; response_index != response_count; ++response_index)
        {
            const SternheimerLCAOOccupiedKPoint* response_kpoint = response_kpoints[response_index];
            const int response_k_index = response_kpoint == nullptr ? 0 : response_kpoint->local_k_index;
            const int response_spin_index = response_kpoint == nullptr ? 0 : response_kpoint->spin_index;
            const SternheimerReducedKPoint response_grid_kpoint
                = response_kpoint == nullptr ? SternheimerReducedKPoint{0.0, 0.0, 0.0}
                                              : sternheimer_lcao_grid_kpoint(*response_kpoint);
            hamiltonians.push_back([&]() {
                if (response_grid.independent)
                {
                    const auto fine_potential
                        = copy_sternheimer_full_local_potential(potential, pw_basis, response_spin_index);
                    const auto response_potential
                        = response_grid.source == SternheimerResponseGridSource::Explicit
                              ? restrict_sternheimer_real_field_rectangular(*response_grid.serial_fine_basis,
                                                                            *response_grid.serial_response_basis,
                                                                            fine_potential)
                              : restrict_sternheimer_real_field(*response_grid.serial_fine_basis,
                                                                *response_grid.serial_response_basis,
                                                                fine_potential);
                    SternheimerABACUSFDGridData response_grid_data = grid_data;
                    response_grid_data.grid.kpoint = response_grid_kpoint;
                    auto nonlocal_projector = make_sternheimer_fd_nonlocal_projector_from_unitcell(
                        ucell, response_grid_data.grid, response_grid_data.volume_element);
                    return make_sternheimer_fd_hamiltonian_from_local_potential(response_grid_data,
                                                                               response_potential,
                                                                               1.0,
                                                                               std::move(nonlocal_projector),
                                                                               PARAM.inp.sternheimer_fd_order);
                }
                return use_frequency_mpi
                    ? make_sternheimer_fd_full_hamiltonian(
                          potential,
                          pw_basis,
                          ucell,
                          response_spin_index,
                          1.0,
                          response_grid_kpoint,
                          PARAM.inp.sternheimer_fd_order)
                    : make_sternheimer_fd_hamiltonian(
                          potential,
                          pw_basis,
                          ucell,
                          response_spin_index,
                          1.0,
                          response_grid_kpoint,
                          PARAM.inp.sternheimer_fd_order);
            }());
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
                                                      fd_num_bands,
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
                if (use_analytic_ao)
                {
                    std::vector<const SternheimerFDHamiltonian::Vector*> ao_states;
                    const auto& wavefunctions = states_by_response[response_index].wavefunctions;
                    for (const auto& state: wavefunctions)
                    {
                        ao_states.push_back(&state);
                    }
                    for (const auto& state: delta_subspaces[response_index].virtual_states)
                    {
                        ao_states.push_back(&state.orbital);
                    }
                    double reconstruction_error = 0.;
                    double gram_reciprocal_condition = 0.;
                    const auto coordinates = recover_delta_grid_ao_coefficients(sampled_ao_functions,
                                                                                ao_states,
                                                                                grid_data.volume_element,
                                                                                1.e-8,
                                                                                reconstruction_error,
                                                                                &gram_reciprocal_condition);
                    analytic_ao_max_reconstruction_error
                        = std::max(analytic_ao_max_reconstruction_error, reconstruction_error);
                    const std::vector<std::vector<SternheimerFDHamiltonian::Complex>> occupied_coordinates(
                        coordinates.begin(),
                        coordinates.begin() + wavefunctions.size());
                    const std::vector<std::vector<SternheimerFDHamiltonian::Complex>> virtual_coordinates(
                        coordinates.begin() + wavefunctions.size(),
                        coordinates.end());
                    analytic_ao_tensors[response_index] = contract_sternheimer_ao_potentials(analytic_ao_potentials,
                                                                                             occupied_coordinates,
                                                                                             virtual_coordinates);
                    if (GlobalV::MY_RANK == 0)
                    {
                        out << "analytic_ao_source " << analytic_ao_path << '\n';
                        out << "analytic_ao_reconstruction_error " << reconstruction_error << '\n';
                        out << "analytic_ao_gram_reciprocal_condition " << gram_reciprocal_condition << '\n';
                        out << "analytic_ao_coupling_unit Hartree\n";
                    }
                }
                delta_fixed_subspaces[response_index]
                    = build_delta_sternheimer_fixed_subspace(occupied_projector_by_response[response_index],
                                                             delta_subspaces[response_index].virtual_states);
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
        solver_options.use_fd_spectral_preconditioner
            = sternheimer_environment_flag(kSpectralPreconditionerEnv, true);
        solver_options.fd_spectral_preconditioner_regularization
            = nonnegative_double_from_env(kSpectralPreconditionerRegularizationEnv, 0.0);

        bool all_converged = true;
        int solved_equations = 0;
        std::int64_t local_iteration_sum = 0;
        double max_solver_relative_residual = 0.0;
        double max_equation_residual_norm = 0.0;
        double max_component_reconstruction_error = 0.0;
        const std::vector<SternheimerRPA::AuxiliaryChannel> auxiliary_channels
            = write_librpa ? make_chi0_auxiliary_channels(channels)
                           : std::vector<SternheimerRPA::AuxiliaryChannel>();

        const std::size_t response_matrix_size
            = write_librpa ? static_cast<std::size_t>(num_channels) * static_cast<std::size_t>(num_channels) : 0;
        std::vector<std::vector<SternheimerRPA::Complex>> chi0_branches(static_cast<std::size_t>(nfreq));
        std::vector<std::vector<SternheimerRPA::Complex>> chi0_sos_branches;
        std::vector<std::vector<SternheimerRPA::Complex>> chi0_pulay_branches;
        std::vector<std::vector<SternheimerRPA::Complex>> chi0_qspace_branches;
        if (write_grid_diagnostics)
        {
            chi0_sos_branches.resize(static_cast<std::size_t>(nfreq));
            chi0_pulay_branches.resize(static_cast<std::size_t>(nfreq));
            chi0_qspace_branches.resize(static_cast<std::size_t>(nfreq));
        }
        std::vector<std::chrono::steady_clock::time_point> frequency_start_times(static_cast<std::size_t>(nfreq));
        for (int ifrequency = 0; ifrequency != nfreq; ++ifrequency)
        {
            const SternheimerRPA::FrequencyMPIAssignment assignment = frequency_assignment(ifrequency);
            const int owner_rank = assignment.frequency_leader_rank;
            if (!assignment.owns_frequency)
            {
                continue;
            }
            if (write_librpa)
            {
                chi0_branches[static_cast<std::size_t>(ifrequency)].assign(response_matrix_size,
                                                                           SternheimerRPA::Complex(0.0, 0.0));
                if (write_grid_diagnostics)
                {
                    chi0_sos_branches[static_cast<std::size_t>(ifrequency)].assign(
                        response_matrix_size, SternheimerRPA::Complex(0.0, 0.0));
                    chi0_pulay_branches[static_cast<std::size_t>(ifrequency)].assign(
                        response_matrix_size, SternheimerRPA::Complex(0.0, 0.0));
                    chi0_qspace_branches[static_cast<std::size_t>(ifrequency)].assign(
                        response_matrix_size, SternheimerRPA::Complex(0.0, 0.0));
                }
            }
            frequency_start_times[static_cast<std::size_t>(ifrequency)] = std::chrono::steady_clock::now();
            const double omega_ha = frequency_grid.omega_ha[static_cast<std::size_t>(ifrequency)];
            append_chi0_progress_event("frequency_start",
                                       ifrequency + 1,
                                       owner_rank,
                                       -1,
                                       -1,
                                       solved_equations,
                                       nullptr,
                                       -1.0,
                                       elapsed_seconds_since(chi0_start_time),
                                       "omega_Ha=" + std::to_string(omega_ha)
                                           + " group_size=" + std::to_string(assignment.frequency_group_size)
                                           + " group_local_rank="
                                           + std::to_string(assignment.frequency_group_local_rank));
        }

        std::vector<int> occupied_state_offsets(response_count, 0);
        int occupied_state_offset = 0;
        for (std::size_t response_index = 0; response_index != response_count; ++response_index)
        {
            occupied_state_offsets[response_index] = occupied_state_offset;
            occupied_state_offset += static_cast<int>(states_by_response[response_index].wavefunctions.size());
        }
        std::vector<siab::ReferenceRow> local_siab_rows;
        const SternheimerMemorySnapshot channel_memory = detect_sternheimer_memory_snapshot();
        SternheimerChannelWorkerPlan channel_worker_plan;
        bool channel_worker_plan_reported = false;

        for (int ifrequency = 0; ifrequency != nfreq; ++ifrequency)
        {
            const SternheimerRPA::FrequencyMPIAssignment assignment = frequency_assignment(ifrequency);
            const int owner_rank = assignment.frequency_leader_rank;
            if (!assignment.owns_frequency)
            {
                continue;
            }
            const double omega_ha = frequency_grid.omega_ha[static_cast<std::size_t>(ifrequency)];
            const double omega_ry = 2.0 * omega_ha;
            std::vector<SternheimerRPA::Complex>* chi0_branch
                = write_librpa ? &chi0_branches[static_cast<std::size_t>(ifrequency)] : nullptr;
            std::vector<SternheimerRPA::Complex>* chi0_sos_branch
                = write_grid_diagnostics ? &chi0_sos_branches[static_cast<std::size_t>(ifrequency)] : nullptr;
            std::vector<SternheimerRPA::Complex>* chi0_pulay_branch
                = write_grid_diagnostics ? &chi0_pulay_branches[static_cast<std::size_t>(ifrequency)] : nullptr;
            std::vector<SternheimerRPA::Complex>* chi0_qspace_branch
                = write_grid_diagnostics ? &chi0_qspace_branches[static_cast<std::size_t>(ifrequency)] : nullptr;

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
                const SternheimerDeltaFixedSubspace& delta_fixed_subspace = delta_fixed_subspaces[response_index];

                for (int ib = 0; ib != static_cast<int>(states.wavefunctions.size()); ++ib)
                {
                    const double occupation = elec_state.wg(spin_index, ib);
                    if (occupation <= 1.0e-8)
                    {
                        continue;
                    }

                    struct ChannelEquationResult
                    {
                        int channel_index = -1;
                        int owner_rank = -1;
                        SternheimerRPA::SolverResult solver;
                        double equation_residual_norm = 0.0;
                        bool has_siab_row = false;
                        siab::ReferenceRow siab_row;
                    };

                    std::vector<int> owned_channels;
                    std::vector<int> equation_owner_ranks;
                    owned_channels.reserve(static_cast<std::size_t>(num_channels));
                    equation_owner_ranks.reserve(static_cast<std::size_t>(num_channels));
                    for (int ichannel = 0; ichannel != num_channels; ++ichannel)
                    {
                        int equation_owner_rank = owner_rank;
                        if (use_global_equation_mpi)
                        {
                            equation_owner_rank = SternheimerRPA::global_equation_owner(
                                occupied_state_offsets[response_index] + ib,
                                ifrequency,
                                ichannel,
                                nfreq,
                                num_channels,
                                GlobalV::NPROC,
                                frequency_rank_shift);
                            if (equation_owner_rank != GlobalV::MY_RANK)
                            {
                                continue;
                            }
                        }
                        else if (use_channel_mpi)
                        {
                            const int group_owner = SternheimerRPA::channel_group_owner(
                                occupied_state_offsets[response_index] + ib,
                                ichannel,
                                num_channels,
                                assignment.frequency_group_size);
                            equation_owner_rank = assignment.frequency_leader_rank + group_owner;
                            if (group_owner != assignment.frequency_group_local_rank)
                            {
                                continue;
                            }
                        }
                        owned_channels.push_back(ichannel);
                        equation_owner_ranks.push_back(equation_owner_rank);
                    }

                    if (owned_channels.empty())
                    {
                        continue;
                    }
                    const SternheimerChannelWorkerPlan local_channel_worker_plan
                        = plan_sternheimer_owned_channel_workers(num_channels,
                                                                static_cast<int>(owned_channels.size()),
                                                                sternheimer_channel_openmp_threads(),
                                                                grid_data.grid.size(),
                                                                channel_worker_user_cap,
                                                                channel_memory,
                                                                channel_batch_width);
                    const int local_channel_batch_width = local_channel_worker_plan.channel_batch_width;
                    if (!channel_worker_plan_reported)
                    {
                        channel_worker_plan = local_channel_worker_plan;
                        channel_worker_plan_reported = true;
                        append_chi0_progress_event(
                            "channel_workers_ready",
                            0,
                            -1,
                            -1,
                            -1,
                            solved_equations,
                            nullptr,
                            -1.0,
                            elapsed_seconds_since(chi0_start_time),
                            "global_channels=" + std::to_string(num_channels)
                                + " local_channels=" + std::to_string(owned_channels.size()) + " "
                                + format_sternheimer_channel_worker_diagnostic(channel_memory,
                                                                               local_channel_worker_plan,
                                                                               grid_data.grid.size(),
                                                                               channel_worker_user_cap));
                    }

                    const auto finalize_channel_result = [&](const int local_task,
                                                             const SternheimerFDHamiltonian::Vector& delta_wavefunction,
                                                             const SternheimerRPA::SolverResult& solver,
                                                             const double equation_residual_norm,
                                                             const SternheimerDeltaPostprocessResult* delta_response) {
                        const int ichannel = owned_channels[static_cast<std::size_t>(local_task)];
                        ChannelEquationResult result;
                        result.channel_index = ichannel;
                        result.owner_rank = equation_owner_ranks[static_cast<std::size_t>(local_task)];
                        result.solver = solver;
                        result.equation_residual_norm = equation_residual_norm;
                        if (write_grid_diagnostics && delta_response != nullptr)
                        {
                            if (use_analytic_ao)
                            {
                                accumulate_sternheimer_analytic_ao_column(analytic_ao_tensors[response_index],
                                                                          ib,
                                                                          delta_response->sos_coefficients,
                                                                          occupation,
                                                                          ichannel,
                                                                          *chi0_sos_branch);
                                accumulate_sternheimer_analytic_ao_column(analytic_ao_tensors[response_index],
                                                                          ib,
                                                                          delta_response->pulay_coefficients,
                                                                          occupation,
                                                                          ichannel,
                                                                          *chi0_pulay_branch);
                            }
                            else
                            {
                                SternheimerRPA::accumulate_chi0_branch_column(potentials,
                                                                              states.wavefunctions[ib],
                                                                              delta_response->in_sos_wavefunction,
                                                                              grid_data.volume_element,
                                                                              occupation,
                                                                              ichannel,
                                                                              *chi0_sos_branch);
                                SternheimerRPA::accumulate_chi0_branch_column(potentials,
                                                                              states.wavefunctions[ib],
                                                                              delta_response->in_pulay_wavefunction,
                                                                              grid_data.volume_element,
                                                                              occupation,
                                                                              ichannel,
                                                                              *chi0_pulay_branch);
                            }
                            SternheimerRPA::accumulate_chi0_branch_column(potentials,
                                                                          states.wavefunctions[ib],
                                                                          delta_response->out_wavefunction,
                                                                          grid_data.volume_element,
                                                                          occupation,
                                                                          ichannel,
                                                                          *chi0_qspace_branch);
                        }
                        if (write_siab && delta_response != nullptr)
                        {
                            const auto& complete_response = delta_response->reconstructed_wavefunction;
                            if (complete_response.size() != static_cast<std::size_t>(grid_data.grid.size()))
                            {
                                throw std::runtime_error(
                                    "Sternheimer SIAB requires each equation owner to hold a complete response grid.");
                            }
                            result.has_siab_row = true;
                            result.siab_row.occupied_state = occupied_state_offsets[response_index] + ib;
                            result.siab_row.auxiliary_channel = ichannel;
                            result.siab_row.frequency_index = ifrequency;
                            result.siab_row.frequency_ha = omega_ha;
                            result.siab_row.occupation = occupation;
                            result.siab_row.frequency_weight
                                = frequency_grid.weights_ha[static_cast<std::size_t>(ifrequency)];
                            result.siab_row.norm = siab::norm(complete_response, grid_data.volume_element);
                            result.siab_row.q
                                = project_siab_response_to_primitives(complete_response, ucell, siab_primitives);
                        }
                        if (write_librpa)
                        {
                            if (use_analytic_ao)
                            {
                                if (delta_response == nullptr)
                                {
                                    throw std::runtime_error("Missing analytic AO Delta response.");
                                }
                                accumulate_sternheimer_analytic_ao_column(analytic_ao_tensors[response_index],
                                                                          ib,
                                                                          delta_response->coefficients,
                                                                          occupation,
                                                                          ichannel,
                                                                          *chi0_branch);
                                SternheimerRPA::accumulate_chi0_branch_column(potentials,
                                                                              states.wavefunctions[ib],
                                                                              delta_response->out_wavefunction,
                                                                              grid_data.volume_element,
                                                                              occupation,
                                                                              ichannel,
                                                                              *chi0_branch);
                            }
                            else
                            {
                                SternheimerRPA::accumulate_chi0_branch_column(potentials,
                                                                              states.wavefunctions[ib],
                                                                              delta_wavefunction,
                                                                              grid_data.volume_element,
                                                                              occupation,
                                                                              ichannel,
                                                                              *chi0_branch);
                            }
                        }
                        return result;
                    };

                    const auto get_delta_couplings = [&](const int ichannel, const std::vector<double>& perturbation) {
                        if (!use_analytic_ao)
                        {
                            return delta_sternheimer_perturbation_matrix_elements(delta_subspace.virtual_states,
                                                                                  perturbation,
                                                                                  states.wavefunctions[ib],
                                                                                  grid_data.volume_element);
                        }
                        const auto& tensor = analytic_ao_tensors[response_index];
                        std::vector<SternheimerFDHamiltonian::Complex> couplings(tensor.virtuals);
                        for (int a = 0; a < tensor.virtuals; ++a)
                        {
                            couplings[a] = kHartreeToRydberg * tensor.at(ib, a, ichannel);
                        }
                        return couplings;
                    };

                    std::vector<ChannelEquationResult> channel_results;
                    if (local_channel_batch_width == 1 || !use_delta_sternheimer)
                    {
                        channel_results = run_sternheimer_channel_tasks<ChannelEquationResult>(
                            static_cast<int>(owned_channels.size()),
                            [&](const int local_task) {
                                const int ichannel = owned_channels[static_cast<std::size_t>(local_task)];
                                const std::size_t channel_index = static_cast<std::size_t>(ichannel);
                                SternheimerFDHamiltonian::Vector rhs;
                                SternheimerRPA::build_rhs_from_hartree_perturbation(perturbations_ry[channel_index],
                                                                                    states.wavefunctions[ib],
                                                                                    rhs);
                                if (use_delta_sternheimer)
                                {
                                    const auto perturbation_matrix_elements
                                        = get_delta_couplings(ichannel, perturbations_ry[channel_index]);
                                    const SternheimerDeltaLinearResponse response
                                        = solve_delta_sternheimer_linear_response(hamiltonian,
                                                                                  delta_fixed_subspace,
                                                                                  states.eigenvalues[ib],
                                                                                  rhs,
                                                                                  delta_subspace.virtual_states,
                                                                                  perturbation_matrix_elements,
                                                                                  omega_ry,
                                                                                  grid_data.volume_element,
                                                                                  solver_options);
                                    return finalize_channel_result(local_task,
                                                                   response.response.reconstructed_wavefunction,
                                                                   response.solver,
                                                                   response.residual_norm,
                                                                   &response.response);
                                }
                                const SternheimerFDLinearResponse response
                                    = solve_sternheimer_fd_linear_response(hamiltonian,
                                                                           occupied,
                                                                           states.eigenvalues[ib],
                                                                           rhs,
                                                                           omega_ry,
                                                                           grid_data.volume_element,
                                                                           solver_options);
                                return finalize_channel_result(local_task,
                                                               response.delta_wavefunction,
                                                               response.solver,
                                                               response.residual_norm,
                                                               nullptr);
                            },
                            local_channel_worker_plan.effective_workers);
                    }
                    else
                    {
                        const std::vector<SternheimerChannelBatch> channel_batches
                            = make_sternheimer_channel_batches(static_cast<int>(owned_channels.size()),
                                                               local_channel_batch_width);
                        std::vector<std::vector<ChannelEquationResult>> grouped_results
                            = run_sternheimer_channel_tasks<std::vector<ChannelEquationResult>>(
                                static_cast<int>(channel_batches.size()),
                                [&](const int batch_task) {
                                    const SternheimerChannelBatch batch
                                        = channel_batches[static_cast<std::size_t>(batch_task)];
                                    SternheimerFDHamiltonian::Matrix rhs_batch(static_cast<std::size_t>(batch.size));
                                    std::vector<std::vector<SternheimerFDHamiltonian::Complex>>
                                        perturbation_matrix_elements(static_cast<std::size_t>(batch.size));
                                    for (int offset = 0; offset != batch.size; ++offset)
                                    {
                                        const int local_task = batch.begin + offset;
                                        const int ichannel = owned_channels[static_cast<std::size_t>(local_task)];
                                        const std::size_t channel_index = static_cast<std::size_t>(ichannel);
                                        SternheimerRPA::build_rhs_from_hartree_perturbation(
                                            perturbations_ry[channel_index],
                                            states.wavefunctions[ib],
                                            rhs_batch[static_cast<std::size_t>(offset)]);
                                        perturbation_matrix_elements[static_cast<std::size_t>(offset)]
                                            = get_delta_couplings(ichannel, perturbations_ry[channel_index]);
                                    }
                                    const std::vector<SternheimerDeltaLinearResponse> responses
                                        = solve_delta_sternheimer_linear_response_batch(hamiltonian,
                                                                                        delta_fixed_subspace,
                                                                                        states.eigenvalues[ib],
                                                                                        rhs_batch,
                                                                                        delta_subspace.virtual_states,
                                                                                        perturbation_matrix_elements,
                                                                                        omega_ry,
                                                                                        grid_data.volume_element,
                                                                                        solver_options);
                                    std::vector<ChannelEquationResult> results;
                                    results.reserve(static_cast<std::size_t>(batch.size));
                                    for (int offset = 0; offset != batch.size; ++offset)
                                    {
                                        const int local_task = batch.begin + offset;
                                        const SternheimerDeltaLinearResponse& response
                                            = responses[static_cast<std::size_t>(offset)];
                                        results.push_back(
                                            finalize_channel_result(local_task,
                                                                    response.response.reconstructed_wavefunction,
                                                                    response.solver,
                                                                    response.residual_norm,
                                                                    &response.response));
                                    }
                                    return results;
                                },
                                local_channel_worker_plan.effective_workers);
                        channel_results.reserve(owned_channels.size());
                        for (std::vector<ChannelEquationResult>& group: grouped_results)
                        {
                            for (ChannelEquationResult& result: group)
                            {
                                channel_results.push_back(std::move(result));
                            }
                        }
                    }

                    for (ChannelEquationResult& result: channel_results)
                    {
                        if (result.has_siab_row)
                        {
                            local_siab_rows.push_back(std::move(result.siab_row));
                        }
                        all_converged = all_converged && result.solver.converged;
                        ++solved_equations;
                        local_iteration_sum += result.solver.iterations;
                        max_solver_relative_residual
                            = std::max(max_solver_relative_residual, result.solver.relative_residual);
                        max_equation_residual_norm
                            = std::max(max_equation_residual_norm, result.equation_residual_norm);
                        append_chi0_progress_event("equation",
                                                   ifrequency + 1,
                                                   result.owner_rank,
                                                   ib,
                                                   result.channel_index,
                                                   solved_equations,
                                                   &result.solver,
                                                   result.equation_residual_norm,
                                                   elapsed_seconds_since(chi0_start_time),
                                                   "spin=" + std::to_string(spin_index + 1)
                                                       + " mode=" + (use_delta_sternheimer ? "delta" : "standard"));
                    }
                }
            }

        }

        std::vector<siab::ReferenceRow> global_siab_rows;
        if (write_siab)
        {
#ifdef __MPI
            global_siab_rows = siab::gather_reference_rows_to_root(local_siab_rows,
                                                                     static_cast<std::size_t>(siab_primitives.primitive_count),
                                                                     0,
                                                                     MPI_COMM_WORLD);
#else
            global_siab_rows = siab::gather_reference_rows_to_root(local_siab_rows,
                                                                     static_cast<std::size_t>(siab_primitives.primitive_count),
                                                                     0);
#endif
            if (GlobalV::MY_RANK == 0)
            {
                std::size_t occupied_total = 0;
                for (const int count: occupied_band_counts)
                {
                    occupied_total += static_cast<std::size_t>(count);
                }
                const std::size_t expected_rows
                    = occupied_total * static_cast<std::size_t>(num_channels) * static_cast<std::size_t>(nfreq);
                if (global_siab_rows.size() != expected_rows)
                {
                    throw std::runtime_error("Sternheimer SIAB global row assembly has missing or duplicate rows.");
                }
                const siab::Provenance provenance = make_siab_production_provenance(ucell,
                                                                                       auxiliary_basis_sha256,
                                                                                       frequency_grid,
                                                                                       pca_threshold,
                                                                                       coulomb_whitening);
                siab::write_v1(join_output_path(output_dir, "sternheimer_matrix.dat"),
                               grid_data.volume_element,
                               siab_primitives.blocks,
                               global_siab_rows,
                               siab_primitives.overlap_s,
                               provenance);
                GlobalV::ofs_running << " Sternheimer SIAB v1 output: "
                                     << join_output_path(output_dir, "sternheimer_matrix.dat") << std::endl;
            }
        }

        for (int ifrequency = 0; ifrequency != nfreq; ++ifrequency)
        {
            const SternheimerRPA::FrequencyMPIAssignment assignment = frequency_assignment(ifrequency);
            const int owner_rank = assignment.frequency_leader_rank;
            if (!assignment.owns_frequency)
            {
                continue;
            }
            const double omega_ha = frequency_grid.omega_ha[static_cast<std::size_t>(ifrequency)];
            if (write_librpa)
            {
#ifdef __MPI
                if (use_channel_mpi && assignment.frequency_group_size > 1)
                {
                    sternheimer_chi0::reduce_branch_to_root(chi0_branches[static_cast<std::size_t>(ifrequency)],
                                                             0,
                                                             chi0_frequency_group_communicator);
                    if (write_grid_diagnostics)
                    {
                        sternheimer_chi0::reduce_branch_to_root(
                            chi0_sos_branches[static_cast<std::size_t>(ifrequency)],
                            0,
                            chi0_frequency_group_communicator);
                        sternheimer_chi0::reduce_branch_to_root(
                            chi0_pulay_branches[static_cast<std::size_t>(ifrequency)],
                            0,
                            chi0_frequency_group_communicator);
                        sternheimer_chi0::reduce_branch_to_root(
                            chi0_qspace_branches[static_cast<std::size_t>(ifrequency)],
                            0,
                            chi0_frequency_group_communicator);
                    }
                }
#endif
                const bool writes_frequency = !use_channel_mpi || assignment.frequency_group_local_rank == 0;
                if (writes_frequency)
                {
                    const std::vector<SternheimerRPA::Complex> chi0
                        = SternheimerRPA::symmetrize_chi0_imaginary_frequency(
                            chi0_branches[static_cast<std::size_t>(ifrequency)], num_channels);
                    const SternheimerRPA::Chi0V1Metadata metadata
                        = make_chi0_v1_metadata(ucell,
                                                channels,
                                                1,
                                                ifrequency + 1,
                                                omega_ha,
                                                frequency_grid.weights_ha[static_cast<std::size_t>(ifrequency)]);
                    const std::string data_file = chi0_v1_filename(metadata.iq, metadata.ifrequency);
                    if (write_grid_diagnostics)
                    {
                        const std::vector<SternheimerRPA::Complex> chi0_sos
                            = SternheimerRPA::symmetrize_chi0_imaginary_frequency(
                                chi0_sos_branches[static_cast<std::size_t>(ifrequency)], num_channels);
                        const std::vector<SternheimerRPA::Complex> chi0_pulay
                            = SternheimerRPA::symmetrize_chi0_imaginary_frequency(
                                chi0_pulay_branches[static_cast<std::size_t>(ifrequency)], num_channels);
                        const std::vector<SternheimerRPA::Complex> chi0_qspace
                            = SternheimerRPA::symmetrize_chi0_imaginary_frequency(
                                chi0_qspace_branches[static_cast<std::size_t>(ifrequency)], num_channels);
                        const double reconstruction_error
                            = relative_component_reconstruction_error(chi0, chi0_sos, chi0_pulay, chi0_qspace);
                        if (reconstruction_error > 1.0e-10)
                        {
                            throw std::runtime_error(
                                "Sternheimer response components do not reconstruct the total response matrix.");
                        }
                        max_component_reconstruction_error
                            = std::max(max_component_reconstruction_error, reconstruction_error);
                        SternheimerRPA::write_chi0_v1_file(data_file, metadata, auxiliary_channels, chi0);
                        const std::string sos_file
                            = sternheimer_component_v1_filename("sos", metadata.iq, metadata.ifrequency, GlobalV::MY_RANK);
                        const std::string pulay_file
                            = sternheimer_component_v1_filename("pulay", metadata.iq, metadata.ifrequency, GlobalV::MY_RANK);
                        const std::string qspace_file
                            = sternheimer_component_v1_filename("qspace", metadata.iq, metadata.ifrequency, GlobalV::MY_RANK);
                        SternheimerRPA::write_chi0_v1_file(sos_file, metadata, auxiliary_channels, chi0_sos);
                        SternheimerRPA::write_chi0_v1_file(pulay_file, metadata, auxiliary_channels, chi0_pulay);
                        SternheimerRPA::write_chi0_v1_file(qspace_file, metadata, auxiliary_channels, chi0_qspace);
                    }
                    else
                    {
                        SternheimerRPA::write_chi0_v1_file(data_file, metadata, auxiliary_channels, chi0);
                    }
                    GlobalV::ofs_running << " Sternheimer chi0 v1 output: " << data_file << std::endl;
                }
            }
            append_chi0_progress_event("frequency_finish",
                                       ifrequency + 1,
                                       owner_rank,
                                       -1,
                                       -1,
                                       solved_equations,
                                       nullptr,
                                       -1.0,
                                       elapsed_seconds_since(chi0_start_time),
                                       "elapsed_freq_s="
                                           + std::to_string(elapsed_seconds_since(
                                               frequency_start_times[static_cast<std::size_t>(ifrequency)])));
        }

        const int local_solved_equations = solved_equations;
        int rank_local_equations_min = local_solved_equations;
        int rank_local_equations_max = local_solved_equations;
        std::int64_t rank_local_iterations_min = local_iteration_sum;
        std::int64_t rank_local_iterations_max = local_iteration_sum;
#ifdef __MPI
        if (use_distributed_mpi && GlobalV::NPROC > 1)
        {
            MPI_Allreduce(&local_solved_equations,
                          &rank_local_equations_min,
                          1,
                          MPI_INT,
                          MPI_MIN,
                          MPI_COMM_WORLD);
            MPI_Allreduce(&local_solved_equations,
                          &rank_local_equations_max,
                          1,
                          MPI_INT,
                          MPI_MAX,
                          MPI_COMM_WORLD);
            MPI_Allreduce(&local_iteration_sum,
                          &rank_local_iterations_min,
                          1,
                          MPI_INT64_T,
                          MPI_MIN,
                          MPI_COMM_WORLD);
            MPI_Allreduce(&local_iteration_sum,
                          &rank_local_iterations_max,
                          1,
                          MPI_INT64_T,
                          MPI_MAX,
                          MPI_COMM_WORLD);
        }
#endif
#ifdef __MPI
        if (write_grid_diagnostics && GlobalV::NPROC > 1)
        {
            MPI_Allreduce(MPI_IN_PLACE,
                          &max_component_reconstruction_error,
                          1,
                          MPI_DOUBLE,
                          MPI_MAX,
                          MPI_COMM_WORLD);
        }
#endif

#ifdef __MPI
        if (owns_chi0_frequency_group_communicator)
        {
            MPI_Comm_free(&chi0_frequency_group_communicator);
        }
#endif

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
        if (write_librpa && GlobalV::MY_RANK == 0)
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
                const int owner_rank = frequency_assignment(ifrequency).frequency_leader_rank;
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
        out << format_rpa_abfs_preorth_report(abfs_data.preorth_report);
        out << "format " << (write_siab ? "siab_v1" : "librpa_v1") << '\n';
        out << "data_files " << (write_siab ? 1 : index_entries.size()) << '\n';
        out << "index_file " << (write_siab ? "none" : "v1_sternheimer_chi0_index.dat") << '\n';
        out << "grid " << grid_data.grid.nx << ' ' << grid_data.grid.ny << ' ' << grid_data.grid.nz << " size "
            << grid_size << " dV " << grid_data.volume_element << '\n';
        out << "nfreq " << nfreq << '\n';
        out << "frequency_grid_source " << frequency_grid_source << '\n';
        if (use_frequency_grid_file)
        {
            out << "frequency_grid_file " << frequency_grid_file << '\n';
        }
        if (transition_window_available)
        {
            out << "transition_window_Ha " << transition_window.emin_ha << ' ' << transition_window.emax_ha << '\n';
        }
        else
        {
            out << "transition_window_Ha unavailable_external_grid\n";
        }
        out << "sternheimer_frequency_mpi " << (use_frequency_mpi ? "yes" : "no") << '\n';
        out << "sternheimer_channel_mpi " << (use_channel_mpi ? "yes" : "no") << '\n';
        out << "sternheimer_mpi_layout " << mpi_layout << '\n';
        out << "equation_owner_formula "
            << (use_global_equation_mpi ? "occupied_frequency_channel_modulo" : "frequency_group_assignment")
            << '\n';
        out << "frequency_group_size " << frequency_group_size << '\n';
        out << "sternheimer_channel_threads " << channel_worker_plan.effective_workers << '\n';
        out << "mpi_ranks " << GlobalV::NPROC << '\n';
        out << "frequency_rank_shift " << frequency_rank_shift << '\n';
        out << "rank_local_equations_min " << rank_local_equations_min << '\n';
        out << "rank_local_equations_max " << rank_local_equations_max << '\n';
        out << "rank_local_iterations_min " << rank_local_iterations_min << '\n';
        out << "rank_local_iterations_max " << rank_local_iterations_max << '\n';
        out << "progress_file_pattern "
            << (write_siab ? "STERNHEIMER_SIAB_PROGRESS_rank*.dat"
                           : "STERNHEIMER_CHI0_PROGRESS_rank*.dat")
            << '\n';
        out << "ifrequency omega_Ha weight_Ha omega_Ry data_file\n";
        for (const auto& entry: index_entries)
        {
            const SternheimerRPA::Chi0V1Metadata& metadata = entry.second;
            out << metadata.ifrequency << ' ' << metadata.omega << ' ' << metadata.weight << ' ' << 2.0 * metadata.omega
                << ' ' << entry.first << '\n';
        }
        out << "pca_threshold " << pca_threshold << '\n';
        out << "ccp_rmesh_times " << ccp_rmesh_times << '\n';
        out << "solver_tolerance " << solver_tolerance << '\n';
        out << "sternheimer_zero_order_source " << (use_lcao_zero_order ? "lcao_ks" : "fd_grid") << '\n';
        out << "sternheimer_lcao_virtual_source " << (use_lcao_zero_order ? "ks_bands" : "none") << '\n';
        out << "analytic_ao_coupling " << (use_analytic_ao ? "yes" : "no") << '\n';
        out << "sternheimer_molecular_coulomb " << PARAM.inp.sternheimer_molecular_coulomb << '\n';
        out << "molecular_coulomb_tail " << (use_analytic_ao ? "isolated_charge_multipole" : "none") << '\n';
        if (use_analytic_ao)
        {
            out << "analytic_ao_reconstruction_max_relative_error " << analytic_ao_max_reconstruction_error << '\n';
            out << "analytic_ao_contract rhs_and_left_projection\n";
            out << "analytic_ao_qspace re_solved_grid_rhs\n";
            out << "analytic_ao_residual hybrid_ao_qspace_equations\n";
        }
        out << "sternheimer_lcao_unoccupied_bands_per_spin";
        for (const SternheimerDeltaSubspace& delta_subspace: delta_subspaces)
        {
            out << ' ' << delta_subspace.virtual_states.size();
        }
        out << '\n';
        out << "sternheimer_response_records " << response_count << '\n';
        out << "sternheimer_response_spin_channels";
        for (const SternheimerLCAOOccupiedKPoint* response_kpoint: response_kpoints)
        {
            out << "sternheimer_lcao_virtual_source "
                << sternheimer_lcao_virtual_source_name(lcao_virtual_source) << '\n';
            out << "sternheimer_lcao_unoccupied_bands_per_spin";
            for (const int count: sternheimer_lcao_unoccupied_bands_per_spin(*lcao_occupied_channels))
            {
                out << ' ' << count;
            }
            out << '\n';
        }
        out << "sternheimer_response_spin_channels " << spin_diagnostics.size() << '\n';
        out << "sternheimer_response_spin_indices";
        for (const SpinResponseDiagnostics& diagnostics: spin_diagnostics)
        {
            out << ' ' << diagnostics.spin_index + 1;
        }
        out << '\n';
        if (spin_diagnostics.size() == 1)
        {
            out << "sternheimer_response_spin_channel " << spin_diagnostics.front().spin_index + 1 << '\n';
        }
        int total_occupied_bands = 0;
        out << "occupied_bands_per_spin";
        for (const SpinResponseDiagnostics& diagnostics: spin_diagnostics)
        {
            out << ' ' << diagnostics.occupied_bands;
            total_occupied_bands += diagnostics.occupied_bands;
        }
        out << '\n';
        out << "occupied_bands " << total_occupied_bands << '\n';
        out << "occupied_projector_dimensions_per_spin";
        for (const SpinResponseDiagnostics& diagnostics: spin_diagnostics)
        {
            out << ' ' << diagnostics.occupied_projector_dimension;
        }
        out << '\n';
        if (spin_diagnostics.size() == 1)
        {
            out << "occupied_projector_dimension " << spin_diagnostics.front().occupied_projector_dimension << '\n';
        }
        out << "abfs_channels " << raw_num_channels << '\n';
        out << "response_channels " << num_channels << '\n';
        out << "abfs_source " << abfs_source << '\n';
        if (write_siab)
        {
            out << "auxiliary_basis_sha256 " << auxiliary_basis_sha256 << '\n';
            out << "raw_auxiliary_dimension " << raw_num_channels << '\n';
            out << "whitened_auxiliary_rank " << coulomb_whitening.retained_rank << '\n';
            out << "discarded_auxiliary_rank " << coulomb_whitening.discarded_rank << '\n';
            out << "coulomb_relative_threshold " << coulomb_whitening.relative_threshold << '\n';
            out << "coulomb_max_orthonormality_error " << coulomb_whitening.max_orthonormality_error << '\n';
            out << "coulomb_transform_sha256 " << hash_coulomb_whitening_transform(coulomb_whitening) << '\n';
            out << "coulomb_whitening_diagnostic STERNHEIMER_SIAB_COULOMB_WHITENING.dat\n";
            out << "primitive_representation serial_reciprocal_pw_v1\n";
            out << "primitive_count " << siab_primitives.primitive_count << '\n';
            out << "primitive_reciprocal_count " << siab_primitives.reciprocal_count << '\n';
            out << "estimated_dense_memory_bytes " << siab_memory_estimate.total_bytes << '\n';
            out << "slurm_memory_per_node_bytes " << siab_slurm_memory_bytes << '\n';
            out << "memory_diagnostic STERNHEIMER_SIAB_MEMORY.dat\n";
            out << "target_file sternheimer_matrix.dat\n";
        }
        out << "occupied_bands " << occupied_bands_total << '\n';
        out << "occupied_bands_per_spin";
        for (const auto& occupied: occupied_by_response)
        {
            out << ' ' << occupied.size();
        }
        out << '\n';
        out << "occupied_projector_dimension " << occupied_projector_total << '\n';
        out << "occupied_projector_dimensions_per_spin";
        for (const auto& occupied_projector: occupied_projector_by_response)
        {
            out << ' ' << occupied_projector.size();
        }
        out << '\n';
        out << "abfs_channels " << num_channels << '\n';
        out << "abfs_source " << abfs_source << '\n';
        out << "sternheimer_delta " << (use_delta_sternheimer ? "yes" : "no") << '\n';
        out << "sternheimer_grid_diagnostics " << (write_grid_diagnostics ? "yes" : "no") << '\n';
        out << "sternheimer_fd_order " << PARAM.inp.sternheimer_fd_order << '\n';
        out << "sternheimer_preconditioner "
            << (solver_options.use_fd_spectral_preconditioner ? "fd_spectral" : "none") << '\n';
        out << "sternheimer_preconditioner_regularization_Ry "
            << solver_options.fd_spectral_preconditioner_regularization << '\n';
        if (write_grid_diagnostics)
        {
            out << "sternheimer_component_reconstruction_error_max " << max_component_reconstruction_error << '\n';
            out << "sternheimer_component_file_pattern "
                   "v1_sternheimer_component_{sos,pulay,qspace}_iq_<iq>_ifreq_<ifreq>_rank<rank>.dat\n";
            out << "sternheimer_grid_matrix_file_pattern STERNHEIMER_DELTA_GRID_MATRICES_spin_<spin>.dat\n";
            out << "sternheimer_perturbation_file_pattern STERNHEIMER_DELTA_PERTURBATION_spin_<spin>.dat\n";
            out << "sternheimer_perturbation_unit Rydberg\n";
        }
        if (use_delta_sternheimer)
        {
            out << "sternheimer_delta_backend " << sternheimer_delta_a_block_mode_name(delta_a_block_mode) << '\n';
            out << "sternheimer_delta_max_states " << PARAM.inp.sternheimer_delta_max_states << '\n';
            out << "sternheimer_delta_norm_tol " << PARAM.inp.sternheimer_delta_norm_tol << '\n';
            out << "sternheimer_delta_virtual_states_per_spin";
            for (const SpinResponseDiagnostics& diagnostics: spin_diagnostics)
            {
                out << ' ' << diagnostics.delta_virtual_states;
            }
            out << '\n';
            out << "sternheimer_delta_accepted_candidates_per_spin";
            for (const SpinResponseDiagnostics& diagnostics: spin_diagnostics)
            {
                out << ' ' << diagnostics.delta_accepted_candidates;
            }
            out << '\n';
            out << "sternheimer_delta_discarded_candidates_per_spin";
            for (const SpinResponseDiagnostics& diagnostics: spin_diagnostics)
            {
                out << ' ' << diagnostics.delta_discarded_candidates;
            }
            out << '\n';
            if (spin_diagnostics.size() == 1)
            {
                out << "sternheimer_delta_virtual_states " << spin_diagnostics.front().delta_virtual_states << '\n';
                out << "sternheimer_delta_accepted_candidates " << spin_diagnostics.front().delta_accepted_candidates
                    << '\n';
                out << "sternheimer_delta_discarded_candidates " << spin_diagnostics.front().delta_discarded_candidates
                    << '\n';
            }
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
        const bool invalid_excitation_spectrum
            = dynamic_cast<const SternheimerExcitationError*>(&error) != nullptr;
        const bool fatal_spectrum_failure = invalid_excitation_spectrum
            || env_is_true("ABACUS_STERNHEIMER_WEAK_RESPONSE_AUDIT")
            || env_is_true("ABACUS_STERNHEIMER_WEAK_MATRIX_AUDIT")
            || env_is_true("ABACUS_STERNHEIMER_SPECTRUM_AUDIT")
            || env_is_true("ABACUS_STERNHEIMER_GALERKIN_AUDIT");
        if (fatal_spectrum_failure)
        {
            out.flush();
            std::cerr << "Fatal Sternheimer spectrum or audit: " << error.what() << std::endl;
        }
#ifdef __MPI
        if ((use_parallel_response_mpi || fatal_spectrum_failure) && GlobalV::NPROC > 1)
        {
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
#endif
        if (fatal_spectrum_failure)
        {
            std::exit(EXIT_FAILURE);
        }
    }
}

void run_sternheimer_abacus_chi0_output(const elecstate::Potential& potential,
                                        const ModulePW::PW_Basis& pw_basis,
                                        const UnitCell& ucell,
                                        const elecstate::ElecState& elec_state,
                                        const std::string& output_dir)
{
    run_sternheimer_abacus_chi0_output_impl(potential,
                                            pw_basis,
                                            ucell,
                                            elec_state,
                                            output_dir,
                                            nullptr,
                                            nullptr,
                                            nullptr,
                                            nullptr,
                                            nullptr,
                                            nullptr);
}

void run_sternheimer_abacus_lcao_chi0_output(const elecstate::Potential& potential,
                                             const ModulePW::PW_Basis& pw_basis,
                                             const UnitCell& ucell,
                                             const elecstate::ElecState& elec_state,
                                             const LCAO_Orbitals& orbitals,
                                             const std::vector<SternheimerLCAOOccupiedKPoint>& occupied_kpoints,
                                             const std::array<int, 3>& kmesh,
                                             const ModulePW::PW_Basis_K* pw_wfc,
                                             const Structure_Factor* structure_factor,
                                             const std::string& output_dir,
                                             const SternheimerOrbitalSet* reusable_rpa_abfs)
{
    run_sternheimer_abacus_chi0_output_impl(potential,
                                            pw_basis,
                                            ucell,
                                            elec_state,
                                            output_dir,
                                            &orbitals,
                                            &occupied_kpoints,
                                            &kmesh,
                                            pw_wfc,
                                            structure_factor,
                                            reusable_rpa_abfs);
}

} // namespace ModuleRI
