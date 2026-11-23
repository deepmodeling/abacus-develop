/** @file nacv_tools.cpp
 * @brief Implement NAC result assembly and eigenstate output.
 */
#include "nacv_tools.h"
#include "source_base/parallel_cell.h"
#include "source_base/timer.h"
#include "utils/lr_util.h"
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace LR
{
NacvTimerScope::NacvTimerScope(const char *category, const char *name) : category_(category), name_(name)
{
    ModuleBase::timer::start(category_, name_);
}
NacvTimerScope::~NacvTimerScope() { ModuleBase::timer::end(category_, name_); }
} // namespace LR

namespace LR
{
void validate_nac_surfaces(const std::vector<int> &pairs, int state_count)
{
    for (const int state : pairs)
    {
        if (state > state_count)
        {
            throw std::invalid_argument("lr_nac_pairs exceeds the resolved number of LR states");
        }
    }
}

void validate_nac_window(int occupied, int virtuals, int basis)
{
    if (occupied <= 0 || virtuals <= 0 || basis <= occupied || virtuals != basis - occupied)
    {
        throw std::invalid_argument("Analytic NAC requires the complete occupied/virtual orbital window");
    }
}

void validate_nac_occupations(
    const ModuleBase::matrix &weights,
    int occupied,
    int basis,
    int spin_channels)
{
    if ((spin_channels != 1 && spin_channels != 2) || occupied <= 0 || occupied >= basis ||
        weights.nr != spin_channels || weights.nc != basis)
    {
        throw std::invalid_argument("Analytic NAC requires complete closed-shell Gamma occupation weights");
    }
    const double full_weight = 2.0 / spin_channels;
    const double roundoff = 64.0 * std::numeric_limits<double>::epsilon() * full_weight;
    for (int spin = 0; spin < spin_channels; ++spin)
    {
        for (int band = 0; band < basis; ++band)
        {
            double expected = 0.0;
            if (band < occupied)
            {
                expected = full_weight;
            }
            const double actual = weights(spin, band);
            const double occupation_error = actual - expected;
            if (!std::isfinite(actual) || std::abs(occupation_error) > roundoff)
            {
                throw std::invalid_argument(
                    "Analytic auxiliary TDA NAC requires integer closed-shell occupations; "
                    "fractional-occupation response is not implemented");
            }
        }
    }
}

void validate_nac_occupations(const ModuleBase::matrix &weights, const std::vector<int> &occupied, int basis)
{
    if (occupied.size() != 2 || weights.nr != 2 || weights.nc != basis)
    {
        throw std::invalid_argument("Unrestricted analytic NAC requires two complete spin occupation rows");
    }
    const double roundoff = 64.0 * std::numeric_limits<double>::epsilon();
    for (int spin = 0; spin < 2; ++spin)
    {
        if (occupied[spin] <= 0 || occupied[spin] >= basis)
        {
            throw std::invalid_argument("Unrestricted analytic NAC requires occupied and virtual bands in both spins");
        }
        for (int band = 0; band < basis; ++band)
        {
            double expected = 0.0;
            if (band < occupied[spin])
            {
                expected = 1.0;
            }
            const double actual = weights(spin, band);
            const double occupation_error = actual - expected;
            if (!std::isfinite(actual) || std::abs(occupation_error) > roundoff)
            {
                throw std::invalid_argument("Unrestricted analytic NAC requires integer spin occupations");
            }
        }
    }
}

void validate_nac_functionals(const std::string &reference_xc, const std::string &response_xc)
{
    const std::string reference = LR_Util::canonical_xc_name(reference_xc);
    const std::string response = LR_Util::canonical_xc_name(response_xc);
    const bool reference_supported = LR_Util::has_local_xc(reference) || reference == "hf" || reference == "none";
    const bool response_supported = LR_Util::has_local_xc(response) || response == "hf" || response == "rpa";
    if (!reference_supported || !response_supported)
    {
        throw std::invalid_argument("Analytic NAC requires implemented KS orbital-response and LR XC derivatives");
    }
}

double nacv_gap_resolution(
    double bra_energy,
    double ket_energy,
    double bra_residual,
    double ket_residual)
{
    const double bra_magnitude = std::abs(bra_energy);
    const double ket_magnitude = std::abs(ket_energy);
    const double energy_scale = std::max(bra_magnitude, ket_magnitude);
    const double roundoff_resolution = 16.0 * std::numeric_limits<double>::epsilon() * energy_scale;
    const double resolution = bra_residual + ket_residual + roundoff_resolution;
    return resolution;
}

/**
 * @brief Form raw/ETF numerators and divide only electronically resolved gaps.
 * @param[in,out] result Pair gap in Ry and AO connection;output matrices are reset.
 * @param[in] force_correction Pair force contraction, atom-major nat-by-3 in Ry/bohr.
 * @param[in] resolution Nonnegative gap-resolution bound in Ry.
 * @throw std::runtime_error For a nonfinite numerator.
 * @throw std::overflow_error For a nonfinite resolved coupling.
 * @details Numerators retain Ry/bohr; resolved vectors have units bohr^-1.
 * Unresolved divided vectors remain zero and their finite numerators are retained.
 */
void assemble_nacv_result(NacvResult &result, const ModuleBase::matrix &force_correction, double resolution)
{
    const double gap = result.gap_rydberg;
    const int atoms = force_correction.nr;
    result.raw.create(atoms, 3);
    result.etf.create(atoms, 3);
    result.numerator_raw.create(atoms, 3);
    result.numerator_etf.create(atoms, 3);
    result.gap_resolved = std::abs(gap) > resolution;
    for (int atom = 0; atom < force_correction.nr; ++atom)
    {
        for (int axis = 0; axis < 3; ++axis)
        {
            const double etf_numerator = -force_correction(atom, axis);
            const double connection = result.ao_connection(atom, axis);
            const double raw_numerator = etf_numerator + gap * connection;
            if (!std::isfinite(connection) || !std::isfinite(etf_numerator) || !std::isfinite(raw_numerator))
            {
                throw std::runtime_error("Analytic auxiliary TDA NAC produced a nonfinite numerator");
            }
            result.numerator_etf(atom, axis) = etf_numerator;
            result.numerator_raw(atom, axis) = raw_numerator;
            if (result.gap_resolved)
            {
                result.etf(atom, axis) = etf_numerator / gap;
                result.raw(atom, axis) = raw_numerator / gap;
                if (!std::isfinite(result.etf(atom, axis)) || !std::isfinite(result.raw(atom, axis)))
                {
                    throw std::overflow_error("Resolved analytic NAC division produced a nonfinite vector");
                }
            }
        }
    }
}

void print_nac_eigenpairs(
    std::ostream &output,
    int ionic_step,
    const std::string &spin_name,
    const double *energies_ry,
    const std::vector<double> &residuals_ry)
{
    const std::streamsize precision = output.precision();
    const std::ios::fmtflags flags = output.flags();
    output << std::scientific << std::setprecision(17);
    for (std::size_t state = 0; state < residuals_ry.size(); ++state)
    {
        output << "LR_EIGENPAIR " << ionic_step << ' ' << spin_name << ' ' << state + 1 << ' '
               << 0.5 * energies_ry[state] << ' ' << 0.5 * residuals_ry[state] << '\n';
    }
    output.flags(flags);
    output.precision(precision);
    if (!output.good())
    {
        throw std::runtime_error("Could not write LR eigenpair diagnostics");
    }
}

void reverse_nac_pair(NacvResult &result)
{
    result.gap_rydberg = -result.gap_rydberg;
    result.raw *= -1.0;
    result.etf *= -1.0;
    result.ao_connection *= -1.0;
}

void print_nac_result(
    std::ostream &output,
    const NacvResult &result,
    int ionic_step,
    int bra_state,
    int ket_state,
    const std::string &spin_name)
{
    const std::streamsize precision = output.precision();
    const std::ios::fmtflags flags = output.flags();
    output << std::scientific << std::setprecision(17);
    const double gap_hartree = 0.5 * result.gap_rydberg;
    output << "NACV_BEGIN " << ionic_step << ' ' << bra_state << ' ' << ket_state << ' ' << spin_name << ' '
           << result.gap_resolved << ' ' << gap_hartree << '\n';
    output << "NACV_DEFINITION auxiliary_TDA_determinant raw_and_ETF\n";
    for (int atom = 0; atom < result.raw.nr; ++atom)
    {
        for (int axis = 0; axis < 3; ++axis)
        {
            output << "NACV_COMPONENT " << atom << ' ' << axis << ' ' << result.raw(atom, axis) << ' '
                   << result.etf(atom, axis) << ' ' << result.ao_connection(atom, axis) << ' '
                   << 0.5 * result.numerator_raw(atom, axis) << ' ' << 0.5 * result.numerator_etf(atom, axis) << '\n';
        }
    }
    output << "NACV_END\n";
    output.flags(flags);
    output.precision(precision);
    if (!output.good())
    {
        throw std::runtime_error("Could not write the analytic NAC record");
    }
}
} // namespace LR

namespace LR
{
namespace
{
/**
 * @brief Serialize unfiltered state rows and confirm successful file completion.
 * @param[in] filename Output path, opened with truncation.
 * @param[in] values Row-major real or complex values in the existing LR convention.
 * @param[in] rows Number of states, or one for the excitation-energy row.
 * @param[in] columns Number of entries per row; may be zero on empty MPI ranks.
 * @return Absolute path of the successfully closed file.
 * @throw std::runtime_error For an open, write, close or path-resolution failure.
 */
template <typename T>
std::string write_nac_state_file(
    const std::string &filename,
    const T *values,
    int rows,
    int columns)
{
    std::ofstream file(filename.c_str());
    if (!file.is_open())
    {
        throw std::runtime_error("Could not open LR state output: " + filename);
    }
    file << std::scientific << std::setprecision(17);
    for (int row = 0; row < rows; ++row)
    {
        for (int column = 0; column < columns; ++column)
        {
            const int index = row * columns + column;
            file << values[index] << ' ';
        }
        file << '\n';
    }
    // close() flushes buffered output and retains write/flush/close failures.
    file.close();
    if (file.fail())
    {
        throw std::runtime_error("Could not complete LR state output: " + filename);
    }
#ifdef _WIN32
    char *absolute_path = ::_fullpath(nullptr, filename.c_str(), 0);
#else
    char *absolute_path = ::realpath(filename.c_str(), nullptr);
#endif
    if (absolute_path == nullptr)
    {
        throw std::runtime_error("Could not resolve LR state output path: " + filename);
    }
    const std::unique_ptr<char, decltype(&std::free)> resolved_path(absolute_path, &std::free);
    return std::string(resolved_path.get());
}
} // namespace

template <typename T>
void write_nac_states(
    const std::string &directory,
    const std::string &label,
    const std::string &suffix,
    int rank,
    const double *energies,
    const T *amplitudes,
    int local_size,
    int state_count,
    std::ostream &output)
{
    std::string confirmations;
    std::string failure;
    try
    {
        if (rank == 0)
        {
            const std::string filename = directory + "Excitation_Energy_" + label + suffix + ".dat";
            const std::string written = write_nac_state_file(
                filename,
                energies,
                1,
                state_count);
            confirmations += "Write file: " + written + "\n";
        }
        const int output_rank = rank + 1;
        const std::string filename =
            directory + "Excitation_Amplitude_" + label + suffix + "_" + std::to_string(output_rank) + ".dat";
        const std::string written = write_nac_state_file(
            filename,
            amplitudes,
            state_count,
            local_size);
        confirmations += "Write file: " + written + "\n";
    }
    catch (const std::runtime_error &error)
    {
        failure = error.what();
    }
    const auto domain = ModuleBase::world_comm_domain();
    const int local_failure = !failure.empty();
    if (domain.max(local_failure) != 0)
    {
        if (!local_failure)
        {
            failure = "LR state output failed on another MPI rank";
        }
        throw std::runtime_error(failure);
    }
    const int maximum_count = std::numeric_limits<int>::max();
    const std::string gathered = domain.gather_text(confirmations, maximum_count);
    output << gathered;
}
template void write_nac_states<double>(
    const std::string &,
    const std::string &,
    const std::string &,
    int,
    const double *,
    const double *,
    int,
    int,
    std::ostream &);
template void write_nac_states<std::complex<double>>(
    const std::string &,
    const std::string &,
    const std::string &,
    int,
    const double *,
    const std::complex<double> *,
    int,
    int,
    std::ostream &);
} // namespace LR
