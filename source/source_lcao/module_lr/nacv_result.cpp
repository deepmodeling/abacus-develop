#include "nacv_cs.h"
#include "nacv_workflow.h"
#include "utils/xc_kernel_names.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <ostream>
#include <stdexcept>

namespace LR
{
namespace
{
/**
 * @brief Check response availability, optionally deferring unresolved INPUT defaults.
 * @param[in] reference_xc Requested or resolved KS functional name.
 * @param[in] response_xc Requested or resolved LR kernel name.
 * @param[in] allow_defaults Whether default names may await pseudopotential resolution.
 * @throw std::invalid_argument For a known unavailable KS or LR response.
 */
void validate_functional_names(const std::string &reference_xc, const std::string &response_xc, bool allow_defaults)
{
    const std::string reference = LR_Util::canonical_xc_name(reference_xc);
    const std::string response = LR_Util::canonical_xc_name(response_xc);
    const bool reference_default = allow_defaults && reference == "default";
    const bool response_default = allow_defaults && response == "default";
    const bool reference_supported =
        reference_default || LR_Util::has_local_xc(reference) || reference == "hf" || reference == "none";
    const bool response_supported =
        response_default || LR_Util::has_local_xc(response) || response == "hf" || response == "rpa";
    if (!reference_supported || !response_supported)
    {
        throw std::invalid_argument("Analytic NAC requires implemented KS orbital-response and LR XC derivatives");
    }
}
} // namespace

bool nac_eigenpairs_requested(const NacvRequest &request, bool open_shell)
{
    return !request.pairs.empty() || (!open_shell && request.output_wavefunctions);
}

void validate_nac_reference(
    const NacvRequest &request,
    bool open_shell,
    bool converged,
    bool data_ready,
    const ModuleBase::matrix &weights,
    const std::vector<int> &occupied,
    int basis,
    int spin_channels)
{
    if (request.pairs.empty())
    {
        return;
    }
    if (!request.embedded_ks || !request.gamma_only || !request.force_response || !converged || !data_ready)
    {
        throw std::logic_error("Analytic NAC requires converged KS/LR data and the force-response window");
    }
    if (request.kpar != 1)
    {
        throw std::logic_error("Analytic NAC requires a real Gamma LR channel and kpar=1");
    }
    if (open_shell)
    {
        for (int spin = 0; spin < 2; ++spin)
        {
            validate_nac_window(occupied.at(spin), basis - occupied.at(spin), basis);
        }
        validate_nac_occupations(weights, occupied, basis);
    }
    else
    {
        const int occupied_bands = occupied.at(0);
        validate_nac_occupations(
            weights,
            occupied_bands,
            basis,
            spin_channels);
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
    validate_functional_names(reference_xc, response_xc, false);
}

void validate_nac_xc_request(const std::string &reference_xc, const std::string &response_xc)
{
    validate_functional_names(reference_xc, response_xc, true);
}

void validate_nac_window(int occupied, int virtuals, int basis)
{
    if (occupied <= 0 || virtuals <= 0 || basis <= occupied || virtuals != basis - occupied)
    {
        throw std::invalid_argument("Analytic NAC requires the complete occupied/virtual orbital window");
    }
}

double nacv_gap_resolution(
    double bra_energy,
    double ket_energy,
    double bra_residual,
    double ket_residual)
{
    if (!std::isfinite(bra_energy) || !std::isfinite(ket_energy) || !std::isfinite(bra_residual) ||
        !std::isfinite(ket_residual) || bra_residual < 0.0 || ket_residual < 0.0)
    {
        throw std::invalid_argument("NAC gap resolution requires finite energies and nonnegative residual norms");
    }
    const double bra_magnitude = std::abs(bra_energy);
    const double ket_magnitude = std::abs(ket_energy);
    const double energy_scale = std::max(bra_magnitude, ket_magnitude);
    const double roundoff_resolution = 16.0 * std::numeric_limits<double>::epsilon() * energy_scale;
    const double resolution = bra_residual + ket_residual + roundoff_resolution;
    if (!std::isfinite(resolution))
    {
        throw std::overflow_error("NAC gap resolution overflowed");
    }
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
    if (!std::isfinite(gap) || !std::isfinite(resolution) || resolution < 0.0 || force_correction.nc != 3 ||
        atoms <= 0 || result.ao_connection.nr != atoms || result.ao_connection.nc != 3)
    {
        throw std::invalid_argument("Invalid NAC result dimensions or gap resolution");
    }
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
    if (ionic_step < 1 || energies_ry == nullptr || residuals_ry.empty())
    {
        throw std::invalid_argument("LR residual output requires a one-based geometry and eigenpairs");
    }
    for (std::size_t state = 0; state < residuals_ry.size(); ++state)
    {
        if (!std::isfinite(energies_ry[state]) || !std::isfinite(residuals_ry[state]) || residuals_ry[state] < 0.0)
        {
            throw std::invalid_argument("LR residual output requires finite eigenpair diagnostics");
        }
    }
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
    if (ionic_step < 1 || bra_state < 0 || ket_state < 0 || bra_state == ket_state)
    {
        throw std::invalid_argument("NAC output requires a one-based geometry and distinct nonnegative states");
    }
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
