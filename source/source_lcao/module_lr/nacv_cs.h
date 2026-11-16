#ifndef ABACUS_SOURCE_LCAO_MODULE_LR_NACV_CS_H
#define ABACUS_SOURCE_LCAO_MODULE_LR_NACV_CS_H

#include "source_base/matrix.h"

#include <iosfwd>
#include <string>
#include <vector>

class TwoCenterBundle;
namespace module_dm
{
template <typename TK, typename TR> class DensityMatrix;
}
namespace LR
{
template <typename T> struct GradientInputs;
struct NacvResponseControls;
template <typename T> class LR_Force;

/**
 * @brief One real Gamma-point auxiliary-TDA derivative-coupling pair.
 * @details Matrices are atom-major nat x 3. The raw vector differentiates the
 * spin-adapted or spin-unrestricted auxiliary TDA determinant wavefunction. The ETF
 * vector omits the antisymmetric moving-AO connection. Neither is a full
 * quadratic-response TDDFT coupling. Unresolved gaps retain their numerators
 * but do not produce divided vectors; no gap clipping is applied.
 */
struct NacvResult
{
    ModuleBase::matrix raw;
    ModuleBase::matrix etf;
    ModuleBase::matrix ao_connection;
    ModuleBase::matrix numerator_raw;
    ModuleBase::matrix numerator_etf;
    double gap_rydberg = 0.0;
    bool gap_resolved = false;
};

/**
 * @brief Require the complete occupied/virtual orbital window for analytic NAC.
 * @param[in] occupied Number of occupied spatial orbitals in the LR window.
 * @param[in] virtuals Number of virtual spatial orbitals in the LR window.
 * @param[in] basis Number of spatial atomic orbitals and canonical MO columns.
 * @throw std::invalid_argument For empty, inconsistent or truncated windows.
 * @details Zero-padding amplitudes cannot restore omitted orbital-response terms.
 */
void validate_nac_window(int occupied, int virtuals, int basis);

/**
 * @brief Require an integer closed-shell determinant rather than fractional occupations.
 * @param[in] weights Gamma occupation weights, spin-major spin_channels by basis.
 * @param[in] occupied Number of fully occupied spatial bands in every spin channel.
 * @param[in] basis Complete spatial orbital count.
 * @param[in] spin_channels One for spin-summed weights or two for separate spins.
 * @throw std::invalid_argument For fractional, nonfinite or inconsistent weights.
 * @details Full weights are 2/spin_channels; only floating-point roundoff is allowed.
 */
void validate_nac_occupations(
    const ModuleBase::matrix &weights,
    int occupied,
    int basis,
    int spin_channels);

/**
 * @brief Require two integer occupied spin spaces for unrestricted Gamma NAC.
 * @param[in] weights Separate alpha/beta weights, two rows by basis columns.
 * @param[in] occupied Occupied alpha and beta counts, each strictly between zero and basis.
 * @param[in] basis Complete canonical MO count in each spin channel.
 * @throw std::invalid_argument For fractional, inconsistent or empty spin windows.
 * @details Each occupied spin orbital has weight one. No occupation normalization is applied.
 */
void validate_nac_occupations(const ModuleBase::matrix &weights, const std::vector<int> &occupied, int basis);

/**
 * @brief Reject unavailable ground-state or LR response kernels before analytic NAC.
 * @param[in] reference_xc Resolved KS functional name; aliases are accepted.
 * @param[in] response_xc Resolved LR kernel name; RPA denotes Hartree-only response.
 * @throw std::invalid_argument For unsupported orbital-response or LR derivatives.
 * @details A KS-recognized name is not evidence of an implemented analytic response.
 */
void validate_nac_functionals(const std::string &reference_xc, const std::string &response_xc);

/**
 * @brief Bound the resolution of a pair gap using certified eigenpair residuals.
 * @param[in] bra_energy Excitation energy in Ry; zero for the ground state.
 * @param[in] ket_energy Excitation energy in Ry.
 * @param[in] bra_residual Normalized operator residual norm in Ry; zero for S0.
 * @param[in] ket_residual Normalized operator residual norm in Ry.
 * @return Sum of residual norms and a floating-point roundoff bound, in Ry.
 * @pre Energies and nonnegative residuals come from certified LR eigenpairs.
 */
double nacv_gap_resolution(
    double bra_energy,
    double ket_energy,
    double bra_residual,
    double ket_residual);

/**
 * @brief Evaluate a real-Gamma auxiliary-TDA analytic derivative coupling.
 * @param[in] inputs Explicit ground-state geometry, orbital layouts, kernel and log.
 * @param[in,out] force_terms Existing analytic force contractions for this geometry.
 * @param[in] ground_density Ground density in the force convention, two channels for unrestricted orbitals.
 * @param[in] two_center_bundle Current analytic two-center integrals.
 * @param[in] spin_channel Closed-shell adaptation index, or zero for the concatenated up/down channel.
 * @param[in] bra_ground True for S0-to-excited coupling; otherwise an excited pair.
 * @param[in] bra Local widened bra amplitudes; ignored for S0. May be null only
 * when bra_ground or inputs.nloc is zero.
 * @param[in] ket Local widened ket amplitudes in the same canonical MO gauge.
 * @param[in] bra_energy_ry Excitation energy of the bra, zero for S0, in Rydberg.
 * @param[in] ket_energy_ry Positive ket excitation energy in Rydberg.
 * @param[in] gap_resolution_ry Nonnegative eigenpair-resolution bound in Rydberg.
 * @param[in] response_solver Existing GS Hessian solver selection.
 * @param[in] controls INPUT-validated NAC response convergence settings.
 * @pre The geometry-local NAC entry has validated the reference and selected roots.
 * @return Owned raw/ETF vectors in bohr^-1, numerators in Ry/bohr and diagnostics.
 * @throw std::invalid_argument For unsupported geometry, spin, energies or layouts.
 * @details Implements direct differentiation of the auxiliary determinant expansion
 * (Ou et al., JCP 141, 024114, 2014). EE uses common-average-energy quadratic
 * polarization and a separate antisymmetric AO connection; GE uses the GS A+B
 * adjoint plus its symmetric overlap correction. Independent overlap-FD validation
 * is required before treating a newly supported functional as production-ready.
 */
NacvResult evaluate_aux_tda_nac(
    const GradientInputs<double> &inputs,
    LR_Force<double> &force_terms,
    const module_dm::DensityMatrix<double, double> &ground_density,
    const TwoCenterBundle &two_center_bundle,
    int spin_channel,
    bool bra_ground,
    const double *bra,
    const double *ket,
    double bra_energy_ry,
    double ket_energy_ry,
    double gap_resolution_ry,
    const std::string &response_solver,
    const NacvResponseControls &controls);

/**
 * @brief Assemble raw/ETF numerators and divide only an electronically resolved gap.
 * @param[in,out] result Gap in Ry and AO connection in bohr^-1;owns returned matrices.
 * @param[in] force_correction Pair force in Ry/bohr,atom-major nat-by-3.
 * @param[in] resolution Nonnegative residual-based gap resolution in Ry.
 * @pre Matrices use the prepared nat-by-3 layout and a certified nonnegative resolution.
 * @throw std::runtime_error For a nonfinite connection or numerator.
 * @throw std::overflow_error For a nonfinite divided coupling.
 * @details Reinitializes output vectors;unresolved gaps retain numerators without division.
 */
void assemble_nacv_result(NacvResult &result, const ModuleBase::matrix &force_correction, double resolution);

/**
 * @brief Log operator-reapplied LR residuals for coupling and overlap validation.
 * @param[in,out] output Open running log on the writer rank;formatting is restored.
 * @param[in] ionic_step One-based geometry index.
 * @param[in] spin_name LR adaptation label, or updown for the unrestricted channel.
 * @param[in] energies_ry Excitation energies in Ry;one entry per residual.
 * @param[in] residuals_ry Normalized global eigenpair residual norms in Ry.
 * @details Writes energy and residual in Hartree to the existing log only.
 */
void print_nac_eigenpairs(
    std::ostream &output,
    int ionic_step,
    const std::string &spin_name,
    const double *energies_ry,
    const std::vector<double> &residuals_ry);

/**
 * @brief Reverse a real NAC pair without changing its gap-weighted numerators.
 * @param[in,out] result Owned coupling result; raw, ETF and AO vectors change sign.
 */
void reverse_nac_pair(NacvResult &result);

/**
 * @brief Write one self-describing NAC record to an already-open output stream.
 * @param[in,out] output Caller-owned running log; precision and flags are restored.
 * @param[in] result Couplings and resolution diagnostics for this geometry.
 * @param[in] ionic_step One-based geometry index used in user-facing output.
 * @param[in] bra_state Surface index with zero representing S0.
 * @param[in] ket_state Distinct surface index in the same spin manifold.
 * @param[in] spin_name LR adaptation label, or updown for the unrestricted channel.
 * @details Vectors are bohr^-1. Printed gaps are Hartree and printed numerators
 * are Hartree/bohr; the internal result remains in Rydberg units.
 */
void print_nac_result(
    std::ostream &output,
    const NacvResult &result,
    int ionic_step,
    int bra_state,
    int ket_state,
    const std::string &spin_name);
} // namespace LR

#endif
