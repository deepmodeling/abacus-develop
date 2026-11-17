/** @file nacv_tools.h
 * @brief Declare NAC response solvers, result handling and state output utilities.
 */
#ifndef ABACUS_LR_NACV_TOOLS_H
#define ABACUS_LR_NACV_TOOLS_H

#include "source_base/matrix.h"
#include <algorithm>
#include <cctype>
#include <complex>
#include <functional>
#include <iosfwd>
#include <string>
#include <vector>

namespace LR
{

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
 * @brief Check deferred surface bounds after the LR root count is resolved.
 * @param[in] pairs INPUT-validated nonnegative surface indices; zero denotes S0.
 * @param[in] state_count Actual number of solved excited states.
 * @throw std::invalid_argument If a requested root exceeds the solved window.
 * @details INPUT may leave lr_nstates nonpositive to request an automatic count.
 */
void validate_nac_surfaces(const std::vector<int> &pairs, int state_count);

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

namespace LR
{
/**
 * @brief Solve a distributed real orbital-response equation and certify its residual.
 * @param[out] Z Local state-major solution, with nstates blocks of ld elements.
 * @param[in] R Local right-hand side, in Rydberg for the physical A+B response operator.
 * @param[in] ld Local vector length, including zero on ranks without entries.
 * @param[in] nstates Number of right-hand sides solved with a combined residual norm.
 * @param[in] f_LZ Collective linear operator; writes all local output entries.
 * @param[in] test_force Whether to print the full final vector for diagnostics.
 * @param[in] tolerance Positive finite absolute global Euclidean norm bound on LZ-R.
 * The physical response residual and tolerance have units Ry; Z is dimensionless.
 * @param[in] max_iterations Positive total update limit, including residual corrections.
 * @pre Dimensions are the prepared NAC response layout; controls passed INPUT validation.
 * @throw std::runtime_error If the explicitly recomputed residual is not below tolerance.
 * @details A converged recurrence is checked against LZ-R. When these differ,
 * retain Z and restart a correction equation with the true residual, using the
 * same tolerance and remaining global update budget. Apply a correction update
 * before another recursive-norm check, so rounding differences cannot cause
 * an unlimited sequence of zero-update restarts. R is never modified.
 */
void solve_nac_cg(
    double *const Z,
    double *R,
    const int &ld,
    const int &nstates,
    std::function<void(const double *const, double *const)> f_LZ,
    const bool test_force,
    const double tolerance,
    const int max_iterations);

/**
 * @brief Reject the unavailable complex CG response solver explicitly.
 * @param[out] Z Reserved local solution storage; not modified.
 * @param[in] R Reserved local right-hand side; not modified.
 * @param[in] ld Local vector length.
 * @param[in] nstates Number of right-hand sides.
 * @param[in] f_LZ Reserved complex operator.
 * @param[in] test_force Reserved diagnostic switch.
 * @param[in] tolerance Reserved absolute residual tolerance.
 * @param[in] max_iterations Reserved iteration bound.
 * @throw std::runtime_error Always: complex CG response is not implemented.
 */
void solve_nac_cg(
    std::complex<double> *const Z,
    std::complex<double> *R,
    const int &ld,
    const int &nstates,
    std::function<void(const std::complex<double> *const, std::complex<double> *const)> f_LZ,
    const bool test_force,
    const double tolerance,
    const int max_iterations);

} // namespace LR

class Parallel_Blkcyclic;
namespace LR
{
/**
 * @brief Own one real NAC response factorization for repeated right-hand sides.
 * @details Matrices are column-major. Distributed layouts are borrowed and
 * must outlive the factorization. No matrix inverse or tolerance change is used.
 */
class NacvFactorization
{
  public:
    /**
     * @param[in] matrix Complete orbital Hessian in Ry; ownership is transferred.
     * @param[in] dimension Global matrix dimension.
     */
    NacvFactorization(std::vector<double> matrix, int dimension);
#ifdef __MPI
    /**
     * @param[in] matrix Local Hessian block in Ry; ownership is transferred.
     * @param[in] layout Square block-cyclic layout, with the existing BLACS context.
     * @param[in] cholesky Use the existing symmetrized positive-definite algorithm.
     */
    NacvFactorization(std::vector<double> matrix, const Parallel_Blkcyclic &layout, bool cholesky);
    /**
     * @param[in,out] rhs Local RHS blocks in Ry, overwritten with dimensionless Z.
     * @param[in] layout RHS layout with the same context and block size as the Hessian.
     */
    void solve(std::vector<double> &rhs, const Parallel_Blkcyclic &layout);
#endif
    /**
     * @param[in,out] rhs Complete RHS columns in Ry, overwritten with dimensionless Z.
     * @param[in] count Number of RHS columns.
     */
    void solve(std::vector<double> &rhs, int count);

  private:
    int dimension_ = 0;
    const Parallel_Blkcyclic *layout_ = nullptr;
    bool cholesky_ = false;
    std::vector<double> factors_;
    std::vector<int> pivots_;
};
} // namespace LR

namespace LR
{
/**
 * @brief Write existing LR eigenstate files with full double precision.
 * @param[in] directory Existing output directory including its trailing separator.
 * @param[in] label Spin label used by existing energy/amplitude filenames.
 * @param[in] suffix Existing optional ionic-step filename suffix.
 * @param[in] rank Zero-based rank; energy output is restricted to rank zero.
 * @param[in] energies Excitation energies in Ry, one per state.
 * @param[in] amplitudes State-major local amplitude array.
 * @param[in] local_size Entries per state on this rank.
 * @param[in] state_count Number of states to write.
 * @param[in,out] output Borrowed running log for successful absolute-path confirmations.
 * @throw std::runtime_error On all ranks if any state output file fails.
 * @pre All ranks in the LR world domain call this collectively; rank is the world rank.
 * @details Values are unfiltered, with 17 digits after the decimal in scientific notation.
 * Output calls execute with NDEBUG as well as in assertion-enabled builds.
 */
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
    std::ostream &output);
} // namespace LR

namespace LR
{
namespace NacvXC
{
/**
 * @brief Resolve built-in ground-state aliases before selecting a response kernel.
 * @param[in] name Case-insensitive functional name from INPUT or a pseudopotential.
 * @return Lowercase canonical name; unknown names are only converted to lowercase.
 * @details This maps names, not physical approximations: NONE does not become RPA.
 */
inline std::string canonical_xc_name(const std::string &name)
{
    std::string canonical = name;
    std::transform(
        canonical.begin(),
        canonical.end(),
        canonical.begin(),
        [](unsigned char character)
        {
            return static_cast<char>(std::tolower(character));
        });
    if (canonical == "pz" || canonical == "slapznogxnogc")
    {
        return "lda";
    }
    if (canonical == "slapwpbxpbc")
    {
        return "pbe";
    }
    if (canonical == "nox+noc")
    {
        return "none";
    }
    return canonical;
}

} // namespace NacvXC
} // namespace LR
namespace LR
{
/** @brief Time a complete NACV operation, including exceptional exits. */
class NacvTimerScope
{
  public:
    /**
     * @param[in] category Owning class or module label, valid for this scope.
     * @param[in] name Enclosing function name, valid for this scope.
     */
    NacvTimerScope(const char *category, const char *name);
    ~NacvTimerScope();
    NacvTimerScope(const NacvTimerScope &) = delete;
    NacvTimerScope &operator=(const NacvTimerScope &) = delete;

  private:
    const char *const category_;
    const char *const name_;
};
} // namespace LR

#endif
