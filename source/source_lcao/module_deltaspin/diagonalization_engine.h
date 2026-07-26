/**
 * @file diagonalization_engine.h
 * @brief Abstract interface and implementations for DeltaSpin diagonalization strategies.
 *
 * @par Purpose
 * Replaces the hard-coded if-else branches in cal_mw_from_lambda() with a
 * strategy-pattern-based architecture. Three diagonalization strategies are supported:
 *
 *   1. FullSpaceDiagonalizer: Full HSolverLCAO diagonalization (ground truth)
 *   2. SubspaceDiagonalizer: Subspace H_sub = H0_sub + dl*P_I_sub, then diag_hegvd
 *   3. FirstOrderResponseEngine: Eigenvalue shifts only (eps_new = eps_old + dl*P_diag)
 *
 * @par Design principles
 * - Strategy pattern: interchangeable diagonalization algorithms via common interface
 * - RAII memory management: SubspaceCache uses std::vector instead of raw new[]/delete[]
 * - State isolation: each engine manages its own temporary variables
 * - Type-safe: template-specialized for std::complex<double>
 *
 * @par Integration
 * SpinConstrain holds a unique_ptr<DiagonalizationEngine> and calls solve()
 * instead of the previous if-else branching logic.
 */

#ifndef DIAGONALIZATION_ENGINE_H
#define DIAGONALIZATION_ENGINE_H

#include <complex>
#include <memory>
#include <string>
#include <vector>

#include "source_base/vector3.h"

namespace spinconstrain
{

// Forward declaration
template <typename TK>
class SpinConstrain;

/**
 * @brief Result of a diagonalization operation.
 */
struct DiagonalizationResult
{
    bool success = false;                ///< Whether diagonalization succeeded
    bool used_subspace_approximation = false;  ///< True if subspace/first-order was used
    double max_eigenvalue_change = 0.0;  ///< Maximum |eps_new - eps_old| across all bands
    int nbands = 0;                      ///< Number of bands
    int nk = 0;                          ///< Number of k-points
};

/**
 * @brief Enum for diagonalization strategy type.
 */
enum class DiagonalizationStrategy
{
    FullSpace,       ///< Full HSolverLCAO diagonalization
    Subspace,        ///< Subspace diagonalization with wavefunction rotation
    FirstOrder       ///< First-order eigenvalue response (no wavefunction change)
};

/**
 * @brief Abstract base class for diagonalization engines in DeltaSpin.
 *
 * @par Lifecycle
 *   1. ESolver creates the appropriate engine based on input parameters
 *   2. For subspace strategies: build_subspace() is called once near convergence
 *   3. Each lambda step: solve() is called to compute eigenvalues and Mi
 *   4. At SCF iteration end: clear_subspace() is called
 *
 * @par Thread safety
 * Not thread-safe. Each SpinConstrain instance owns its own engine.
 */
class DiagonalizationEngine
{
public:
    virtual ~DiagonalizationEngine() = default;

    /**
     * @brief Solve the Hamiltonian for the current lambda values.
     *
     * @details This is the core method called repeatedly during lambda optimization.
     * The engine updates pelec->ekb, pelec->wg, and computes Mi via SpinConstrain.
     *
     * @param i_step Inner lambda step index (-1=init, 0+=optimization, -2=build subspace)
     * @return DiagonalizationResult with success status and metadata
     */
    virtual DiagonalizationResult solve(int i_step) = 0;

    /**
     * @brief Build the subspace cache at the current lambda reference point.
     *
     * @details Performs a full diagonalization, then extracts H0_sub, S_sub, P_I_sub.
     * This is called once when acceleration is activated (RMS drops below threshold).
     *
     * @param lambda_ref Lambda values at the reference point
     * @return true if subspace was successfully built
     */
    virtual bool build_subspace(const std::vector<ModuleBase::Vector3<double>>& lambda_ref) = 0;

    /**
     * @brief Check if subspace cache is valid and ready for accelerated solves.
     */
    virtual bool has_subspace() const = 0;

    /**
     * @brief Clear the subspace cache and release associated memory.
     *
     * @details Called at the end of a lambda loop or when switching strategies.
     */
    virtual void clear_subspace() = 0;

    /**
     * @brief Get the human-readable name of this strategy.
     */
    virtual std::string name() const = 0;

    /**
     * @brief Get the strategy type enum.
     */
    virtual DiagonalizationStrategy type() const = 0;
};

/**
 * @brief RAII container for LCAO subspace cache data.
 *
 * @par Memory layout
 * All matrices are stored in column-major order (LAPACK convention).
 * - H0_sub: [nk * nbands * nbands], flattened 1D array
 * - S_sub:  [nk * nbands * nbands], flattened 1D array
 * - P_I_sub: nested vector [nk][nat][nbands*nbands] for per-atom projector matrices
 *
 * @par RAII guarantee
 * All memory is managed via std::vector. No manual new[]/delete[] is needed.
 * The destructor automatically releases all memory when the cache goes out of scope.
 */
class SubspaceCache
{
public:
    SubspaceCache();
    ~SubspaceCache() = default;

    // Non-copyable, movable
    SubspaceCache(const SubspaceCache&) = delete;
    SubspaceCache& operator=(const SubspaceCache&) = delete;
    SubspaceCache(SubspaceCache&&) = default;
    SubspaceCache& operator=(SubspaceCache&&) = default;

    /**
     * @brief Build the subspace cache from computed data.
     *
     * @param nk Number of k-points
     * @param nbands Number of bands
     * @param nat Number of atoms
     * @param H0_sub_raw H0_sub for all k-points [nk * nbands^2] (column-major)
     * @param S_sub_raw S_sub for all k-points [nk * nbands^2] (column-major)
     * @param P_I_sub_all Per-k, per-atom projector matrices [nk][iat][nbands^2]
     * @param ekb_ref_all Eigenvalues from full diagonalization [nk * nbands]
     * @param lambda_ref Lambda values at the reference point [nat]
     */
    void build(int nk,
               int nbands,
               int nat,
               const std::complex<double>* H0_sub_raw,
               const std::complex<double>* S_sub_raw,
               std::vector<std::vector<std::vector<std::complex<double>>>> P_I_sub_all,
               const std::vector<double>& ekb_ref_all,
               const std::vector<ModuleBase::Vector3<double>>& lambda_ref);

    /// @brief Clear all cached data and reset to empty state.
    void clear();

    /// @brief Check if the cache has been built and is valid.
    bool is_valid() const { return valid_; }

    /// @brief Number of k-points in the cache.
    int nk() const { return nk_; }

    /// @brief Number of bands in the cache.
    int nbands() const { return nbands_; }

    /// @brief Number of atoms in the cache.
    int nat() const { return nat_; }

    // Accessors (const references to prevent accidental modification)

    /// @brief H0_sub for k-point ik [nbands * nbands, column-major]
    const std::complex<double>* H0_sub(int ik) const;

    /// @brief S_sub for k-point ik [nbands * nbands, column-major]
    const std::complex<double>* S_sub(int ik) const;

    /// @brief P_I_sub for k-point ik, atom iat [nbands * nbands, column-major]
    /// Returns nullptr if the atom is unconstrained (empty matrix).
    const std::complex<double>* P_I_sub(int ik, int iat) const;

    /// @brief Reference eigenvalues for k-point ik [nbands]
    const double* ekb_ref(int ik) const;

    /// @brief Reference lambda values at which the subspace was built [nat]
    const std::vector<ModuleBase::Vector3<double>>& lambda_ref() const;

    /// @brief Full P_I_sub data [nk][nat][nbands^2] for cal_mi_lcao_subspace
    const std::vector<std::vector<std::vector<std::complex<double>>>>& P_I_sub_all() const
    {
        return P_I_sub_;
    }

private:
    int nk_ = 0;
    int nbands_ = 0;
    int nat_ = 0;
    bool valid_ = false;

    // Flat arrays for H0_sub and S_sub (column-major, nk * nbands^2)
    std::vector<std::complex<double>> H0_sub_;
    std::vector<std::complex<double>> S_sub_;

    // P_I_sub: [nk][nat][nbands^2], empty vectors for unconstrained atoms
    std::vector<std::vector<std::vector<std::complex<double>>>> P_I_sub_;

    // Reference eigenvalues [nk * nbands]
    std::vector<double> ekb_ref_;

    // Reference lambda values [nat]
    std::vector<ModuleBase::Vector3<double>> lambda_ref_;
};

/**
 * @brief Full-space diagonalizer using HSolverLCAO.
 *
 * @par Algorithm
 *   1. Update DeltaSpin operator with current lambda
 *   2. Call HSolverLCAO::solve() (full O(N^3) diagonalization)
 *   3. Compute weights from eigenvalues
 *   4. Compute Mi via cal_mi_lcao()
 *
 * @par Characteristics
 * - Always correct (ground truth)
 * - Most expensive: O(N^3) per step
 * - No subspace cache needed
 * - Works for both nspin=2 and nspin=4
 */
class FullSpaceDiagonalizer : public DiagonalizationEngine
{
public:
    explicit FullSpaceDiagonalizer(SpinConstrain<std::complex<double>>& sc);

    DiagonalizationResult solve(int i_step) override;
    bool build_subspace(const std::vector<ModuleBase::Vector3<double>>& lambda_ref) override;
    bool has_subspace() const override;
    void clear_subspace() override;
    std::string name() const override { return "FullSpace"; }
    DiagonalizationStrategy type() const override { return DiagonalizationStrategy::FullSpace; }

private:
    SpinConstrain<std::complex<double>>& sc_;
};

/**
 * @brief Subspace diagonalizer with wavefunction rotation.
 *
 * @par Algorithm
 *   Phase 1 (build_subspace):
 *     1. Full diagonalization to get correct psi at lambda_ref
 *     2. Compute H0_sub = C^dag H C, S_sub = C^dag S C
 *     3. Compute P_I_sub = C^dag D_I C for each constrained atom
 *
 *   Phase 2 (solve):
 *     1. H_sub(lambda) = H0_sub + sum_I (lambda_I - lambda_ref_I) * P_I_sub
 *     2. diag_hegvd: H_sub V = S_sub V eps
 *     3. Rotate psi: C_new = C_original * V
 *     4. Build full-space DM from rotated psi -> compute Mi
 *     5. Restore original psi (BFGS loop expects unmodified state)
 *
 * @par Characteristics
 * - Valid only for small delta-lambda (near convergence)
 * - Much faster than full diag: O(nbands^3) << O(NLOCAL^3)
 * - Requires subspace cache (build_subspace must be called first)
 * - Currently only supports nspin=2
 */
class SubspaceDiagonalizer : public DiagonalizationEngine
{
public:
    explicit SubspaceDiagonalizer(SpinConstrain<std::complex<double>>& sc);

    DiagonalizationResult solve(int i_step) override;
    bool build_subspace(const std::vector<ModuleBase::Vector3<double>>& lambda_ref) override;
    bool has_subspace() const override;
    void clear_subspace() override;
    std::string name() const override { return "Subspace"; }
    DiagonalizationStrategy type() const override { return DiagonalizationStrategy::Subspace; }

private:
    SpinConstrain<std::complex<double>>& sc_;
    SubspaceCache cache_;

    // Internal temporary buffers (reused across solve calls)
    std::vector<std::complex<double>> h_tmp_;
    std::vector<std::complex<double>> s_tmp_;
    std::vector<std::complex<double>> s_copy_;
    std::vector<double> eigenvalues_;
    std::vector<std::complex<double>> eigenvectors_;
};

/**
 * @brief First-order eigenvalue response engine.
 *
 * @par Algorithm
 *   Phase 1 (build_subspace): same as SubspaceDiagonalizer
 *
 *   Phase 2 (solve):
 *     1. delta_eps_ib = sum_I (lambda_I - lambda_ref_I) * diag(P_I_sub)_ib
 *     2. eps_new = eps_ref + delta_eps
 *     3. Wavefunctions unchanged (V = identity)
 *     4. Compute Mi from original psi with updated Fermi weights
 *
 * @par Characteristics
 * - Fastest approximation: O(nbands * nat) per step
 * - Valid only for very small delta-lambda
 * - Wavefunction mixing is neglected
 * - Currently only supports nspin=2
 */
class FirstOrderResponseEngine : public DiagonalizationEngine
{
public:
    explicit FirstOrderResponseEngine(SpinConstrain<std::complex<double>>& sc);

    DiagonalizationResult solve(int i_step) override;
    bool build_subspace(const std::vector<ModuleBase::Vector3<double>>& lambda_ref) override;
    bool has_subspace() const override;
    void clear_subspace() override;
    std::string name() const override { return "FirstOrder"; }
    DiagonalizationStrategy type() const override { return DiagonalizationStrategy::FirstOrder; }

private:
    SpinConstrain<std::complex<double>>& sc_;
    SubspaceCache cache_;
};

/**
 * @brief Factory function to create the appropriate diagonalization engine.
 *
 * @param strategy The desired diagonalization strategy
 * @param sc Reference to the SpinConstrain instance
 * @return unique_ptr to the created engine
 */
std::unique_ptr<DiagonalizationEngine> create_diagonalization_engine(
    DiagonalizationStrategy strategy,
    SpinConstrain<std::complex<double>>& sc);

/**
 * @brief Convert string to DiagonalizationStrategy (case-insensitive).
 *
 * @param s Strategy name: "full", "fullspace", "subspace", "subspace_diag",
 *          "first_order", "firstorder", "linear_response"
 * @return Corresponding enum value
 */
DiagonalizationStrategy strategy_from_string(const std::string& s);

/**
 * @brief Convert DiagonalizationStrategy to string.
 */
std::string strategy_to_string(DiagonalizationStrategy s);

} // namespace spinconstrain

#endif // DIAGONALIZATION_ENGINE_H
