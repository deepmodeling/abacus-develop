/**
 * @file i_scf_controller.hpp
 * @brief Abstract interface for SCF controller
 *
 * Defines the interface for controlling SCF iterations,
 * allowing different implementations and Python-controlled workflows.
 */

#ifndef PYABACUS_ESOLVER_I_SCF_CONTROLLER_HPP
#define PYABACUS_ESOLVER_I_SCF_CONTROLLER_HPP

#include <functional>
#include <string>

namespace pyabacus {
namespace esolver {

/**
 * @brief SCF convergence status
 */
enum class SCFStatus
{
    NotStarted,     ///< SCF has not been started
    Running,        ///< SCF is in progress
    Converged,      ///< SCF converged successfully
    MaxIterReached, ///< Maximum iterations reached without convergence
    Failed          ///< SCF failed due to error
};

/**
 * @brief SCF convergence criteria
 */
struct SCFConvergenceCriteria
{
    double drho_threshold = 1e-6;   ///< Charge density difference threshold
    double energy_threshold = 1e-6; ///< Energy difference threshold (Ry)
    int max_iterations = 100;       ///< Maximum number of iterations
    bool check_energy = true;       ///< Whether to check energy convergence
    bool check_drho = true;         ///< Whether to check charge density convergence
};

/**
 * @brief Callback type for SCF iteration events
 *
 * Called after each SCF iteration with:
 * - iter: current iteration number
 * - drho: charge density difference
 * - energy: current total energy
 *
 * Return true to continue, false to stop SCF.
 */
using SCFIterationCallback = std::function<bool(int iter, double drho, double energy)>;

// Forward declarations for component interfaces
template <typename TK, typename TR> class IHamiltonianBuilder;
class IChargeMixer;
template <typename TK> class IDiagonalizer;

/**
 * @brief Abstract interface for SCF controller
 *
 * This interface defines the contract for SCF loop control,
 * allowing different implementations (LCAO, PW, etc.) and
 * enabling Python-controlled SCF workflows with breakpoints.
 */
class ISCFController
{
public:
    virtual ~ISCFController() = default;

    // ==================== Lifecycle ====================

    /**
     * @brief Initialize SCF calculation
     * @param istep Ion step index (for MD/relaxation)
     */
    virtual void initialize(int istep = 0) = 0;

    /**
     * @brief Finalize SCF calculation
     * @param istep Ion step index
     */
    virtual void finalize(int istep = 0) = 0;

    // ==================== Iteration Control ====================

    /**
     * @brief Run a single SCF iteration
     * @param iter Iteration number (1-based)
     * @return Current SCF status
     */
    virtual SCFStatus run_iteration(int iter) = 0;

    /**
     * @brief Run complete SCF loop
     * @param criteria Convergence criteria
     * @param callback Optional callback for each iteration
     * @return Final SCF status
     */
    virtual SCFStatus run_scf(const SCFConvergenceCriteria& criteria,
                              SCFIterationCallback callback = nullptr) = 0;

    /**
     * @brief Check if SCF is converged
     * @return true if converged
     */
    virtual bool is_converged() const = 0;

    /**
     * @brief Get current SCF status
     * @return Current status
     */
    virtual SCFStatus get_status() const = 0;

    // ==================== State Queries ====================

    /**
     * @brief Get current iteration number
     * @return Iteration number (0 if not started)
     */
    virtual int get_iteration() const = 0;

    /**
     * @brief Get current charge density difference
     * @return drho value
     */
    virtual double get_drho() const = 0;

    /**
     * @brief Get current total energy
     * @return Total energy in Ry
     */
    virtual double get_energy() const = 0;

    // ==================== Component Access ====================

    /**
     * @brief Get Hamiltonian builder component
     * @return Pointer to Hamiltonian builder (may be nullptr)
     */
    virtual void* get_hamiltonian_builder() = 0;

    /**
     * @brief Get charge mixer component
     * @return Pointer to charge mixer (may be nullptr)
     */
    virtual IChargeMixer* get_charge_mixer() = 0;

    /**
     * @brief Get diagonalizer component
     * @return Pointer to diagonalizer (may be nullptr)
     */
    virtual void* get_diagonalizer() = 0;
};

} // namespace esolver
} // namespace pyabacus

#endif // PYABACUS_ESOLVER_I_SCF_CONTROLLER_HPP
