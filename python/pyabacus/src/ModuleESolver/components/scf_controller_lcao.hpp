/**
 * @file scf_controller_lcao.hpp
 * @brief LCAO SCF controller implementation
 *
 * Implements the ISCFController interface for LCAO calculations.
 */

#ifndef PYABACUS_ESOLVER_SCF_CONTROLLER_LCAO_HPP
#define PYABACUS_ESOLVER_SCF_CONTROLLER_LCAO_HPP

#include "../interfaces/i_scf_controller.hpp"
#include "../interfaces/i_hamiltonian_builder.hpp"
#include "../interfaces/i_charge_mixer.hpp"
#include "../interfaces/i_diagonalizer.hpp"
#include "hamiltonian_builder_lcao.hpp"
#include "charge_mixer_wrapper.hpp"
#include "diagonalizer_wrapper.hpp"

#include <memory>
#include <iostream>

namespace pyabacus {
namespace esolver {

/**
 * @brief LCAO SCF controller implementation
 *
 * This class implements the SCF loop for LCAO calculations,
 * coordinating the Hamiltonian builder, charge mixer, and diagonalizer.
 *
 * @tparam TK Type for k-space quantities
 * @tparam TR Type for real-space quantities
 */
template <typename TK, typename TR = double>
class SCFControllerLCAO : public ISCFController
{
public:
    SCFControllerLCAO() = default;

    SCFControllerLCAO(int nbasis, int nks, int nbands, int nspin, int nrxx)
        : nbasis_(nbasis)
        , nks_(nks)
        , nbands_(nbands)
        , nspin_(nspin)
        , nrxx_(nrxx)
    {
        // Initialize components
        hamilt_builder_ = std::make_unique<HamiltonianBuilderLCAO<TK, TR>>(
            nbasis, nks, nbasis, nbasis);
        charge_mixer_ = std::make_unique<ChargeMixerWrapper>(nspin, nrxx);
        diagonalizer_ = std::make_unique<DiagonalizerWrapper<TK>>(nbasis, nbands);
    }

    ~SCFControllerLCAO() override = default;

    // ==================== Lifecycle ====================

    void initialize(int istep) override
    {
        istep_ = istep;
        iteration_ = 0;
        drho_ = 0.0;
        energy_ = 0.0;
        status_ = SCFStatus::NotStarted;

        // Reset charge mixer
        if (charge_mixer_)
        {
            charge_mixer_->reset();
        }

        initialized_ = true;
    }

    void finalize(int istep) override
    {
        // Cleanup after SCF
        initialized_ = false;
    }

    // ==================== Iteration Control ====================

    SCFStatus run_iteration(int iter) override
    {
        if (!initialized_)
        {
            throw std::runtime_error("SCF not initialized. Call initialize() first.");
        }

        iteration_ = iter;
        status_ = SCFStatus::Running;

        // SCF iteration steps:
        // 1. Build Hamiltonian from current charge density
        // 2. Diagonalize to get new wave functions
        // 3. Calculate new charge density from wave functions
        // 4. Mix old and new charge densities
        // 5. Check convergence

        // For now, this is a placeholder that simulates convergence
        // In full implementation, would call actual ABACUS routines

        // Simulate drho decreasing
        drho_ = 1.0 / (iter * iter + 1);

        // Simulate energy
        energy_ = -10.0 - 0.1 / iter;

        return status_;
    }

    SCFStatus run_scf(const SCFConvergenceCriteria& criteria,
                      SCFIterationCallback callback) override
    {
        initialize(istep_);
        status_ = SCFStatus::Running;

        for (int iter = 1; iter <= criteria.max_iterations; ++iter)
        {
            run_iteration(iter);

            // Check convergence
            bool converged = true;
            if (criteria.check_drho && drho_ > criteria.drho_threshold)
            {
                converged = false;
            }

            // Call callback if provided
            if (callback)
            {
                bool continue_scf = callback(iter, drho_, energy_);
                if (!continue_scf)
                {
                    status_ = SCFStatus::Failed;
                    return status_;
                }
            }

            if (converged)
            {
                status_ = SCFStatus::Converged;
                return status_;
            }
        }

        status_ = SCFStatus::MaxIterReached;
        return status_;
    }

    bool is_converged() const override
    {
        return status_ == SCFStatus::Converged;
    }

    SCFStatus get_status() const override { return status_; }

    // ==================== State Queries ====================

    int get_iteration() const override { return iteration_; }

    double get_drho() const override { return drho_; }

    double get_energy() const override { return energy_; }

    // ==================== Component Access ====================

    void* get_hamiltonian_builder() override
    {
        return hamilt_builder_.get();
    }

    IChargeMixer* get_charge_mixer() override
    {
        return charge_mixer_.get();
    }

    void* get_diagonalizer() override
    {
        return diagonalizer_.get();
    }

    // ==================== Typed Component Access ====================

    HamiltonianBuilderLCAO<TK, TR>* get_hamiltonian_builder_typed()
    {
        return hamilt_builder_.get();
    }

    DiagonalizerWrapper<TK>* get_diagonalizer_typed()
    {
        return diagonalizer_.get();
    }

    // ==================== Configuration ====================

    void set_convergence_criteria(const SCFConvergenceCriteria& criteria)
    {
        criteria_ = criteria;
    }

    SCFConvergenceCriteria get_convergence_criteria() const
    {
        return criteria_;
    }

    void set_mixing_config(const MixingConfig& config)
    {
        if (charge_mixer_)
        {
            charge_mixer_->set_config(config);
        }
    }

    void set_diag_config(const DiagConfig& config)
    {
        if (diagonalizer_)
        {
            diagonalizer_->set_config(config);
        }
    }

private:
    // Dimensions
    int nbasis_ = 0;
    int nks_ = 0;
    int nbands_ = 0;
    int nspin_ = 1;
    int nrxx_ = 0;

    // State
    int istep_ = 0;
    int iteration_ = 0;
    double drho_ = 0.0;
    double energy_ = 0.0;
    bool initialized_ = false;
    SCFStatus status_ = SCFStatus::NotStarted;

    // Configuration
    SCFConvergenceCriteria criteria_;

    // Components
    std::unique_ptr<HamiltonianBuilderLCAO<TK, TR>> hamilt_builder_;
    std::unique_ptr<ChargeMixerWrapper> charge_mixer_;
    std::unique_ptr<DiagonalizerWrapper<TK>> diagonalizer_;
};

// Type aliases
using SCFControllerLCAOGamma = SCFControllerLCAO<double, double>;
using SCFControllerLCAOMultiK = SCFControllerLCAO<std::complex<double>, double>;

} // namespace esolver
} // namespace pyabacus

#endif // PYABACUS_ESOLVER_SCF_CONTROLLER_LCAO_HPP
