#include "spin_constrain.h"

#include "source_base/formatter.h"
#include "source_lcao/module_operator_lcao/dspin_lcao.h"

#include <cmath>
#include <complex>

namespace spinconstrain
{

/**
 * @brief Singleton instance accessor.
 *
 * @details Uses Meyers' Singleton pattern (local static variable).
 * Guaranteed thread-safe initialization in C++11 and later.
 * Each template instantiation (complex<double>, double) gets its own instance.
 */
template <typename TK>
SpinConstrain<TK>& SpinConstrain<TK>::getScInstance()
{
    static SpinConstrain<TK> instance; // Guaranteed to be created and destroyed only once
    return instance;
}

/**
 * @brief Get spin sign for k-point: determines whether this k-point is
 * spin-up (+1) or spin-down (-1) in collinear (nspin=2) calculations.
 *
 * @details In collinear spin, the wavefunction is split into two k-point pools:
 * - isk[ik] == 0: spin-up channel (majority spin) -> sign = +1
 * - isk[ik] == 1: spin-down channel (minority spin) -> sign = -1
 * For non-collinear (npol=2), always returns +1 since both components
 * are handled together.
 *
 * @return +1 for spin-up, -1 for spin-down, +1 for non-collinear
 */
template <typename TK>
int SpinConstrain<TK>::get_spin_sign(int ik) const
{
    if (this->state_.get_npol() == 2) return 1;
    // npol == 1 (nspin == 2): isk[ik]==0 => spin-up (+1), isk[ik]==1 => spin-down (-1)
    return (this->pelec->klist->isk[ik] == 0) ? 1 : -1;
}

template <typename TK>
void SpinConstrain<TK>::set_solver_parameters(const K_Vectors& kv_in,
                                                  void* p_hamilt_in,
                                                  void* psi_in,
                                                  elecstate::ElecState* pelec_in)
{
    this->kv_ = kv_in;
    this->p_hamilt = p_hamilt_in;
    this->psi = psi_in;
    this->pelec = pelec_in;
}

/// @brief  set ParaV
template <typename TK>
void SpinConstrain<TK>::set_ParaV(Parallel_Orbitals* ParaV_in)
{
    this->ParaV = ParaV_in;
    int nloc = this->ParaV->nloc;
    if (nloc <= 0)
    {
        ModuleBase::WARNING_QUIT("SpinConstrain::set_ParaV", "nloc <= 0");
    }
}

/**
 * @brief Reset DeltaSpin operator initialization state.
 *
 * @details The DeltaSpin operator caches internal state (projector matrices, etc.)
 * from a previous SCF iteration. When the constraint parameters change (e.g., new
 * target moments or lambda values), the cached state may be invalid. This function
 * forces the operator to reinitialize on the next call.
 *
 * @par When to call
 * - After changing target_mag_ or constrain_ arrays
 * - When restarting from a previous SCF calculation with different constraints
 * - When switching between LCAO and PW basis sets
 */
template <typename TK>
void SpinConstrain<TK>::reset_dspin_operator()
{
#ifdef __LCAO
    if (this->p_operator == nullptr)
    {
        return;
    }
    if (this->state_.get_nspin() == 4)
    {
        auto* dspin = dynamic_cast<hamilt::DeltaSpin<hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>>*>(this->p_operator);
        if (dspin)
        {
            dspin->reset_initialized();
        }
    }
    else if (this->state_.get_nspin() == 2)
    {
        auto* dspin = dynamic_cast<hamilt::DeltaSpin<hamilt::OperatorLCAO<std::complex<double>, double>>*>(this->p_operator);
        if (dspin)
        {
            dspin->reset_initialized();
        }
    }
#endif
}

template class SpinConstrain<std::complex<double>>;
template class SpinConstrain<double>;

} // namespace spinconstrain
