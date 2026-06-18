#include <iostream>

#include "source_base/matrix.h"
#include "source_base/name_angular.h"
#include "source_base/parallel_reduce.h"
#include "source_base/tool_title.h"
#include "source_base/timer.h"
#include "source_pw/module_pwdft/onsite_proj.h"
#include "spin_constrain.h"
#include "source_io/module_parameter/parameter.h"
#ifdef __LCAO
#include "source_estate/elecstate_lcao.h"
#include "source_lcao/hamilt_lcao.h"
#include "source_lcao/module_operator_lcao/dspin_lcao.h"

/**
 * @file cal_mw.cpp
 * @brief Magnetic moment calculation for LCAO and PW basis sets.
 *
 * @par cal_mi_lcao (LCAO)
 * Uses the DeltaSpin operator to compute magnetic moments from the density
 * matrix via real-space projection. For nspin=2, only the z-component is
 * extracted. For nspin=4, all three components are extracted from the
 * interleaved 4-component spinor density matrix.
 *
 * @par cal_mi_pw (PW)
 * Uses the OnsiteProjector to compute atomic projections <alpha_{l,m}|psi_{k,i}>
 * (becp coefficients), then decomposes these into magnetic moments using
 * Pauli matrix traces (accumulate_Mi_from_becp).
 *
 * @par Error conditions
 * - Dynamic cast failure: p_operator is not the correct DeltaSpin type.
 *   This happens if set_operator() was not called with the correct type.
 *   Solution: Ensure set_operator() is called before cal_mi_lcao().
 */

/**
 * @brief Calculate atomic magnetic moments using real-space projection (LCAO basis).
 *
 * @details The DeltaSpin operator computes magnetic moments by projecting the
 * density matrix onto atomic orbitals. For each constrained atom:
 *   M_i = Tr[P_at * (rho_up - rho_dn)]  (nspin=2)
 *   M_i = Tr[P_at * rho_spinor]          (nspin=4, decomposed via Pauli matrices)
 *
 * @param step Current SCF iteration number (for logging)
 * @param print Whether to print moments (unused in this implementation)
 */
template <>
void spinconstrain::SpinConstrain<std::complex<double>>::cal_mi_lcao(const int& step, bool print)
{
    ModuleBase::TITLE("module_deltaspin", "cal_mi_lcao");
    ModuleBase::timer::start("spinconstrain::SpinConstrain", "cal_mi_lcao");
    // Reset Mi before calculation
    this->zero_Mi();
    const hamilt::HContainer<double>* dmr = this->dm_->get_DMR_pointer(1);
    std::vector<double> moments;
    if(this->nspin_==2)
    {
        // Switch to spin-difference density matrix (rho_up - rho_dn)
        this->dm_->switch_dmr(2);

        // Compute moments via DeltaSpin operator
        moments = static_cast<hamilt::DeltaSpin<hamilt::OperatorLCAO<std::complex<double>, double>>*>(this->p_operator)->cal_moment(dmr, this->get_constrain());

        // Switch back to total density matrix
        this->dm_->switch_dmr(0);

        // For nspin=2, only z-component is meaningful
        for(int iat=0;iat<this->Mi_.size();iat++)
        {
            this->Mi_[iat].x = 0.0;
            this->Mi_[iat].y = 0.0;
            this->Mi_[iat].z = moments[iat];
        }
    }
    else if(this->nspin_==4)
    {
        // For nspin=4, moments array contains interleaved [Mx, My, Mz] per atom
        moments = static_cast<hamilt::DeltaSpin<hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>>*>(this->p_operator)->cal_moment(dmr, this->get_constrain());
        for(int iat=0;iat<this->Mi_.size();iat++)
        {
            this->Mi_[iat].x = moments[iat*3];
            this->Mi_[iat].y = moments[iat*3+1];
            this->Mi_[iat].z = moments[iat*3+2];
        }
    }

    ModuleBase::timer::end("spinconstrain::SpinConstrain", "cal_mi_lcao");
}

#endif

/**
 * @brief Calculate atomic magnetic moments using projector overlap (PW basis).
 *
 * @details For each k-point:
 *   1. Tabulate atomic projectors: set up |alpha_{l,m}> for each atom
 *   2. Compute becp = <alpha_{l,m}|psi_{k,i}> via overlap_proj_psi
 *   3. Decompose becp into magnetic moments via accumulate_Mi_from_becp
 *
 * The magnetic moment is computed as:
 *   Mi = sum_{k,i} w_{k,i} * <psi_{k,i}|P_at|sigma|psi_{k,i}>
 * where P_at is the atomic projector and sigma are the Pauli matrices.
 *
 * Finally, Mi is summed across all MPI k-pool ranks since each pool only
 * has a subset of k-points.
 */
template <>
void spinconstrain::SpinConstrain<std::complex<double>>::cal_mi_pw()
{
    ModuleBase::TITLE("module_deltaspin", "cal_mi_pw");
    ModuleBase::timer::start("spinconstrain::SpinConstrain", "cal_mi_pw");

    this->zero_Mi();
    if(PARAM.inp.device == "cpu")
    {
        auto* onsite_p = projectors::OnsiteProjector<double, base_device::DEVICE_CPU>::get_instance();
        // Loop over k-points to calculate Mi of sum_{k,i,l,m}<Psi_{k,i}|alpha_{l,m}><alpha_{l,m}|Psi_{k,i}>
        std::complex<double>* psi_pointer = nullptr;
        psi::Psi<std::complex<double>, base_device::DEVICE_CPU>* psi_t = static_cast<psi::Psi<std::complex<double>, base_device::DEVICE_CPU>*>(this->psi);
        const int nbands = psi_t->get_nbands();
        const int nks = psi_t->get_nk();
        const int npol = psi_t->get_npol();
        for(int ik = 0; ik < nks; ik++)
        {
            psi_t->fix_k(ik);
            psi_pointer = psi_t->get_pointer();
            onsite_p->tabulate_atomic(ik); // Set up atomic projectors for this k-point
            onsite_p->overlap_proj_psi(nbands * npol, psi_pointer); // Compute becp = <alpha|psi>
            const std::complex<double>* becp = onsite_p->get_h_becp();
            int nkb = onsite_p->get_tot_nproj();
            this->accumulate_Mi_from_becp(becp, nkb, nbands, npol, ik,
                &this->pelec->wg(ik, 0), &onsite_p->get_nh(0));
        }
    }
#if ((defined __CUDA) || (defined __ROCM))
    else
    {
        auto* onsite_p = projectors::OnsiteProjector<double, base_device::DEVICE_GPU>::get_instance();
        std::complex<double>* psi_pointer = nullptr;
        psi::Psi<std::complex<double>, base_device::DEVICE_GPU>* psi_t = static_cast<psi::Psi<std::complex<double>, base_device::DEVICE_GPU>*>(this->psi);
        const int nbands = psi_t->get_nbands();
        const int nks = psi_t->get_nk();
        const int npol = psi_t->get_npol();
        for(int ik = 0; ik < nks; ik++)
        {
            psi_t->fix_k(ik);
            psi_pointer = psi_t->get_pointer();
            onsite_p->tabulate_atomic(ik);
            onsite_p->overlap_proj_psi(nbands * npol, psi_pointer);
            const std::complex<double>* becp = onsite_p->get_h_becp();
            int nkb = onsite_p->get_size_becp() / nbands / npol;
            this->accumulate_Mi_from_becp(becp, nkb, nbands, npol, ik,
                &this->pelec->wg(ik, 0), &onsite_p->get_nh(0));
        }
    }
#endif
    // MPI reduction: sum Mi across all k-pool ranks
    Parallel_Reduce::reduce_double_allpool(PARAM.inp.kpar, GlobalV::NPROC_IN_POOL, &(this->Mi_[0][0]), 3 * this->Mi_.size());

    ModuleBase::timer::end("spinconstrain::SpinConstrain", "cal_mi_pw");
}

/// @brief Set the DeltaSpin operator pointer for LCAO magnetic moment calculation
template <>
void spinconstrain::SpinConstrain<std::complex<double>>::set_operator(
    hamilt::Operator<std::complex<double>>* op_in)
{
    this->p_operator = op_in;
}

/// @brief Set the DeltaSpin operator pointer (double specialization)
template <>
void spinconstrain::SpinConstrain<double>::set_operator(
    hamilt::Operator<double>* op_in)
{
    this->p_operator = op_in;
}
