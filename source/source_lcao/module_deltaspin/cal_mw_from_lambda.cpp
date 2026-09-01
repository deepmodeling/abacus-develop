#include "mi_tools.h"
#include "source_base/global_variable.h"
#include "source_base/parallel_comm.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_base/tool_title.h"
#include "source_estate/elecstate_tools.h"
#include "source_hsolver/diag_comm_info.h"
#include "source_hsolver/diago_iter_assist.h"
#include "source_hsolver/hsolver_lcao.h"
#include "source_io/module_parameter/parameter.h"
#include "source_pw/module_pwdft/onsite_proj.h"
#include "spin_constrain.h"

#ifdef __LCAO
#include "source_estate/elecstate_lcao.h"
#include "source_estate/module_dm/cal_dm_psi.h"
#include "source_lcao/module_operator_lcao/dspin_lcao.h"
#endif

/**
 * @file cal_mw_from_lambda.cpp
 * @brief Core computational functions for DeltaSpin.
 *
 * @par calculate_delta_hcc
 * Computes the DeltaSpin correction to the subspace Hamiltonian:
 *   H_corrected = H_original + becp^† * delta_lambda * becp
 *
 * For npol=2 (non-collinear), the 2x2 Pauli matrix coefficients are:
 *   coeff0 = (lambda_z, 0)        coeff1 = (lambda_x, lambda_y)
 *   coeff2 = (lambda_x, -lambda_y) coeff3 = (-lambda_z, 0)
 * Applied as: ps_up = coeff0 * becp_up + coeff2 * becp_dn
 *             ps_dn = coeff1 * becp_up + coeff3 * becp_dn
 *
 * For npol=1 (collinear), only the z-component:
 *   ps = lambda_z * spin_sign * becp
 *
 * @par update_psi_charge_pw_cpu/gpu
 * Two-stage process for PW basis:
 *   1. Subspace diagonalization: apply DeltaSpin correction, rotate psi
 *   2. Full-space update: either run HSolverPW (pw_solve=true) or update weights (pw_solve=false)
 *
 * @par cal_mw_from_lambda
 * The central workflow function called repeatedly during lambda optimization:
 *   LCAO: update lambda in operator -> solve HSolverLCAO -> compute Mi
 *   PW: save subspace data (first call) -> apply H correction -> diagonalize in subspace -> compute Mi from becp
 *
 * @par Error conditions
 * - assert(sub_h_save != nullptr): cal_mw_from_lambda() must be called before
 *   update_psi_charge_pw(). Failure means the workflow order is wrong.
 *   Solution: Ensure cal_mw_from_lambda() is called at the start of each SCF step.
 */

// calculate_delta_hcc() has been moved to
// source/source_pw/module_pwdft/deltaspin_pw_impl.cpp
// because it is PW-specific (operates on PW projector subspace Hamiltonian
// and is only called by update_psi_charge_pw_cpu/gpu).

// update_psi_charge_pw_cpu/gpu() have been moved to
// source/source_pw/module_pwdft/deltaspin_pw_impl.cpp
// because they depend on PW-specific HSolverPW, ElecStatePW, and
// OnsiteProjector.

/**
 * @brief Core workflow: apply lambda -> solve Hamiltonian -> compute magnetic moments.
 *
 * @par LCAO path:
 *   1. Update lambda in DeltaSpin operator (dspin->update_lambda())
 *   2. Solve HSolverLCAO with charge update disabled (last param = true means no charge update)
 *   3. Calculate weights from new eigenvalues
 *   4. Call cal_mi_lcao() to compute moments from density matrix
 *
 * @par PW path:
 *   1. [First call only, i_step==-1] Save subspace H, S, becp from Hamiltonian
 *      This captures the "unperturbed" state before any lambda is applied.
 *   2. [i_step!=-1] Apply DeltaSpin correction via calculate_delta_hcc()
 *      For the first call (i_step==-1), no correction is applied (lambda=0).
 *   3. Diagonalize in subspace via diag_responce(), update becp coefficients
 *   4. Calculate weights from new eigenvalues
 *   5. Call accumulate_Mi_from_becp() for each k-point to compute Mi
 *   6. MPI reduce Mi across k-pools (each pool has a partial sum)
 *
 * @param i_step Current inner lambda step (-1 = initialization, 0+ = optimization)
 * @param delta_lambda Change in lambda from previous step (unused in this function,
 *                     the full lambda_ is used for H correction)
 */
template <>
void spinconstrain::SpinConstrain<std::complex<double>>::cal_mw_from_lambda(
		int i_step,
		const ModuleBase::Vector3<double>* delta_lambda)
{
    ModuleBase::TITLE("spinconstrain::SpinConstrain", "cal_mw_from_lambda");
    ModuleBase::timer::start("spinconstrain::SpinConstrain", "cal_mw_from_lambda");
#ifdef __MPI
    const hsolver::diag_comm_info diag_comm(POOL_WORLD, GlobalV::RANK_IN_POOL, GlobalV::NPROC_IN_POOL);
#else
    const hsolver::diag_comm_info diag_comm(0, 1);
#endif

#ifdef __LCAO
    if (PARAM.inp.basis_type == "lcao")
    {
        // =============================================================
        // LCAO PATH: Update lambda in operator, solve, compute Mi
        // =============================================================
        psi::Psi<std::complex<double>>* psi_t = static_cast<psi::Psi<std::complex<double>>*>(this->psi);
        hamilt::Hamilt<std::complex<double>>* hamilt_t = static_cast<hamilt::Hamilt<std::complex<double>>*>(this->p_hamilt);
        hsolver::HSolverLCAO<std::complex<double>> hsolver_t(this->ParaV,
                                                             PARAM.inp.ks_solver,
                                                             PARAM.globalv.kpar_lcao,
                                                             PARAM.globalv.nlocal,
                                                             PARAM.inp.nbands,
                                                             PARAM.inp.nelec,
                                                             PARAM.inp.device == "gpu",
                                                             GlobalV::NPROC,
                                                             GlobalV::MY_RANK);
        if (this->nspin_ == 2)
        {
            dynamic_cast<hamilt::DeltaSpin<hamilt::OperatorLCAO<std::complex<double>, double>>*>(this->p_operator)
                ->update_lambda();
        }
        else if (this->nspin_ == 4)
        {
            dynamic_cast<hamilt::DeltaSpin<hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>>*>(
                this->p_operator)
                ->update_lambda();
        }
        // Diagonalization without updating charge density (last param = true means skip charge update)
        hsolver_t.solve(hamilt_t, psi_t[0], this->pelec, *this->dm_, *this->pelec->charge, this->nspin_, true);
        elecstate::calculate_weights(this->pelec->ekb,
                                     this->pelec->wg,
                                     this->pelec->klist,
                                     this->pelec->eferm,
                                     this->pelec->f_en,
                                     this->pelec->nelec_spin,
                                     PARAM.inp.nbands,
                                     this->pelec->skip_weights);
        elecstate::calEBand(this->pelec->ekb,this->pelec->wg,this->pelec->f_en);

        // Note: although update_lambda() modifies lambda in-place above,
        // solve() unconditionally recomputes DM and DMR (via cal_dm_psi +
        // cal_DMR) from the psi obtained by diagonalizing with the new
        // lambda. Therefore the DMR used inside cal_mi_lcao() is consistent
        // with the updated lambda and is NOT stale.
        this->cal_mi_lcao(i_step);
    }
    else
#endif
    {
        {
            this->zero_Mi();
            int size_becp = 0;
            std::vector<std::complex<double>> becp_tmp;
            int nk = 0;
            int nkb = 0;
            int nbands = 0;
            int npol = 0;
            const int* nh_iat = nullptr;
            if (PARAM.inp.device == "cpu")
            {
                // =============================================================
                // PW PATH (CPU): Subspace diagonalization + Mi from becp
                // =============================================================
                psi::Psi<std::complex<double>>* psi_t = static_cast<psi::Psi<std::complex<double>>*>(this->psi);
                hamilt::Hamilt<std::complex<double>, base_device::DEVICE_CPU>* hamilt_t = static_cast<hamilt::Hamilt<std::complex<double>, base_device::DEVICE_CPU>*>(this->p_hamilt);
                auto* onsite_p = projectors::OnsiteProjector<double, base_device::DEVICE_CPU>::get_instance();
                nbands = psi_t->get_nbands();
                npol = psi_t->get_npol();
                nkb = onsite_p->get_tot_nproj();
                nk = psi_t->get_nk();
                nh_iat = &onsite_p->get_nh(0);
                size_becp = nbands * nkb * npol;
                becp_tmp.resize(size_becp * nk);
                std::vector<std::complex<double>> h_tmp(nbands * nbands), s_tmp(nbands * nbands);
                int initial_hs = 0;
                if(this->sub_h_save == nullptr)
                {
                    // FIRST CALL: save subspace data for reuse across lambda steps
                    initial_hs = 1;
                    this->sub_h_save = new std::complex<double>[nbands * nbands * nk];
                    this->sub_s_save = new std::complex<double>[nbands * nbands * nk];
                    this->becp_save = new std::complex<double>[size_becp * nk];
                    this->lambda_in_sub_ = this->lambda_;
                }
                for (int ik = 0; ik < nk; ++ik)
                {

                    psi_t->fix_k(ik);

                    std::complex<double>* h_k = this->sub_h_save + ik * nbands * nbands;
                    std::complex<double>* s_k = this->sub_s_save + ik * nbands * nbands;
                    std::complex<double>* becp_k = this->becp_save + ik * size_becp;
                    if(initial_hs)
                    {
                        /// Compute H(k) and extract subspace matrices for this k-point
                        hamilt_t->updateHk(ik);
                        hsolver::DiagoIterAssist<std::complex<double>>::cal_hs_subspace(hamilt_t,
                                                                                        psi_t[0],
                                                                                        h_k,
                                                                                        s_k,
                                                                                        diag_comm);
                        memcpy(becp_k, onsite_p->get_becp(), sizeof(std::complex<double>) * size_becp);
                    }
                    memcpy(h_tmp.data(), h_k, sizeof(std::complex<double>) * nbands * nbands);
                    memcpy(s_tmp.data(), s_k, sizeof(std::complex<double>) * nbands * nbands);
                    // Apply DeltaSpin correction (skip for initialization step i_step=-1)
                    if (i_step != -1) this->calculate_delta_hcc(h_tmp.data(), becp_k, this->lambda_.data(), nbands, nkb, nh_iat, ik, true);

                    // Diagonalize in subspace, update becp (response wavefunctions)
                    hsolver::DiagoIterAssist<std::complex<double>>::diag_responce(h_tmp.data(),
                                                                                  s_tmp.data(),
                                                                                  nbands,
                                                                                  becp_k,
                                                                                  &becp_tmp[ik * size_becp],
                                                                                  nkb * npol,
                                                                                  &this->pelec->ekb(ik, 0));
                }
            }
#if ((defined __CUDA) || (defined __ROCM))
            else
            {
                // =============================================================
                // PW PATH (GPU): Same as CPU but with GPU memory management
                // =============================================================
                psi::Psi<std::complex<double>, base_device::DEVICE_GPU>* psi_t = static_cast<psi::Psi<std::complex<double>, base_device::DEVICE_GPU>*>(this->psi);
                hamilt::Hamilt<std::complex<double>, base_device::DEVICE_GPU>* hamilt_t = static_cast<hamilt::Hamilt<std::complex<double>, base_device::DEVICE_GPU>*>(this->p_hamilt);
                auto* onsite_p = projectors::OnsiteProjector<double, base_device::DEVICE_GPU>::get_instance();
                nbands = psi_t->get_nbands();
                npol = psi_t->get_npol();
                nkb = onsite_p->get_tot_nproj();
                nk = psi_t->get_nk();
                nh_iat = &onsite_p->get_nh(0);
                size_becp = nbands * nkb * npol;
                becp_tmp.resize(size_becp * nk);
                std::complex<double>* h_tmp = nullptr;
                std::complex<double>* s_tmp = nullptr;
                base_device::memory::resize_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(h_tmp, nbands * nbands);
                base_device::memory::resize_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(s_tmp, nbands * nbands);
                int initial_hs = 0;
                if(this->sub_h_save == nullptr)
                {
                    initial_hs = 1;
                    base_device::memory::resize_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(this->sub_h_save, nbands * nbands * nk);
                    base_device::memory::resize_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(this->sub_s_save, nbands * nbands * nk);
                    base_device::memory::resize_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(this->becp_save, size_becp * nk);
                    this->lambda_in_sub_ = this->lambda_;
                }
                std::complex<double>* becp_pointer = nullptr;
                base_device::memory::resize_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(becp_pointer, size_becp);
                for (int ik = 0; ik < nk; ++ik)
                {
                    psi_t->fix_k(ik);

                    std::complex<double>* h_k = this->sub_h_save + ik * nbands * nbands;
                    std::complex<double>* s_k = this->sub_s_save + ik * nbands * nbands;
                    std::complex<double>* becp_k = this->becp_save + ik * size_becp;
                    if(initial_hs)
                    {
                        hamilt_t->updateHk(ik);
                        hsolver::DiagoIterAssist<std::complex<double>, base_device::DEVICE_GPU>::cal_hs_subspace(
                            hamilt_t,
                            psi_t[0],
                            h_k,
                            s_k,
                            diag_comm);
                        base_device::memory::synchronize_memory_op<std::complex<double>, base_device::DEVICE_GPU, base_device::DEVICE_GPU>()(becp_k, onsite_p->get_becp(), size_becp);
                    }
                    base_device::memory::synchronize_memory_op<std::complex<double>, base_device::DEVICE_GPU, base_device::DEVICE_GPU>()(h_tmp, h_k, nbands * nbands);
                    base_device::memory::synchronize_memory_op<std::complex<double>, base_device::DEVICE_GPU, base_device::DEVICE_GPU>()(s_tmp, s_k, nbands * nbands);
                    if (i_step != -1) this->calculate_delta_hcc(h_tmp, becp_k, this->lambda_.data(), nbands, nkb, nh_iat, ik, true);

                    hsolver::DiagoIterAssist<std::complex<double>, base_device::DEVICE_GPU>::diag_responce(h_tmp,
                                                                                  s_tmp,
                                                                                  nbands,
                                                                                  becp_k,
                                                                                  becp_pointer,
                                                                                  nkb * npol,
                                                                                  &this->pelec->ekb(ik, 0));
                    base_device::memory::synchronize_memory_op<std::complex<double>, base_device::DEVICE_CPU, base_device::DEVICE_GPU>()(&becp_tmp[ik * size_becp], becp_pointer, size_becp);
                }

                base_device::memory::delete_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(becp_pointer);
                base_device::memory::delete_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(h_tmp);
                base_device::memory::delete_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(s_tmp);
            }
#endif

            // Calculate weights from eigenvalues to update occupation
            elecstate::calculate_weights(this->pelec->ekb,
                                         this->pelec->wg,
                                         this->pelec->klist,
                                         this->pelec->eferm,
                                         this->pelec->f_en,
                                         this->pelec->nelec_spin,
                                         PARAM.inp.nbands,
                                         this->pelec->skip_weights);
            // Calculate Mi from becp coefficients for each k-point
            for (int ik = 0; ik < nk; ik++)
            {
                const std::complex<double>* becp = &becp_tmp[ik * size_becp];
                const int spin_sign = (this->npol_ == 2) ? 1 : this->get_spin_sign(ik);
                accumulate_Mi_from_becp(becp, nkb, nbands, this->npol_, spin_sign,
                    &this->pelec->wg(ik, 0), nh_iat, this->Mi_);
            }
            // MPI reduction: sum Mi across all k-pool ranks
            Parallel_Reduce::reduce_double_allpool(PARAM.inp.kpar,
                                                    GlobalV::NPROC_IN_POOL,
                                                    &(this->Mi_[0][0]),
                                                    3 * this->Mi_.size());
        }
    }
    ModuleBase::timer::end("spinconstrain::SpinConstrain", "cal_mw_from_lambda");
}

/**
 * @brief Dispatcher: route to LCAO or PW (CPU/GPU) wavefunction/charge update.
 *
 * @details For LCAO: simply calls psiToRho() since the Hamiltonian already
 * includes the DeltaSpin correction.
 * For PW: calls update_psi_charge_pw_cpu or update_psi_charge_pw_gpu
 * which perform subspace diagonalization and optional full-space refinement.
 */
template <>
void spinconstrain::SpinConstrain<std::complex<double>>::update_psi_charge(const ModuleBase::Vector3<double>* delta_lambda, bool pw_solve, bool full_update)
{
    ModuleBase::TITLE("spinconstrain::SpinConstrain", "update_psi_charge");
    ModuleBase::timer::start("spinconstrain::SpinConstrain", "update_psi_charge");
#ifdef __LCAO
    if (PARAM.inp.basis_type == "lcao")
    {
        // TODO: Known issue — the base-class psiToRho() is a no-op for
        // LCAO, so the charge density rho is NOT recomputed here after the
        // final lambda update. After the lambda loop converges, rho remains
        // from the last solve() call with skip_charge=false, which was
        // computed with a different lambda. To fix this, update_psi_charge
        // should recalculate DM from the current psi/weights, then DMR,
        // then rho via dm2rho (similar to the PW path which performs a
        // subspace diagonalization + optional full solve). The DMR inside
        // cal_mi_lcao() itself is fresh (see comment above), but the
        // charge density fed back into the next SCF iteration is stale.
        psi::Psi<std::complex<double>>* psi_t = static_cast<psi::Psi<std::complex<double>>*>(this->psi);
        this->pelec->psiToRho(*psi_t);
    }
    else
#endif
    {
        if (PARAM.inp.device == "cpu")
        {
            this->update_psi_charge_pw_cpu(delta_lambda, pw_solve, full_update);
        }
#if ((defined __CUDA) || (defined __ROCM))
        else
        {
            this->update_psi_charge_pw_gpu(delta_lambda, pw_solve, full_update);
        }
#endif
    }
    ModuleBase::timer::end("spinconstrain::SpinConstrain", "update_psi_charge");
}
