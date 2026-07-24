#include "source_base/timer.h"
#include "source_base/tool_title.h"
#include "source_base/global_variable.h"
#include "source_hsolver/diago_iter_assist.h"
#include "source_io/module_parameter/parameter.h"
#include "spin_constrain.h"
#include "source_pw/module_pwdft/onsite_proj.h"
#include "source_base/parallel_reduce.h"
#include "source_base/kernels/math_kernel_op.h"
#include "source_hsolver/hsolver_lcao.h"
#include "source_hsolver/hsolver_pw.h"
#include "source_estate/elecstate_pw.h"
#include "source_estate/elecstate_tools.h"
#include "source_basis/module_ao/parallel_orbitals.h"

#ifdef __LCAO
#include "source_estate/elecstate_lcao.h"
#include "source_estate/elecstate_tools.h"
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

namespace spinconstrain {
extern void gather_sub_matrix_to_all(const std::complex<double>* local_data,
                                     const Parallel_Orbitals* ParaV,
                                     std::complex<double>* full_data,
                                     int nbands);

extern void scatter_sub_matrix_to_local(const std::complex<double>* full_data,
                                        const Parallel_Orbitals* ParaV,
                                        std::complex<double>* local_data,
                                        int nbands);

extern void extract_diagonal_from_local_block(const std::complex<double>* local_data,
                                              const Parallel_Orbitals* ParaV,
                                              double* diag,
                                              int nbands);
}

/**
 * @brief Compute DeltaSpin correction to the subspace Hamiltonian.
 *
 * @details Adds the constraint term to H in the projector subspace:
 *   H += becp^† * ps, where ps = delta_lambda * becp
 *
 * Code paths by spin type:
 * ---------------------------------------------------------------
 * nspin=4 (npol=2, non-collinear): full Pauli matrix treatment
 *   H_delta = | lambda_z      lambda_x + i*lambda_y |
 *             | lambda_x - i*lambda_y   -lambda_z   |
 *   The 2x2 Pauli matrix couples spin-up and spin-down channels,
 *   allowing rotation of the magnetization direction in 3D space.
 *   projector coefficients are mixed: ps_up = c0*becp_up + c2*becp_dn
 *
 * nspin=2 (npol=1, collinear): diagonal z-component only
 *   H_delta = lambda_z * spin_sign
 *   where spin_sign = +1 for spin-up k-points, -1 for spin-down k-points.
 *   This is a simple scalar correction: different sign per spin channel.
 *   No spin-flip terms exist in collinear mode.
 *
 * @param h_tmp Subspace Hamiltonian (nbands x nbands, modified in place)
 * @param becp_k Projector coefficients for k-point ik
 * @param delta_lambda Lambda change per atom (or full lambda if full_update)
 * @param nbands Number of bands
 * @param nkb Total number of projectors
 * @param nh_iat Number of projectors per atom
 * @param ik K-point index (for spin_sign lookup in collinear mode)
 * @param full_update If true, compute delta = lambda_current - lambda_at_save
 */
template <>
void spinconstrain::SpinConstrain<std::complex<double>>::calculate_delta_hcc(std::complex<double>* h_tmp, const std::complex<double>* becp_k, const ModuleBase::Vector3<double>* delta_lambda, const int nbands, const int nkb, const int* nh_iat, const int ik, bool full_update)
{
    ModuleBase::TITLE("spinconstrain::SpinConstrain", "calculate_delta_hcc");
    ModuleBase::timer::start("spinconstrain::SpinConstrain", "calculate_delta_hcc");

    // If full_update, compute actual delta = lambda_current - lambda_at_save
    // This applies only the CHANGE in lambda, not the full lambda value
    std::vector<ModuleBase::Vector3<double>> actual_delta;
    const ModuleBase::Vector3<double>* effective_lambda = delta_lambda;
    if (full_update)
    {
        int nat = this->get_nat();
        actual_delta.resize(nat);
        for (int iat = 0; iat < nat; iat++)
        {
            actual_delta[iat] = delta_lambda[iat] - this->lambda_in_sub_[iat];
        }
        effective_lambda = actual_delta.data();
    }

    int sum = 0; // Running sum of projectors across atoms
    int size_ps = nkb * this->npol_ * nbands; // Total size of ps array
    std::complex<double>* becp_cpu = nullptr;

    // Handle GPU/CPU memory for becp
    if(PARAM.inp.device == "gpu")
    {
#if ((defined __CUDA) || (defined __ROCM))
        base_device::memory::resize_memory_op<std::complex<double>, base_device::DEVICE_CPU>()(becp_cpu, size_ps);
        base_device::memory::synchronize_memory_op<std::complex<double>, base_device::DEVICE_CPU, base_device::DEVICE_GPU>()(becp_cpu, becp_k, size_ps);
#endif
    }
    else if (PARAM.inp.device == "cpu")
    {
        becp_cpu = const_cast<std::complex<double>*>(becp_k);
    }

    // Compute modified projector coefficients: ps = delta_lambda * becp
    std::vector<std::complex<double>> ps(size_ps, 0.0);
    if(this->npol_ == 2)
    {
        // =============================================================
        // nspin=4 (non-collinear): full Pauli matrix treatment
        // =============================================================
        // For each atom, construct 2x2 coefficients:
        //   | lambda_z      lambda_x + i*lambda_y |
        //   | lambda_x - i*lambda_y   -lambda_z   |
        // Then: ps_up = coeff0 * becp_up + coeff2 * becp_dn
        //        ps_dn = coeff1 * becp_up + coeff3 * becp_dn
        //
        // The spin-flip terms (coeff1, coeff2) couple spin-up and spin-down
        // projector coefficients, enabling rotation of the magnetization
        // direction in the full 3D space.
        for (int iat = 0; iat < this->Mi_.size(); iat++)
        {
            const int nproj = nh_iat[iat];
            const std::complex<double> coefficients0(effective_lambda[iat][2], 0.0);
            const std::complex<double> coefficients1(effective_lambda[iat][0] , effective_lambda[iat][1]);
            const std::complex<double> coefficients2(effective_lambda[iat][0] , -1 * effective_lambda[iat][1]);
            const std::complex<double> coefficients3(-1 * effective_lambda[iat][2], 0.0);
            for (int ib = 0; ib < nbands * this->npol_; ib += this->npol_)
            {
                for (int ip = 0; ip < nproj; ip++)
                {
                    const int becpind = ib * nkb + sum + ip;
                    const std::complex<double> becp1 = becp_cpu[becpind];
                    const std::complex<double> becp2 = becp_cpu[becpind + nkb];
                    ps[becpind] += coefficients0 * becp1
                                    + coefficients2 * becp2;
                    ps[becpind + nkb] += coefficients1 * becp1
                                        + coefficients3 * becp2;
                }
            }
            sum += nproj;
        }
    }
    else if(this->npol_ == 1)
    {
        // =============================================================
        // nspin=2 (collinear): only z-component with spin_sign
        // =============================================================
        // In collinear mode, spins are constrained along the z-axis only.
        // The DeltaSpin correction is a simple scalar:
        //   H_delta = lambda_z * spin_sign
        // where spin_sign = +1 for spin-up k-points, -1 for spin-down.
        //
        // This means:
        //   - For spin-up bands:   H += +lambda_z (lowers energy)
        //   - For spin-down bands: H += -lambda_z (raises energy)
        //
        // The sign difference creates an energy splitting between spin-up
        // and spin-down states, which drives the magnetization toward
        // the target value. Positive lambda_z increases M_z (more spin-up
        // occupation, less spin-down).
        //
        // Note: effective_lambda[iat][0] (lambda_x) and [1] (lambda_y)
        // are zero for collinear mode (constrained by init_sc.cpp).
        for (int iat = 0; iat < this->Mi_.size(); iat++)
        {
            const int nproj = nh_iat[iat];
            double coefficients0 = effective_lambda[iat][2] * this->get_spin_sign(ik);
            for (int ib = 0; ib < nbands; ib++)
            {
                for (int ip = 0; ip < nproj; ip++)
                {
                    const int becpind = ib * nkb + sum + ip;
                    const std::complex<double> becp1 = becp_cpu[becpind];
                    ps[becpind] += coefficients0 * becp1;
                }
            }
            sum += nproj;
        }
    }

    // Copy ps to GPU if needed
    std::complex<double>* ps_pointer = nullptr;
    if(PARAM.inp.device == "gpu")
    {
#if ((defined __CUDA) || (defined __ROCM))
        base_device::memory::resize_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(ps_pointer, size_ps);
        base_device::memory::synchronize_memory_op<std::complex<double>, base_device::DEVICE_GPU, base_device::DEVICE_CPU>()(ps_pointer, ps.data(), size_ps);
#endif
    }
    else if (PARAM.inp.device == "cpu")
    {
        ps_pointer = ps.data();
    }

    // =============================================================
    // H += becp^† * ps (GEMM: C = alpha * A^† * B + beta * C)
    // A = becp_k (npm x nbands), B = ps (npm x nbands), C = h_tmp (nbands x nbands)
    // =============================================================
    char transa = 'C'; // Conjugate transpose of becp
    char transb = 'N'; // Normal ps
    const int npm = nkb * this->npol_;
    if (PARAM.inp.device == "gpu")
    {
#if ((defined __CUDA) || (defined __ROCM))
        ModuleBase::gemm_op<std::complex<double>, base_device::DEVICE_GPU>()(
            transa,
            transb,
            nbands,
            nbands,
            npm,
            &ModuleBase::ONE,
            becp_k,
            npm,
            ps_pointer,
            npm,
            &ModuleBase::ONE,
            h_tmp,
            nbands
        );
        base_device::memory::delete_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(ps_pointer);
        base_device::memory::delete_memory_op<std::complex<double>, base_device::DEVICE_CPU>()(becp_cpu);
#endif

    }
    else if (PARAM.inp.device == "cpu")
    {
        ModuleBase::gemm_op<std::complex<double>, base_device::DEVICE_CPU>()(
            transa,
            transb,
            nbands,
            nbands,
            npm,
            &ModuleBase::ONE,
            becp_k,
            npm,
            ps_pointer,
            npm,
            &ModuleBase::ONE,
            h_tmp,
            nbands
        );
    }
    ModuleBase::timer::end("spinconstrain::SpinConstrain", "calculate_delta_hcc");
}

/**
 * @brief CPU implementation of PW wavefunction and charge density update.
 *
 * @par Two-stage process:
 * Stage 1 - Subspace diagonalization:
 *   For each k-point, apply DeltaSpin correction to the saved subspace H,
 *   then diagonalize to rotate the wavefunctions. This is a cheap operation
 *   in the reduced subspace (nbands x nbands).
 *
 * Stage 2 - Full-space update:
 *   Option A (pw_solve=true): Run HSolverPW for iterative refinement in the
 *     full plane-wave space. This is more accurate but expensive.
 *   Option B (pw_solve=false): Update weights from new eigenvalues and call
 *     psiToRho() to build the charge density from current psi. Faster but
 *     may be less accurate if the subspace rotation was not sufficient.
 *
 * @par Memory management
 * Frees sub_h_save, sub_s_save, becp_save after use. These are allocated
 * on the first cal_mw_from_lambda() call and should only be freed here.
 *
 * @param delta_lambda Lambda change for incremental H correction
 * @param pw_solve If true, run full PW solver; if false, just update weights
 * @param full_update If true, apply full lambda (not delta) to H correction
 */
template <>
void spinconstrain::SpinConstrain<std::complex<double>>::update_psi_charge_pw_cpu(const ModuleBase::Vector3<double>* delta_lambda, bool pw_solve, bool full_update)
{
    ModuleBase::TITLE("spinconstrain::SpinConstrain", "update_psi_charge_pw_cpu");
    ModuleBase::timer::start("spinconstrain::SpinConstrain", "update_psi_charge_pw_cpu");

    psi::Psi<std::complex<double>>* psi_t = static_cast<psi::Psi<std::complex<double>>*>(this->psi);
    hamilt::Hamilt<std::complex<double>, base_device::DEVICE_CPU>* hamilt_t = static_cast<hamilt::Hamilt<std::complex<double>, base_device::DEVICE_CPU>*>(this->p_hamilt);
    auto* onsite_p = projectors::OnsiteProjector<double, base_device::DEVICE_CPU>::get_instance();

    int nbands = psi_t->get_nbands();
    int npol = psi_t->get_npol();
    int nkb = onsite_p->get_tot_nproj();
    int nk = psi_t->get_nk();
    int size_becp = nbands * nkb * npol;
    const int* nh_iat = &onsite_p->get_nh(0);

    std::vector<std::complex<double>> h_tmp(nbands * nbands), s_tmp(nbands * nbands);

    // CRITICAL: subspace data must have been saved by cal_mw_from_lambda()
    assert(this->sub_h_save != nullptr);
    assert(this->sub_s_save != nullptr);
    assert(this->becp_save != nullptr);

    // Determine which lambda to use for H correction
    const ModuleBase::Vector3<double>* lambda_for_hcc = delta_lambda;
    std::vector<ModuleBase::Vector3<double>> computed_delta;
    if (full_update)
    {
        lambda_for_hcc = this->lambda_.data();
    }

    // =============================================================
    // STAGE 1: Subspace diagonalization for each k-point
    // =============================================================
    for (int ik = 0; ik < nk; ++ik)
    {
        std::complex<double>* h_k = this->sub_h_save + ik * nbands * nbands;
        std::complex<double>* s_k = this->sub_s_save + ik * nbands * nbands;
        std::complex<double>* becp_k = this->becp_save + ik * size_becp;

        psi_t->fix_k(ik);

        // Copy saved subspace matrices to temp
        memcpy(h_tmp.data(), h_k, sizeof(std::complex<double>) * nbands * nbands);
        memcpy(s_tmp.data(), s_k, sizeof(std::complex<double>) * nbands * nbands);

        // Apply DeltaSpin correction: H += becp^† * lambda * becp
        this->calculate_delta_hcc(h_tmp.data(), becp_k, lambda_for_hcc, nbands, nkb, nh_iat, ik, full_update);

        // Diagonalize in subspace to update wavefunction coefficients and eigenvalues
        hsolver::DiagoIterAssist<std::complex<double>>::diag_subspace_psi(h_tmp.data(),
                                                                        s_tmp.data(),
                                                                        nbands,
                                                                        psi_t[0],
                                                                        &this->pelec->ekb(ik, 0));
    }

    // Free saved subspace data (allocated in cal_mw_from_lambda)
    delete[] this->sub_h_save;
    delete[] this->sub_s_save;
    delete[] this->becp_save;
    this->sub_h_save = nullptr;
    this->sub_s_save = nullptr;
    this->becp_save = nullptr;

    // =============================================================
    // STAGE 2: Full-space update
    // =============================================================
    if (pw_solve)
    {
        // Full PW diagonalization: subspace rotation provides a good initial guess,
        // then HSolverPW iteratively refines psi in the full plane-wave space and calls psiToRho.
        hsolver::HSolverPW<std::complex<double>, base_device::DEVICE_CPU> hsolver_pw_obj(
            this->pw_wfc_,
            PARAM.inp.calculation,
            PARAM.inp.basis_type,
            PARAM.inp.ks_solver,
            PARAM.globalv.use_uspp,
            PARAM.inp.nspin,
            hsolver::DiagoIterAssist<std::complex<double>>::SCF_ITER,
            hsolver::DiagoIterAssist<std::complex<double>>::PW_DIAG_NMAX,
            hsolver::DiagoIterAssist<std::complex<double>>::PW_DIAG_THR,
            hsolver::DiagoIterAssist<std::complex<double>>::need_subspace,
            PARAM.inp.use_k_continuity);

        hsolver_pw_obj.solve(hamilt_t, psi_t[0], this->pelec, this->pelec->ekb.c,
            GlobalV::RANK_IN_POOL, GlobalV::NPROC_IN_POOL, false, this->tpiba, this->get_nat());
    }
    else
    {
        // No full solver: update weights from new eigenvalues, then build rho from current psi
        elecstate::calculate_weights(this->pelec->ekb,
                                     this->pelec->wg,
                                     this->pelec->klist,
                                     this->pelec->eferm,
                                     this->pelec->f_en,
                                     this->pelec->nelec_spin,
                                     this->pelec->skip_weights);
        elecstate::calEBand(this->pelec->ekb, this->pelec->wg, this->pelec->f_en);
        reinterpret_cast<elecstate::ElecStatePW<std::complex<double>, base_device::DEVICE_CPU>*>(this->pelec)->psiToRho(*psi_t);
    }
    ModuleBase::timer::end("spinconstrain::SpinConstrain", "update_psi_charge_pw_cpu");
}

#if ((defined __CUDA) || (defined __ROCM))
/**
 * @brief GPU implementation of PW wavefunction and charge density update.
 *
 * @details Same algorithm as update_psi_charge_pw_cpu(), but with GPU memory
 * management (device allocation, host-device synchronization).
 */
template <>
void spinconstrain::SpinConstrain<std::complex<double>>::update_psi_charge_pw_gpu(const ModuleBase::Vector3<double>* delta_lambda, bool pw_solve, bool full_update)
{
    ModuleBase::TITLE("spinconstrain::SpinConstrain", "update_psi_charge_pw_gpu");
    ModuleBase::timer::start("spinconstrain::SpinConstrain", "update_psi_charge_pw_gpu");

    psi::Psi<std::complex<double>, base_device::DEVICE_GPU>* psi_t = static_cast<psi::Psi<std::complex<double>, base_device::DEVICE_GPU>*>(this->psi);
    hamilt::Hamilt<std::complex<double>, base_device::DEVICE_GPU>* hamilt_t = static_cast<hamilt::Hamilt<std::complex<double>, base_device::DEVICE_GPU>*>(this->p_hamilt);
    auto* onsite_p = projectors::OnsiteProjector<double, base_device::DEVICE_GPU>::get_instance();

    int nbands = psi_t->get_nbands();
    int npol = psi_t->get_npol();
    int nkb = onsite_p->get_tot_nproj();
    int nk = psi_t->get_nk();
    int size_becp = nbands * nkb * npol;
    const int* nh_iat = &onsite_p->get_nh(0);

    std::complex<double>* h_tmp = nullptr;
    std::complex<double>* s_tmp = nullptr;
    base_device::memory::resize_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(h_tmp, nbands * nbands);
    base_device::memory::resize_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(s_tmp, nbands * nbands);

    assert(this->sub_h_save != nullptr);
    assert(this->sub_s_save != nullptr);
    assert(this->becp_save != nullptr);

    const ModuleBase::Vector3<double>* lambda_for_hcc = delta_lambda;
    std::vector<ModuleBase::Vector3<double>> computed_delta;
    if (full_update)
    {
        lambda_for_hcc = this->lambda_.data();
    }

    // STAGE 1: Subspace diagonalization for each k-point (GPU)
    for (int ik = 0; ik < nk; ++ik)
    {
        std::complex<double>* h_k = this->sub_h_save + ik * nbands * nbands;
        std::complex<double>* s_k = this->sub_s_save + ik * nbands * nbands;
        std::complex<double>* becp_k = this->becp_save + ik * size_becp;

        psi_t->fix_k(ik);

        base_device::memory::synchronize_memory_op<std::complex<double>, base_device::DEVICE_GPU, base_device::DEVICE_GPU>()(h_tmp, h_k, nbands * nbands);
        base_device::memory::synchronize_memory_op<std::complex<double>, base_device::DEVICE_GPU, base_device::DEVICE_GPU>()(s_tmp, s_k, nbands * nbands);

        this->calculate_delta_hcc(h_tmp, becp_k, lambda_for_hcc, nbands, nkb, nh_iat, ik, full_update);

        hsolver::DiagoIterAssist<std::complex<double>, base_device::DEVICE_GPU>::diag_subspace_psi(h_tmp,
                                                                                s_tmp,
                                                                                nbands,
                                                                                psi_t[0],
                                                                                &this->pelec->ekb(ik, 0));
    }

    // Free GPU memory for saved subspace data
    base_device::memory::delete_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(sub_h_save);
    base_device::memory::delete_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(sub_s_save);
    base_device::memory::delete_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(becp_save);
    this->sub_h_save = nullptr;
    this->sub_s_save = nullptr;
    this->becp_save = nullptr;

    // STAGE 2: Full-space update (GPU)
    if (pw_solve)
    {
        hsolver::HSolverPW<std::complex<double>, base_device::DEVICE_GPU> hsolver_pw_obj(
            this->pw_wfc_,
            PARAM.inp.calculation,
            PARAM.inp.basis_type,
            PARAM.inp.ks_solver,
            PARAM.globalv.use_uspp,
            PARAM.inp.nspin,
            hsolver::DiagoIterAssist<std::complex<double>, base_device::DEVICE_GPU>::SCF_ITER,
            hsolver::DiagoIterAssist<std::complex<double>, base_device::DEVICE_GPU>::PW_DIAG_NMAX,
            hsolver::DiagoIterAssist<std::complex<double>, base_device::DEVICE_GPU>::PW_DIAG_THR,
            hsolver::DiagoIterAssist<std::complex<double>, base_device::DEVICE_GPU>::need_subspace,
            PARAM.inp.use_k_continuity);

        hsolver_pw_obj.solve(hamilt_t, psi_t[0], this->pelec, this->pelec->ekb.c,
            GlobalV::RANK_IN_POOL, GlobalV::NPROC_IN_POOL, false, this->tpiba, this->get_nat());
    }
    else
    {
        elecstate::calculate_weights(this->pelec->ekb,
                                     this->pelec->wg,
                                     this->pelec->klist,
                                     this->pelec->eferm,
                                     this->pelec->f_en,
                                     this->pelec->nelec_spin,
                                     this->pelec->skip_weights);
        elecstate::calEBand(this->pelec->ekb, this->pelec->wg, this->pelec->f_en);
        reinterpret_cast<elecstate::ElecStatePW<std::complex<double>, base_device::DEVICE_GPU>*>(this->pelec)->psiToRho(*psi_t);
    }
    ModuleBase::timer::end("spinconstrain::SpinConstrain", "update_psi_charge_pw_gpu");
}
#endif

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

#ifdef __LCAO
    if (PARAM.inp.basis_type == "lcao")
    {
        // =============================================================
        // LCAO PATH: Update lambda in operator, solve, compute Mi
        // =============================================================
        //
        // This function is called from the BFGS lambda optimization loop
        // (run_lambda_loop) to compute magnetic moments Mi for a given lambda.
        //
        // i_step semantics:
        //   i_step = -1: Initialization step from main loop (compute initial Mi)
        //   i_step =  0: First BFGS iteration (delta_lambda != 0)
        //   i_step =  1,2,...: Subsequent BFGS iterations
        //   i_step = -2: SPECIAL - build subspace cache at current lambda
        //                (only called by lambda_loop when acceleration activates,
        //                 NOT from the main loop)
        //
        // Three execution paths:
        //   1. i_step == -2: Full diagonalization + subspace cache build
        //   2. accel_enabled + subspace_built: Accelerated (subspace/first-order)
        //   3. else (default): Full HSolverLCAO diagonalization (original behavior)
        //
        // IMPORTANT: Path 3 (full diagonalization) is the DEFAULT and ONLY path
        // used when acceleration is disabled (sc_acceleration_mode = "off").
        // This ensures the original correctness is preserved. The accelerated
        // paths are opt-in and only activate once RMS drops below the threshold.
        // =============================================================
        psi::Psi<std::complex<double>>* psi_t = static_cast<psi::Psi<std::complex<double>>*>(this->psi);
        hamilt::Hamilt<std::complex<double>>* hamilt_t = static_cast<hamilt::Hamilt<std::complex<double>>*>(this->p_hamilt);

        // Check if acceleration is configured (user must set sc_acceleration_mode != "off"
        // and sc_acceleration_rms_thr > 0 in the input file).
        // Note: this only checks configuration; actual activation requires
        // acceleration_active_ to be set by lambda_loop when RMS crosses threshold.
        const bool accel_enabled = (this->sc_acceleration_mode_ != "off") &&
                                   (this->sc_acceleration_rms_thr_ > 0.0) &&
                                   this->acceleration_active_;

        // =================================================================
        // BRANCH 1: i_step == -2 → Build subspace cache at current lambda
        // =================================================================
        //
        // PURPOSE: Called by lambda_loop when acceleration is first activated
        // (RMS drops below sc_acceleration_rms_thr). Builds a reference subspace
        // at the current converged lambda for subsequent accelerated steps.
        //
        // WHY -2? The main BFGS loop uses i_step = -1, 0, 1, ..., nsc-1.
        // Using -2 ensures this branch is NEVER reached from the main loop,
        // only from the explicit cal_mw_from_lambda(-2, ...) call in lambda_loop.
        // This avoids side effects on the main optimization state.
        //
        // WHY full diagonalization? Subspace methods cannot replace the full
        // solver for large lambda steps because they don't correctly update
        // the wavefunction and density matrix. Subspace is only valid for
        // small perturbations near convergence.
        // =================================================================
        if (i_step == -2 && accel_enabled)
        {
            // Step 1: Update DeltaSpin operator with current lambda
            if (this->nspin_ == 2)
            {
                dynamic_cast<hamilt::DeltaSpin<hamilt::OperatorLCAO<std::complex<double>, double>>*>(this->p_operator)
                    ->update_lambda();
            }
            else if (this->nspin_ == 4)
            {
                dynamic_cast<hamilt::DeltaSpin<hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>>* >(
                    this->p_operator)
                    ->update_lambda();
            }

            // Step 2: Full diagonalization (same as default path)
            // Must use full solver to get correct wavefunctions at this lambda
            hsolver::HSolverLCAO<std::complex<double>> hsolver_t(this->ParaV, PARAM.inp.ks_solver);
            hsolver_t.solve(hamilt_t, psi_t[0], this->pelec, *this->dm_, *this->pelec->charge, this->nspin_, true);
            elecstate::calculate_weights(this->pelec->ekb,
                                         this->pelec->wg,
                                         this->pelec->klist,
                                         this->pelec->eferm,
                                         this->pelec->f_en,
                                         this->pelec->nelec_spin,
                                         this->pelec->skip_weights);
            elecstate::calEBand(this->pelec->ekb, this->pelec->wg, this->pelec->f_en);

            // Step 3: Cache subspace data for subsequent accelerated steps
            this->free_lcao_subspace_cache();
            const int nk = psi_t->get_nk();
            const int nbands = PARAM.inp.nbands;
            const int nlocal = this->ParaV->get_global_row_size();
            const int nn = nbands * nbands;
            this->lcao_nbands_ = nbands;
            this->lcao_nk_ = nk;
            this->lcao_nlocal_ = nlocal;

            this->lcao_sub_h_save = new std::complex<double>[nk * nn];
            this->lcao_sub_s_save = new std::complex<double>[nk * nn];
            this->lcao_PI_sub_save_.resize(nk);
            this->lcao_PI_sub_diag_.resize(nk);
            this->lcao_ekb_save_.resize(nk * nbands);

            const int nloc_eij = this->ParaV->nrow * this->ParaV->ncol_bands;
            const int nat = this->get_nat();

            for (int ik = 0; ik < nk; ik++)
            {
                psi_t->fix_k(ik);

                this->calculate_lcao_sub_hs(
                    this->p_hamilt, psi_t[0], this->ParaV,
                    this->lcao_sub_h_save + ik * nn,
                    this->lcao_sub_s_save + ik * nn,
                    ik, nbands, nlocal);

                if (this->nspin_ == 2)
                {
                    auto* dspin_op = dynamic_cast<hamilt::DeltaSpin<hamilt::OperatorLCAO<std::complex<double>, double>>*>(
                        this->p_operator);
                    for (int iat = 0; iat < nat; iat++)
                    {
                        if (!dspin_op->get_constraint_atom_list()[iat])
                        {
                            continue;
                        }
                        this->lcao_PI_sub_save_[ik][iat].resize(nloc_eij, {0.0, 0.0});
                        this->calculate_PI_sub_from_hr(
                            dspin_op->get_pre_hr(iat),
                            psi_t[0], this->ParaV,
                            this->kv_.kvec_d[ik],
                            this->lcao_PI_sub_save_[ik][iat].data(),
                            nbands, nlocal);

                        this->lcao_PI_sub_diag_[ik][iat].resize(nbands, 0.0);
                        extract_diagonal_from_local_block(
                            this->lcao_PI_sub_save_[ik][iat].data(),
                            this->ParaV,
                            this->lcao_PI_sub_diag_[ik][iat].data(),
                            nbands);
                    }
                }
                else if (this->nspin_ == 4)
                {
                    auto* dspin_op = dynamic_cast<hamilt::DeltaSpin<hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>>* >(
                        this->p_operator);
                    for (int iat = 0; iat < nat; iat++)
                    {
                        if (!dspin_op->get_constraint_atom_list()[iat])
                        {
                            continue;
                        }
                        this->lcao_PI_sub_save_[ik][iat].resize(nloc_eij, {0.0, 0.0});
                        this->calculate_PI_sub_from_hr(
                            dspin_op->get_pre_hr(iat),
                            psi_t[0], this->ParaV,
                            this->kv_.kvec_d[ik],
                            this->lcao_PI_sub_save_[ik][iat].data(),
                            nbands, nlocal);

                        this->lcao_PI_sub_diag_[ik][iat].resize(nbands, 0.0);
                        extract_diagonal_from_local_block(
                            this->lcao_PI_sub_save_[ik][iat].data(),
                            this->ParaV,
                            this->lcao_PI_sub_diag_[ik][iat].data(),
                            nbands);
                    }
                }

                // Save eigenvalues for Fermi weight calculation
                for (int ib = 0; ib < nbands; ib++)
                    this->lcao_ekb_save_[ik * nbands + ib] = this->pelec->ekb(ik, ib);
            }
            // Record the lambda at which subspace was built
            this->lcao_lambda_in_sub_ = this->lambda_;
            this->lcao_subspace_initialized_ = true;

            // Compute Mi via subspace method at lambda_ref (V=I).
            // This ensures subsequent subspace Mi values are consistent with
            // the reference, avoiding mixed full/subspace gradient artifacts.
            {
                std::vector<std::vector<std::complex<double>>> vcc_identity(nk,
                    std::vector<std::complex<double>>(nn, {0.0, 0.0}));
                for (int ik = 0; ik < nk; ik++)
                    for (int ib = 0; ib < nbands; ib++)
                        vcc_identity[ik][ib + ib * nbands] = {1.0, 0.0};

                this->cal_mi_lcao_subspace(vcc_identity, nbands, nk, this->get_npol());
            }

            // Signal BFGS to reset history on next iteration
            this->subspace_just_activated_ = true;
        }
        // =================================================================
        // BRANCH 2: Accelerated path (opt-in, only after RMS drops below threshold)
        // =================================================================
        //
        // PURPOSE: When the BFGS loop is near convergence (RMS < sc_acceleration_rms_thr),
        // lambda changes are small. Subspace methods can approximate the response
        // much faster than full diagonalization while maintaining accuracy.
        //
        // Two modes:
        //   "first_order": Only eigenvalues change (ε_new = ε_old + Δλ · P_sub).
        //                  Wavefunctions are unchanged. Fastest, but least accurate.
        //                  Uses trace formula for Mi from subspace density matrix.
        //
        //   "subspace": Full subspace diagonalization with wavefunction rotation.
        //               H_sub(λ) = H₀_sub + (λ - λ_ref) · P_sub is diagonalized
        //               in the nbands × nbands subspace. Wavefunctions are rotated
        //               to build the full-space density matrix (DMR) for Mi.
        //               More accurate than first-order, still faster than full diag.
        //
        // KEY DESIGN: Both modes compute Mi using the FULL-space density matrix
        // (via cal_dm_psi + cal_DMR + cal_mi_lcao), NOT the subspace trace formula.
        // This is critical because subspace Mi (cal_mi_lcao_subspace) was found
        // to give biased results when wavefunction rotation is not handled correctly.
        //
        // WHY restore psi after subspace rotation? The accelerated path should
        // not modify the wavefunction state that the BFGS loop expects. The
        // rotation is only temporary for Mi calculation. The actual wavefunction
        // update happens in update_psi_charge() after convergence.
        // =================================================================
        else if (accel_enabled && this->acceleration_subspace_built_)
        {
            const int nk = psi_t->get_nk();
            const int nbands = PARAM.inp.nbands;
            const int nrow = this->ParaV->nrow;
            const int nloc_wfc = this->ParaV->nloc_wfc;
            const int nn = nbands * nbands;

            if (this->sc_acceleration_mode_ == "first_order")
            {
                // ---------------------------------------------------------
                // First-order response mode
                // ---------------------------------------------------------
                // Δε_ib = Σ_I Δλ_I · ⟨ψ_ib | D_I | ψ_ib⟩
                // Only eigenvalues shift; wavefunctions are frozen at reference.
                // Valid when Δλ is small enough that wavefunction mixing is negligible.
                //
                // nspin=2: Uses σ_z sign (±1) for up/down spin channels.
                // nspin=4: Uses Δλ_z only (z-component dominates in collinear-like states).
                // ---------------------------------------------------------
                for (int ik = 0; ik < nk; ik++)
                {
                    for (int ib = 0; ib < nbands; ib++)
                    {
                        double delta_epsilon = 0.0;
                        for (const auto& [iat, diag] : this->lcao_PI_sub_diag_[ik])
                        {
                            double p_diag = diag[ib];
                            double dl = this->lambda_[iat].z - this->lcao_lambda_in_sub_[iat].z;
                            delta_epsilon += dl * p_diag;
                        }
                        if (this->nspin_ == 2)
                        {
                            int spin_sign = this->get_spin_sign(ik);
                            this->pelec->ekb(ik, ib) = this->lcao_ekb_save_[ik * nbands + ib]
                                                       - spin_sign * delta_epsilon;
                        }
                        else
                        {
                            this->pelec->ekb(ik, ib) = this->lcao_ekb_save_[ik * nbands + ib]
                                                       - delta_epsilon;
                        }
                    }
                }

                elecstate::calculate_weights(this->pelec->ekb,
                                             this->pelec->wg,
                                             this->pelec->klist,
                                             this->pelec->eferm,
                                             this->pelec->f_en,
                                             this->pelec->nelec_spin,
                                             this->pelec->skip_weights);
                elecstate::calEBand(this->pelec->ekb, this->pelec->wg, this->pelec->f_en);

                // Compute Mi from original (unrotated) psi with updated Fermi weights
                if (this->subspace_exec_precision_ == ModuleGint::GintPrecision::fp32)
                {
                    elecstate::cal_dm_psi_mixed(this->ParaV, this->pelec->wg, *psi_t, *this->dm_);
                }
                else
                {
                    elecstate::cal_dm_psi(this->ParaV, this->pelec->wg, *psi_t, *this->dm_);
                }
                this->dm_->cal_DMR();
                this->cal_mi_lcao(i_step);
            }
            else // sc_acceleration_mode_ == "subspace"
            {
                // ---------------------------------------------------------
                // Subspace diagonalization mode (distributed P_I_sub cache)
                // ---------------------------------------------------------
                // 1. Scatter H₀_sub(k) to local 2D block
                // 2. Accumulate correction: h_sub_local += Σ_I Δλ_I · P_I_sub_local
                // 3. Gather h_sub_local back to full matrix for diagonalization
                // 4. Same for S_sub (unchanged, just copy)
                // ---------------------------------------------------------
                std::vector<std::vector<std::complex<double>>> vcc_all(nk);
                const int nloc_eij = this->ParaV->nrow * this->ParaV->ncol_bands;

                if (this->h_sub_local_buf_.size() != static_cast<std::size_t>(nloc_eij))
                {
                    this->h_sub_local_buf_.resize(nloc_eij, {0.0, 0.0});
                    this->h_tmp_buf_.resize(nn);
                    this->s_tmp_buf_.resize(nn);
                    this->vcc_buf_.resize(nn);
                    this->s_copy_buf_.resize(nn);
                    this->eigenvalues_buf_.resize(nbands, 0.0);
                }

                for (int ik = 0; ik < nk; ik++)
                {
                    std::fill(this->h_sub_local_buf_.begin(), this->h_sub_local_buf_.end(), std::complex<double>(0.0, 0.0));
                    scatter_sub_matrix_to_local(this->lcao_sub_h_save + ik * nn, this->ParaV,
                                                this->h_sub_local_buf_.data(), nbands);

                    this->calculate_delta_hcc_lcao(this->h_sub_local_buf_.data(), this->lcao_PI_sub_save_[ik],
                                                   this->lambda_.data(), nbands, ik, true, this->ParaV);

                    gather_sub_matrix_to_all(this->h_sub_local_buf_.data(), this->ParaV, this->h_tmp_buf_.data(), nbands);

                    std::memcpy(this->s_tmp_buf_.data(), this->lcao_sub_s_save + ik * nn, sizeof(std::complex<double>) * nn);

                    std::memcpy(this->s_copy_buf_.data(), this->s_tmp_buf_.data(), sizeof(std::complex<double>) * nn);

                    hsolver::DiagoIterAssist<std::complex<double>>::diag_hegvd(
                        nbands, nbands, this->h_tmp_buf_.data(), this->s_copy_buf_.data(), nbands,
                        this->eigenvalues_buf_.data(), this->vcc_buf_.data());

                    vcc_all[ik].assign(this->vcc_buf_.data(), this->vcc_buf_.data() + nn);
                    for (int ib = 0; ib < nbands; ib++)
                        this->pelec->ekb(ik, ib) = this->eigenvalues_buf_[ib];
                }

                elecstate::calculate_weights(this->pelec->ekb,
                                             this->pelec->wg,
                                             this->pelec->klist,
                                             this->pelec->eferm,
                                             this->pelec->f_en,
                                             this->pelec->nelec_spin,
                                             this->pelec->skip_weights);
                elecstate::calEBand(this->pelec->ekb, this->pelec->wg, this->pelec->f_en);

                // Rotate psi, build DMR, compute Mi, restore psi
                this->cal_mi_lcao_subspace(vcc_all, nbands, nk, this->get_npol());
            }
        }
        // =================================================================
        // BRANCH 3: DEFAULT PATH — Full HSolverLCAO diagonalization
        // =================================================================
        //
        // PURPOSE: This is the ORIGINAL behavior of the code. It is ALWAYS used
        // when acceleration is disabled (sc_acceleration_mode = "off"), which
        // is the default. It is also used during early BFGS iterations before
        // RMS drops below the acceleration threshold.
        //
        // WHY full diagonalization is mandatory:
        //   - Subspace methods only approximate small perturbations around a
        //     reference lambda. They do NOT correctly handle large lambda steps
        //     where the Hamiltonian changes significantly.
        //   - The full solver rebuilds the Hamiltonian with the updated lambda,
        //     correctly updates the wavefunctions, and produces a consistent
        //     density matrix for the next SCF step.
        //   - Replacing the full solver with subspace diagonalization was the
        //     ROOT CAUSE of the convergence failure (bias ~0.12-0.39 uB) observed
        //     in earlier experiments.
        // =================================================================
        else
        {
            hsolver::HSolverLCAO<std::complex<double>> hsolver_t(this->ParaV, PARAM.inp.ks_solver);
            if (this->nspin_ == 2)
            {
                dynamic_cast<hamilt::DeltaSpin<hamilt::OperatorLCAO<std::complex<double>, double>>*>(this->p_operator)
                    ->update_lambda();
            }
            else if (this->nspin_ == 4)
            {
                dynamic_cast<hamilt::DeltaSpin<hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>>* >(
                    this->p_operator)
                    ->update_lambda();
            }
            // Diagonalization without updating charge density (last param = true)
            // Charge density is only updated after convergence in update_psi_charge()
            hsolver_t.solve(hamilt_t, psi_t[0], this->pelec, *this->dm_, *this->pelec->charge, this->nspin_, true);
            elecstate::calculate_weights(this->pelec->ekb,
                                         this->pelec->wg,
                                         this->pelec->klist,
                                         this->pelec->eferm,
                                         this->pelec->f_en,
                                         this->pelec->nelec_spin,
                                         this->pelec->skip_weights);
            elecstate::calEBand(this->pelec->ekb, this->pelec->wg, this->pelec->f_en);

            this->cal_mi_lcao(i_step);
        }
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
                        hsolver::DiagoIterAssist<std::complex<double>>::cal_hs_subspace(hamilt_t, psi_t[0], h_k, s_k);
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
                        hsolver::DiagoIterAssist<std::complex<double>, base_device::DEVICE_GPU>::cal_hs_subspace(hamilt_t, psi_t[0], h_k, s_k);
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
            }
#endif

            // Calculate weights from eigenvalues to update occupation
            elecstate::calculate_weights(this->pelec->ekb,
                                         this->pelec->wg,
                                         this->pelec->klist,
                                         this->pelec->eferm,
                                         this->pelec->f_en,
                                         this->pelec->nelec_spin,
                                         this->pelec->skip_weights);
            // Calculate Mi from becp coefficients for each k-point
            for (int ik = 0; ik < nk; ik++)
            {
                const std::complex<double>* becp = &becp_tmp[ik * size_becp];
                this->accumulate_Mi_from_becp(becp, nkb, nbands, this->npol_, ik,
                    &this->pelec->wg(ik, 0), nh_iat);
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
        psi::Psi<std::complex<double>>* psi_t = static_cast<psi::Psi<std::complex<double>>*>(this->psi);

        // If subspace acceleration was active, psi is still at the reference lambda,
        // not at the converged lambda. Must do a final full diag to get correct psi
        // for charge density construction.
        if (this->acceleration_active_ && this->acceleration_subspace_built_)
        {
            hamilt::Hamilt<std::complex<double>>* hamilt_t
                = static_cast<hamilt::Hamilt<std::complex<double>>*>(this->p_hamilt);
            if (this->nspin_ == 2)
            {
                dynamic_cast<hamilt::DeltaSpin<hamilt::OperatorLCAO<std::complex<double>, double>>*>(
                    this->p_operator)->update_lambda();
            }
            else if (this->nspin_ == 4)
            {
                dynamic_cast<hamilt::DeltaSpin<hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>>* >(
                    this->p_operator)->update_lambda();
            }
            hsolver::HSolverLCAO<std::complex<double>> hsolver_t(this->ParaV, PARAM.inp.ks_solver);
            hsolver_t.solve(hamilt_t, psi_t[0], this->pelec, *this->dm_,
                            *this->pelec->charge, this->nspin_, true);
            elecstate::calculate_weights(this->pelec->ekb, this->pelec->wg,
                                         this->pelec->klist, this->pelec->eferm,
                                         this->pelec->f_en, this->pelec->nelec_spin,
                                         this->pelec->skip_weights);
            elecstate::calEBand(this->pelec->ekb, this->pelec->wg, this->pelec->f_en);
        }

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
