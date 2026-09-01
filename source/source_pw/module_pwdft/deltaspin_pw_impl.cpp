#include "source_base/matrix.h"
#include "source_base/parallel_reduce.h"
#include "source_base/tool_title.h"
#include "source_base/timer.h"
#include "source_base/kernels/math_kernel_op.h"
#include "source_pw/module_pwdft/onsite_proj.h"
#include "source_lcao/module_deltaspin/spin_constrain.h"
#include "source_lcao/module_deltaspin/mi_tools.h"
#include "source_io/module_parameter/parameter.h"
#include "source_hsolver/diago_iter_assist.h"
#include "source_hsolver/hsolver_pw.h"
#include "source_estate/elecstate_pw.h"
#include "source_estate/elecstate_tools.h"

namespace spinconstrain {

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
void SpinConstrain<std::complex<double>>::cal_mi_pw()
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
            const int spin_sign = (npol == 2) ? 1 : this->get_spin_sign(ik);
            accumulate_Mi_from_becp(becp, nkb, nbands, npol, spin_sign,
                &this->pelec->wg(ik, 0), &onsite_p->get_nh(0), this->Mi_);
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
            const int spin_sign = (npol == 2) ? 1 : this->get_spin_sign(ik);
            accumulate_Mi_from_becp(becp, nkb, nbands, npol, spin_sign,
                &this->pelec->wg(ik, 0), &onsite_p->get_nh(0), this->Mi_);
        }
    }
#endif
    // MPI reduction: sum Mi across all k-pool ranks
    Parallel_Reduce::reduce_double_allpool(PARAM.inp.kpar, GlobalV::NPROC_IN_POOL, &(this->Mi_[0][0]), 3 * this->Mi_.size());

    ModuleBase::timer::end("spinconstrain::SpinConstrain", "cal_mi_pw");
}

/**
 * @brief Compute DeltaSpin correction to the subspace Hamiltonian.
 *
 * @details Adds the constraint term to H in the projector subspace:
 *   H += becp^† * ps, where ps = delta_lambda * becp
 *
 * For non-collinear (npol=2), this implements the full 2x2 Pauli matrix:
 *   H_delta = | lambda_z     lambda_x + i*lambda_y |
 *             | lambda_x - i*lambda_y   -lambda_z  |
 *
 * For collinear (npol=1), only the diagonal z-component with spin_sign:
 *   H_delta = lambda_z * spin_sign
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
void SpinConstrain<std::complex<double>>::calculate_delta_hcc(std::complex<double>* h_tmp, const std::complex<double>* becp_k, const ModuleBase::Vector3<double>* delta_lambda, const int nbands, const int nkb, const int* nh_iat, const int ik, bool full_update)
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
        // ps = lambda_z * spin_sign * becp
        // spin_sign = +1 for spin-up k-points, -1 for spin-down
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
void SpinConstrain<std::complex<double>>::update_psi_charge_pw_cpu(const ModuleBase::Vector3<double>* delta_lambda, bool pw_solve, bool full_update)
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
            PARAM.inp.nbands,
            PARAM.inp.diago_smooth_ethr,
            PARAM.inp.pw_diag_ndim,
            PARAM.inp.diag_subspace,
            PARAM.inp.nb2d,
            PARAM.inp.use_k_continuity);

        hsolver_pw_obj.solve(hamilt_t,
                             psi_t[0],
                             this->pelec,
                             this->pelec->ekb.c,
                             GlobalV::RANK_IN_POOL,
                             GlobalV::NPROC_IN_POOL,
                             GlobalV::ofs_running,
                             false,
                             this->tpiba,
                             this->get_nat());
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
                                     PARAM.inp.nbands,
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
void SpinConstrain<std::complex<double>>::update_psi_charge_pw_gpu(const ModuleBase::Vector3<double>* delta_lambda, bool pw_solve, bool full_update)
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

    base_device::memory::delete_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(h_tmp);
    base_device::memory::delete_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(s_tmp);

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
            PARAM.inp.nbands,
            PARAM.inp.diago_smooth_ethr,
            PARAM.inp.pw_diag_ndim,
            PARAM.inp.diag_subspace,
            PARAM.inp.nb2d,
            PARAM.inp.use_k_continuity);

        hsolver_pw_obj.solve(hamilt_t,
                             psi_t[0],
                             this->pelec,
                             this->pelec->ekb.c,
                             GlobalV::RANK_IN_POOL,
                             GlobalV::NPROC_IN_POOL,
                             GlobalV::ofs_running,
                             false,
                             this->tpiba,
                             this->get_nat());
    }
    else
    {
        elecstate::calculate_weights(this->pelec->ekb,
                                     this->pelec->wg,
                                     this->pelec->klist,
                                     this->pelec->eferm,
                                     this->pelec->f_en,
                                     this->pelec->nelec_spin,
                                     PARAM.inp.nbands,
                                     this->pelec->skip_weights);
        elecstate::calEBand(this->pelec->ekb, this->pelec->wg, this->pelec->f_en);
        reinterpret_cast<elecstate::ElecStatePW<std::complex<double>, base_device::DEVICE_GPU>*>(this->pelec)->psiToRho(*psi_t);
    }
    ModuleBase::timer::end("spinconstrain::SpinConstrain", "update_psi_charge_pw_gpu");
}
#endif

} // namespace spinconstrain
