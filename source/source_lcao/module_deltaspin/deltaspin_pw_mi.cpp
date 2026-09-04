/**
 * @file deltaspin_pw_mi.cpp
 * @brief PW-basis DeltaSpin computation path (free functions).
 *
 * @details Implementation moved from source/source_pw/module_pwdft/deltaspin_pw_impl.cpp.
 * The former SpinConstrain<std::complex<double>> member functions
 * (cal_mi_pw / calculate_delta_hcc / update_psi_charge_pw_cpu / update_psi_charge_pw_gpu)
 * are now spinconstrain::pw free functions; their dependencies (state, cache, psi,
 * hamiltonian, electronic state, PW basis) are passed explicitly.
 *
 * The PW path is always instantiated on complex wavefunctions, so these functions
 * are not templated on TK. This file must NOT be wrapped in #ifdef __LCAO so that
 * PW DeltaSpin also compiles when ENABLE_LCAO=off (matching the previous behaviour
 * in module_pwdft, which was compiled unconditionally).
 */
#include "deltaspin_pw_mi.h"

#include <cassert>
#include <cstring>
#include <vector>

#include "source_base/matrix.h"
#include "source_base/parallel_reduce.h"
#include "source_base/tool_title.h"
#include "source_base/timer.h"
#include "source_base/kernels/math_kernel_op.h"
#include "source_base/module_device/device.h"
#include "source_pw/module_pwdft/onsite_proj.h"
#include "deltaspin_pw_cache.h"
#include "mi_tools.h"
#include "source_io/module_parameter/parameter.h"
#include "source_hsolver/diago_iter_assist.h"
#include "source_hsolver/hsolver_pw.h"
#include "source_estate/elecstate.h"
#include "source_estate/elecstate_pw.h"
#include "source_estate/elecstate_tools.h"
#include "source_psi/psi.h"

namespace spinconstrain
{
namespace pw
{

namespace
{
/// Collinear spin sign for k-point ik: +1 for spin-up, -1 for spin-down.
/// Returns 1 for non-collinear (npol == 2).
inline int spin_sign_at(const ScState& state, const elecstate::ElecState* pelec, int ik)
{
    if (state.get_npol() == 2)
    {
        return 1;
    }
    return (pelec->klist->isk[ik] == 0) ? 1 : -1;
}
} // namespace

void cal_mi_pw(ScState& state,
               void* psi,
               elecstate::ElecState* pelec)
{
    ModuleBase::TITLE("module_deltaspin", "cal_mi_pw");
    ModuleBase::timer::start("spinconstrain::SpinConstrain", "cal_mi_pw");

    state.zero_Mi();
    if (PARAM.inp.device == "cpu")
    {
        auto* onsite_p = projectors::OnsiteProjector<double, base_device::DEVICE_CPU>::get_instance();
        // Loop over k-points to calculate Mi of sum_{k,i,l,m}<Psi_{k,i}|alpha_{l,m}><alpha_{l,m}|Psi_{k,i}>
        std::complex<double>* psi_pointer = nullptr;
        psi::Psi<std::complex<double>, base_device::DEVICE_CPU>* psi_t = static_cast<psi::Psi<std::complex<double>, base_device::DEVICE_CPU>*>(psi);
        const int nbands = psi_t->get_nbands();
        const int nks = psi_t->get_nk();
        const int npol = psi_t->get_npol();
        for (int ik = 0; ik < nks; ik++)
        {
            psi_t->fix_k(ik);
            psi_pointer = psi_t->get_pointer();
            onsite_p->tabulate_atomic(ik); // Set up atomic projectors for this k-point
            onsite_p->overlap_proj_psi(nbands * npol, psi_pointer); // Compute becp = <alpha|psi>
            const std::complex<double>* becp = onsite_p->get_h_becp();
            int nkb = onsite_p->get_tot_nproj();
            const int spin_sign = (npol == 2) ? 1 : spin_sign_at(state, pelec, ik);
            accumulate_Mi_from_becp(becp, nkb, nbands, npol, spin_sign,
                &pelec->wg(ik, 0), &onsite_p->get_nh(0), state.Mi_);
        }
    }
#if ((defined __CUDA) || (defined __ROCM))
    else
    {
        auto* onsite_p = projectors::OnsiteProjector<double, base_device::DEVICE_GPU>::get_instance();
        std::complex<double>* psi_pointer = nullptr;
        psi::Psi<std::complex<double>, base_device::DEVICE_GPU>* psi_t = static_cast<psi::Psi<std::complex<double>, base_device::DEVICE_GPU>*>(psi);
        const int nbands = psi_t->get_nbands();
        const int nks = psi_t->get_nk();
        const int npol = psi_t->get_npol();
        for (int ik = 0; ik < nks; ik++)
        {
            psi_t->fix_k(ik);
            psi_pointer = psi_t->get_pointer();
            onsite_p->tabulate_atomic(ik);
            onsite_p->overlap_proj_psi(nbands * npol, psi_pointer);
            const std::complex<double>* becp = onsite_p->get_h_becp();
            int nkb = onsite_p->get_size_becp() / nbands / npol;
            const int spin_sign = (npol == 2) ? 1 : spin_sign_at(state, pelec, ik);
            accumulate_Mi_from_becp(becp, nkb, nbands, npol, spin_sign,
                &pelec->wg(ik, 0), &onsite_p->get_nh(0), state.Mi_);
        }
    }
#endif
    // MPI reduction: sum Mi across all k-pool ranks
    Parallel_Reduce::reduce_double_allpool(PARAM.inp.kpar, GlobalV::NPROC_IN_POOL, &(state.Mi_[0][0]), 3 * state.Mi_.size());

    ModuleBase::timer::end("spinconstrain::SpinConstrain", "cal_mi_pw");
}

void calculate_delta_hcc(ScState& state,
                         const SubspaceCache& cache,
                         elecstate::ElecState* pelec,
                         std::complex<double>* h_tmp,
                         const std::complex<double>* becp_k,
                         const ModuleBase::Vector3<double>* delta_lambda,
                         const int nbands,
                         const int nkb,
                         const int* nh_iat,
                         const int ik,
                         const bool full_update)
{
    ModuleBase::TITLE("spinconstrain::SpinConstrain", "calculate_delta_hcc");
    ModuleBase::timer::start("spinconstrain::SpinConstrain", "calculate_delta_hcc");

    // If full_update, compute actual delta = lambda_current - lambda_at_save
    // This applies only the CHANGE in lambda, not the full lambda value
    std::vector<ModuleBase::Vector3<double>> actual_delta;
    const ModuleBase::Vector3<double>* effective_lambda = delta_lambda;
    if (full_update)
    {
        int nat = state.get_nat();
        actual_delta.resize(nat);
        for (int iat = 0; iat < nat; iat++)
        {
            actual_delta[iat] = delta_lambda[iat] - cache.lambda_in_sub()[iat];
        }
        effective_lambda = actual_delta.data();
    }

    int sum = 0; // Running sum of projectors across atoms
    int size_ps = nkb * state.npol_ * nbands; // Total size of ps array
    std::complex<double>* becp_cpu = nullptr;

    // Handle GPU/CPU memory for becp
    if (PARAM.inp.device == "gpu")
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
    if (state.npol_ == 2)
    {
        // =============================================================
        // nspin=4 (non-collinear): full Pauli matrix treatment
        // =============================================================
        // For each atom, construct 2x2 coefficients:
        //   | lambda_z      lambda_x + i*lambda_y |
        //   | lambda_x - i*lambda_y   -lambda_z   |
        // Then: ps_up = coeff0 * becp_up + coeff2 * becp_dn
        //        ps_dn = coeff1 * becp_up + coeff3 * becp_dn
        for (size_t iat = 0; iat < state.Mi_.size(); iat++)
        {
            const int nproj = nh_iat[iat];
            const std::complex<double> coefficients0(effective_lambda[iat][2], 0.0);
            const std::complex<double> coefficients1(effective_lambda[iat][0], effective_lambda[iat][1]);
            const std::complex<double> coefficients2(effective_lambda[iat][0], -1 * effective_lambda[iat][1]);
            const std::complex<double> coefficients3(-1 * effective_lambda[iat][2], 0.0);
            for (int ib = 0; ib < nbands * state.npol_; ib += state.npol_)
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
    else if (state.npol_ == 1)
    {
        // =============================================================
        // nspin=2 (collinear): only z-component with spin_sign
        // =============================================================
        // ps = lambda_z * spin_sign * becp
        // spin_sign = +1 for spin-up k-points, -1 for spin-down
        const int spin_sign = spin_sign_at(state, pelec, ik);
        for (size_t iat = 0; iat < state.Mi_.size(); iat++)
        {
            const int nproj = nh_iat[iat];
            double coefficients0 = effective_lambda[iat][2] * spin_sign;
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
    if (PARAM.inp.device == "gpu")
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
    // H += becp^dagger * ps (GEMM: C = alpha * A^dagger * B + beta * C)
    // A = becp_k (npm x nbands), B = ps (npm x nbands), C = h_tmp (nbands x nbands)
    // =============================================================
    char transa = 'C'; // Conjugate transpose of becp
    char transb = 'N'; // Normal ps
    const int npm = nkb * state.npol_;
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

void update_psi_charge_pw_cpu(ScState& state,
                              SubspaceCache& cache,
                              void* psi,
                              void* p_hamilt,
                              elecstate::ElecState* pelec,
                              ModulePW::PW_Basis_K* pw_wfc,
                              const ModuleBase::Vector3<double>* delta_lambda,
                              bool pw_solve,
                              bool full_update)
{
    ModuleBase::TITLE("spinconstrain::SpinConstrain", "update_psi_charge_pw_cpu");
    ModuleBase::timer::start("spinconstrain::SpinConstrain", "update_psi_charge_pw_cpu");

    psi::Psi<std::complex<double>>* psi_t = static_cast<psi::Psi<std::complex<double>>*>(psi);
    hamilt::Hamilt<std::complex<double>, base_device::DEVICE_CPU>* hamilt_t = static_cast<hamilt::Hamilt<std::complex<double>, base_device::DEVICE_CPU>*>(p_hamilt);
    auto* onsite_p = projectors::OnsiteProjector<double, base_device::DEVICE_CPU>::get_instance();

    int nbands = psi_t->get_nbands();
    int npol = psi_t->get_npol();
    int nkb = onsite_p->get_tot_nproj();
    int nk = psi_t->get_nk();
    int size_becp = nbands * nkb * npol;
    const int* nh_iat = &onsite_p->get_nh(0);

    std::vector<std::complex<double>> h_tmp(nbands * nbands), s_tmp(nbands * nbands);

    // CRITICAL: subspace data must have been saved by cal_mw_from_lambda()
    assert(cache.allocated());

    // Determine which lambda to use for H correction
    const ModuleBase::Vector3<double>* lambda_for_hcc = delta_lambda;
    if (full_update)
    {
        lambda_for_hcc = state.lambda_.data();
    }

    // =============================================================
    // STAGE 1: Subspace diagonalization for each k-point
    // =============================================================
    for (int ik = 0; ik < nk; ++ik)
    {
        std::complex<double>* h_k = cache.h_k(ik, nbands);
        std::complex<double>* s_k = cache.s_k(ik, nbands);
        std::complex<double>* becp_k = cache.becp_k(ik, size_becp);

        psi_t->fix_k(ik);

        // Copy saved subspace matrices to temp
        memcpy(h_tmp.data(), h_k, sizeof(std::complex<double>) * nbands * nbands);
        memcpy(s_tmp.data(), s_k, sizeof(std::complex<double>) * nbands * nbands);

        // Apply DeltaSpin correction: H += becp^dagger * lambda * becp
        calculate_delta_hcc(state, cache, pelec, h_tmp.data(), becp_k, lambda_for_hcc, nbands, nkb, nh_iat, ik, full_update);

        // Diagonalize in subspace to update wavefunction coefficients and eigenvalues
        hsolver::DiagoIterAssist<std::complex<double>>::diag_subspace_psi(h_tmp.data(),
                                                                        s_tmp.data(),
                                                                        nbands,
                                                                        psi_t[0],
                                                                        &pelec->ekb(ik, 0));
    }

    // Free saved subspace data (allocated in cal_mw_from_lambda)
    cache.release_cpu();

    // =============================================================
    // STAGE 2: Full-space update
    // =============================================================
    if (pw_solve)
    {
        // Full PW diagonalization: subspace rotation provides a good initial guess,
        // then HSolverPW iteratively refines psi in the full plane-wave space and calls psiToRho.
        hsolver::HSolverPW<std::complex<double>, base_device::DEVICE_CPU> hsolver_pw_obj(
            pw_wfc,
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

        hsolver_pw_obj.solve(hamilt_t, psi_t[0], pelec, pelec->ekb.c,
            GlobalV::RANK_IN_POOL, GlobalV::NPROC_IN_POOL, false, state.tpiba, state.get_nat());
    }
    else
    {
        // No full solver: update weights from new eigenvalues, then build rho from current psi
        elecstate::calculate_weights(pelec->ekb,
                                     pelec->wg,
                                     pelec->klist,
                                     pelec->eferm,
                                     pelec->f_en,
                                     pelec->nelec_spin,
                                     PARAM.inp.nbands,
                                     pelec->skip_weights);
        elecstate::calEBand(pelec->ekb, pelec->wg, pelec->f_en);
        reinterpret_cast<elecstate::ElecStatePW<std::complex<double>, base_device::DEVICE_CPU>*>(pelec)->psiToRho(*psi_t);
    }
    ModuleBase::timer::end("spinconstrain::SpinConstrain", "update_psi_charge_pw_cpu");
}

#if ((defined __CUDA) || (defined __ROCM))
void update_psi_charge_pw_gpu(ScState& state,
                              SubspaceCache& cache,
                              void* psi,
                              void* p_hamilt,
                              elecstate::ElecState* pelec,
                              ModulePW::PW_Basis_K* pw_wfc,
                              const ModuleBase::Vector3<double>* delta_lambda,
                              bool pw_solve,
                              bool full_update)
{
    ModuleBase::TITLE("spinconstrain::SpinConstrain", "update_psi_charge_pw_gpu");
    ModuleBase::timer::start("spinconstrain::SpinConstrain", "update_psi_charge_pw_gpu");

    psi::Psi<std::complex<double>, base_device::DEVICE_GPU>* psi_t = static_cast<psi::Psi<std::complex<double>, base_device::DEVICE_GPU>*>(psi);
    hamilt::Hamilt<std::complex<double>, base_device::DEVICE_GPU>* hamilt_t = static_cast<hamilt::Hamilt<std::complex<double>, base_device::DEVICE_GPU>*>(p_hamilt);
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

    assert(cache.allocated());

    const ModuleBase::Vector3<double>* lambda_for_hcc = delta_lambda;
    if (full_update)
    {
        lambda_for_hcc = state.lambda_.data();
    }

    // STAGE 1: Subspace diagonalization for each k-point (GPU)
    for (int ik = 0; ik < nk; ++ik)
    {
        std::complex<double>* h_k = cache.h_k(ik, nbands);
        std::complex<double>* s_k = cache.s_k(ik, nbands);
        std::complex<double>* becp_k = cache.becp_k(ik, size_becp);

        psi_t->fix_k(ik);

        base_device::memory::synchronize_memory_op<std::complex<double>, base_device::DEVICE_GPU, base_device::DEVICE_GPU>()(h_tmp, h_k, nbands * nbands);
        base_device::memory::synchronize_memory_op<std::complex<double>, base_device::DEVICE_GPU, base_device::DEVICE_GPU>()(s_tmp, s_k, nbands * nbands);

        calculate_delta_hcc(state, cache, pelec, h_tmp, becp_k, lambda_for_hcc, nbands, nkb, nh_iat, ik, full_update);

        hsolver::DiagoIterAssist<std::complex<double>, base_device::DEVICE_GPU>::diag_subspace_psi(h_tmp,
                                                                                s_tmp,
                                                                                nbands,
                                                                                psi_t[0],
                                                                                &pelec->ekb(ik, 0));
    }

    base_device::memory::delete_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(h_tmp);
    base_device::memory::delete_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(s_tmp);

    // Free GPU memory for saved subspace data
    cache.release_gpu();

    // STAGE 2: Full-space update (GPU)
    if (pw_solve)
    {
        hsolver::HSolverPW<std::complex<double>, base_device::DEVICE_GPU> hsolver_pw_obj(
            pw_wfc,
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

        hsolver_pw_obj.solve(hamilt_t, psi_t[0], pelec, pelec->ekb.c,
            GlobalV::RANK_IN_POOL, GlobalV::NPROC_IN_POOL, false, state.tpiba, state.get_nat());
    }
    else
    {
        elecstate::calculate_weights(pelec->ekb,
                                     pelec->wg,
                                     pelec->klist,
                                     pelec->eferm,
                                     pelec->f_en,
                                     pelec->nelec_spin,
                                     PARAM.inp.nbands,
                                     pelec->skip_weights);
        elecstate::calEBand(pelec->ekb, pelec->wg, pelec->f_en);
        reinterpret_cast<elecstate::ElecStatePW<std::complex<double>, base_device::DEVICE_GPU>*>(pelec)->psiToRho(*psi_t);
    }
    ModuleBase::timer::end("spinconstrain::SpinConstrain", "update_psi_charge_pw_gpu");
}
#endif // __CUDA || __ROCM

} // namespace pw
} // namespace spinconstrain
