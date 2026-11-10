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
#include "source_pw/module_proj/onsite_proj.h"
#include "source_pw/module_pwdft/kernels/onsite_op.h"
#if defined(__CUDA)
#include "source_hsolver/kernels/fixed_overlap.h"
#endif
#include "deltaspin_pw_cache.h"
#include "hpsi_transfer.h"
#include "mi_tools.h"
#include "source_io/module_parameter/parameter.h"
#include "source_hsolver/diago_iter_assist.h"
#include "source_hamilt/hamilt_hs_adapter.h"
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
                &pelec->wg(ik, 0), &onsite_p->get_nh(0), state.get_mi());
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
        const int atoms = static_cast<int>(state.get_mi().size());
        const int* nh_iat = &onsite_p->get_nh(0);
        std::vector<int> offsets(atoms + 1, 0);
        for (int iat = 0; iat < atoms; ++iat)
        {
            offsets[iat + 1] = offsets[iat] + nh_iat[iat];
        }
        int* d_offsets = nullptr;
        double* d_mi = nullptr;
        double* d_weights = nullptr;
        base_device::memory::resize_memory_op<int, base_device::DEVICE_GPU>()(d_offsets, atoms + 1);
        base_device::memory::resize_memory_op<double, base_device::DEVICE_GPU>()(d_mi, 3 * atoms);
        base_device::memory::resize_memory_op<double, base_device::DEVICE_GPU>()(d_weights, nbands);
        base_device::memory::synchronize_memory_op<int, base_device::DEVICE_GPU, base_device::DEVICE_CPU>()(
            d_offsets, offsets.data(), offsets.size());
        base_device::memory::set_memory_op<double, base_device::DEVICE_GPU>()(d_mi, 0, 3 * atoms);
        for (int ik = 0; ik < nks; ik++)
        {
            psi_t->fix_k(ik);
            psi_pointer = psi_t->get_pointer();
            onsite_p->tabulate_atomic(ik);
            onsite_p->overlap_proj_psi(nbands * npol, psi_pointer);
            const std::complex<double>* becp = onsite_p->get_becp();
            const int nkb = onsite_p->get_tot_nproj();
            const int spin_sign = (npol == 2) ? 1 : spin_sign_at(state, pelec, ik);
            base_device::memory::synchronize_memory_op<double, base_device::DEVICE_GPU, base_device::DEVICE_CPU>()(
                d_weights, &pelec->wg(ik, 0), nbands);
            hamilt::accumulate_projector_moments_gpu(nbands, npol, atoms, nkb, spin_sign,
                                                     d_offsets, d_weights, becp, d_mi);
        }
        base_device::memory::synchronize_memory_op<double, base_device::DEVICE_CPU, base_device::DEVICE_GPU>()(
            &(state.get_mi()[0][0]), d_mi, 3 * atoms);
        base_device::memory::delete_memory_op<int, base_device::DEVICE_GPU>()(d_offsets);
        base_device::memory::delete_memory_op<double, base_device::DEVICE_GPU>()(d_mi);
        base_device::memory::delete_memory_op<double, base_device::DEVICE_GPU>()(d_weights);
    }
#endif
    // MPI reduction: sum Mi across all k-pool ranks
    Parallel_Reduce::reduce_double_allpool(PARAM.inp.kpar, GlobalV::NPROC_IN_POOL, &(state.get_mi()[0][0]), 3 * state.get_mi().size());

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

#if defined(__CUDA)
    if (PARAM.inp.device == "gpu")
    {
        using Complex = std::complex<double>;
        const int rows = nkb * state.npol_;
        const int atoms = static_cast<int>(state.get_mi().size());
        // Each band shares the same three constraint coefficients per projector.
        std::vector<double> coefficients(3 * nkb, 0.0);
        int offset = 0;
        for (int iat = 0; iat < atoms; ++iat)
        {
            for (int ip = 0; ip < nh_iat[iat]; ++ip)
            {
                coefficients[3 * (offset + ip)] = effective_lambda[iat].x;
                coefficients[3 * (offset + ip) + 1] = effective_lambda[iat].y;
                coefficients[3 * (offset + ip) + 2] = effective_lambda[iat].z;
            }
            offset += nh_iat[iat];
        }
        const int sign = state.npol_ == 2 ? 1 : spin_sign_at(state, pelec, ik);
        double* d_lambda = nullptr;
        Complex* d_ps = nullptr;
        base_device::memory::resize_memory_op<double, base_device::DEVICE_GPU>()(d_lambda, 3 * nkb);
        base_device::memory::resize_memory_op<Complex, base_device::DEVICE_GPU>()(d_ps, rows * nbands);
        base_device::memory::synchronize_memory_op<double, base_device::DEVICE_GPU, base_device::DEVICE_CPU>()(
            d_lambda, coefficients.data(), coefficients.size());
        hamilt::spin_constraint_projectors_gpu(nbands, state.npol_, nkb, sign, d_lambda, becp_k, d_ps);
        char transa = 'C';
        char transb = 'N';
        ModuleBase::gemm_op<Complex, base_device::DEVICE_GPU>()(transa,
                                                               transb,
                                                               nbands,
                                                               nbands,
                                                               rows,
                                                               &ModuleBase::ONE,
                                                               becp_k,
                                                               rows,
                                                               d_ps,
                                                               rows,
                                                               &ModuleBase::ONE,
                                                               h_tmp,
                                                               nbands);
        base_device::memory::delete_memory_op<double, base_device::DEVICE_GPU>()(d_lambda);
        base_device::memory::delete_memory_op<Complex, base_device::DEVICE_GPU>()(d_ps);
        ModuleBase::timer::end("spinconstrain::SpinConstrain", "calculate_delta_hcc");
        return;
    }
#endif

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
        // For each atom, construct lambda dot sigma:
        //   | lambda_z                  lambda_x - i*lambda_y |
        //   | lambda_x + i*lambda_y    -lambda_z              |
        // The coefficient array uses {up-up, down-up, up-down, down-down}.
        // Then: ps_up = coeff0 * becp_up + coeff2 * becp_dn
        //        ps_dn = coeff1 * becp_up + coeff3 * becp_dn
        for (size_t iat = 0; iat < state.get_mi().size(); iat++)
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
        for (size_t iat = 0; iat < state.get_mi().size(); iat++)
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

#if defined(__CUDA)
/// Add the DeltaSpin constraint action vkb * (delta_lambda * becp) to output.
/// output is (npwx*npol, nbands) with the solver's leading dimension and
/// becp is the cached projector coefficient matrix for the current k point.
void apply_constraint_action_gpu(ScState& state,
                                 elecstate::ElecState* pelec,
                                 const std::complex<double>* becp,
                                 const ModuleBase::Vector3<double>* delta_lambda,
                                 std::complex<double>* output,
                                 const int leading,
                                 const int npwx,
                                 const int ik,
                                 const int nbands,
                                 const int nkb,
                                 const int* nh_iat)
{
    using Complex = std::complex<double>;
    const int npol = state.npol_;
    const int rows = nkb * npol;
    const int atoms = static_cast<int>(state.get_mi().size());
    std::vector<double> coefficients(3 * nkb, 0.0);
    int offset = 0;
    for (int iat = 0; iat < atoms; ++iat)
    {
        for (int ip = 0; ip < nh_iat[iat]; ++ip)
        {
            coefficients[3 * (offset + ip)] = delta_lambda[iat].x;
            coefficients[3 * (offset + ip) + 1] = delta_lambda[iat].y;
            coefficients[3 * (offset + ip) + 2] = delta_lambda[iat].z;
        }
        offset += nh_iat[iat];
    }
    const int sign = npol == 2 ? 1 : spin_sign_at(state, pelec, ik);
    double* d_lambda = nullptr;
    Complex* d_ps = nullptr;
    base_device::memory::resize_memory_op<double, base_device::DEVICE_GPU>()(d_lambda, 3 * nkb);
    base_device::memory::resize_memory_op<Complex, base_device::DEVICE_GPU>()(d_ps, rows * nbands);
    base_device::memory::synchronize_memory_op<double, base_device::DEVICE_GPU, base_device::DEVICE_CPU>()(
        d_lambda, coefficients.data(), coefficients.size());
    hamilt::spin_constraint_projectors_gpu(nbands, npol, nkb, sign, d_lambda, becp, d_ps);
    auto* onsite = projectors::OnsiteProjector<double, base_device::DEVICE_GPU>::get_instance();
    const Complex* vkb = onsite->get_fs_tools()->get_vkb();
    for (int spin = 0; spin < npol; ++spin)
    {
        ModuleBase::gemm_op<Complex, base_device::DEVICE_GPU>()(
            'N', 'N', npwx, nbands, nkb,
            &ModuleBase::ONE, vkb, npwx, d_ps + spin * nkb, rows,
            &ModuleBase::ONE, output + spin * npwx, leading);
    }
    base_device::memory::delete_memory_op<double, base_device::DEVICE_GPU>()(d_lambda);
    base_device::memory::delete_memory_op<Complex, base_device::DEVICE_GPU>()(d_ps);
}
#endif

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
        lambda_for_hcc = state.get_sc_lambda().data();
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

        hamilt::HamiltHSOperator<std::complex<double>, base_device::DEVICE_CPU> op(hamilt_t, pw_wfc);
        hsolver_pw_obj.solve(op,
                             psi_t[0],
                             pelec,
                             pelec->ekb.c,
                             GlobalV::RANK_IN_POOL,
                             GlobalV::NPROC_IN_POOL,
                             GlobalV::ofs_running,
                             false);
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

    assert(cache.allocated());

    const ModuleBase::Vector3<double>* lambda_for_hcc = delta_lambda;
    if (full_update)
    {
        lambda_for_hcc = state.get_sc_lambda().data();
    }

    // Davidson diagonalizes the complete initial subspace itself. An invertible
    // rotation before that solve leaves the span unchanged, so skip the explicit
    // subspace update. Density-only updates and other solvers still need it.
    const bool rotate_basis = !pw_solve || PARAM.inp.ks_solver != "dav_subspace";
    std::complex<double>* h_tmp = nullptr;
    std::complex<double>* s_tmp = nullptr;
    if (rotate_basis)
    {
        base_device::memory::resize_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(h_tmp, nbands * nbands);
        base_device::memory::resize_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(s_tmp, nbands * nbands);
    }

    // STAGE 1: Subspace diagonalization for each k-point (GPU)
    for (int ik = 0; rotate_basis && ik < nk; ++ik)
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

    if (rotate_basis)
    {
        base_device::memory::delete_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(h_tmp);
        base_device::memory::delete_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(s_tmp);
    }

#if defined(__CUDA)
    // The response loop is complete; wavefunction updates use the original
    // H/S/projector matrices and have no use for the orthonormal copies.
    delete static_cast<hsolver::FixedOverlapSubspace*>(cache.fixed_overlap);
    cache.fixed_overlap = nullptr;
    if (cache.hpsi_host != nullptr)
    {
        CHECK_CUDA(cudaFreeHost(cache.hpsi_host));
        cache.hpsi_host = nullptr;
        cache.hpsi_rows = cache.hpsi_bands = cache.hpsi_nk = 0;
    }
#endif

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

        hamilt::HamiltHSOperator<std::complex<double>, base_device::DEVICE_GPU> op(hamilt_t, pw_wfc);
        bool use_action_cache = cache.hpsi_host != nullptr && cache.hpsi_nk == nk
                                && cache.hpsi_bands == nbands;
        if (use_action_cache)
        {
            const int leading = psi_t->get_current_ngk();
            use_action_cache = cache.hpsi_rows == leading;
        }
        if (use_action_cache)
        {
            std::vector<ModuleBase::Vector3<double>> delta_lambda(nh_iat != nullptr ? state.get_nat() : 0);
            for (int iat = 0; iat < state.get_nat(); ++iat)
            {
                delta_lambda[iat] = state.get_sc_lambda()[iat] - cache.lambda_in_sub()[iat];
            }
            const int npwx = onsite_p->get_npwx();
            const int action_bands = nbands;
            const int action_rows = cache.hpsi_rows;
            spinconstrain::HPsiTransfer transfer(cache.hpsi_host, cache.hpsi_nk, action_bands, action_rows);
            typename hsolver::Diago_DavSubspace<std::complex<double>, base_device::DEVICE_GPU>::InitialHPsiFunc
                initial_hpsi = [&](int ik, std::complex<double>* output, int output_leading) {
                    transfer.restore(ik, output, output_leading);
                    // Prepare the next projector block while the pinned upload
                    // is in flight; the default stream orders the following GEMM.
                    onsite_p->tabulate_atomic(ik);
                    apply_constraint_action_gpu(state, pelec, cache.becp_k(ik, nbands * nkb * npol),
                                                delta_lambda.data(), output, output_leading, npwx, ik,
                                                nbands, nkb, nh_iat);
                };
            typename hsolver::Diago_DavSubspace<std::complex<double>,
                                                base_device::DEVICE_GPU>::BasisObserver observer;
            const int becp_rows = nkb * npol;
            observer.rotate = [&](int observer_ik, int count, const std::complex<double>* coefficients,
                                  int coefficient_leading, bool /*final*/) {
                if (observer_ik < 0 || observer_ik >= nk || count <= 0)
                {
                    return;
                }
                std::complex<double>* becp = cache.becp_k(observer_ik, nbands * nkb * npol);
                std::complex<double>* rotated = nullptr;
                base_device::memory::resize_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(
                    rotated, size_t(becp_rows) * nbands);
                ModuleBase::gemm_op<std::complex<double>, base_device::DEVICE_GPU>()(
                    'N', 'N', becp_rows, nbands, count,
                    &ModuleBase::ONE, becp, becp_rows, coefficients, coefficient_leading,
                    &ModuleBase::ZERO, rotated, becp_rows);
                base_device::memory::synchronize_memory_op<std::complex<double>,
                                                           base_device::DEVICE_GPU,
                                                           base_device::DEVICE_GPU>()(
                    becp, rotated, size_t(becp_rows) * nbands);
                base_device::memory::delete_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(rotated);
            };
            hsolver_pw_obj.solve(op,
                                 psi_t[0],
                                 pelec,
                                 pelec->ekb.c,
                                 GlobalV::RANK_IN_POOL,
                                 GlobalV::NPROC_IN_POOL,
                                 GlobalV::ofs_running,
                                 false,
                                 initial_hpsi,
                                 &observer);
        }
        else
        {
            hsolver_pw_obj.solve(op,
                                 psi_t[0],
                                 pelec,
                                 pelec->ekb.c,
                                 GlobalV::RANK_IN_POOL,
                                 GlobalV::NPROC_IN_POOL,
                                 GlobalV::ofs_running,
                                 false);
        }
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
