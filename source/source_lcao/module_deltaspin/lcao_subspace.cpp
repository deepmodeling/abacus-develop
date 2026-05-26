/**
 * @file lcao_subspace.cpp
 * @brief LCAO subspace diagonalization for DeltaSpin lambda loop optimization.
 *
 * @par Purpose
 * Implements the subspace diagonalization strategy for LCAO basis in DeltaSpin,
 * mirroring the PW path. The key idea:
 *
 *   First call (i_step = -1): Full diagonalization -> cache H0_sub, S_sub, P_I_sub
 *   Subsequent calls (i_step >= 0): H_sub = H0_sub + Sum lambda_I P_I_sub -> subspace diag
 *
 * @par Parallel Strategy
 * - H_sub = Ct.H.C and S_sub = Ct.S.C are computed via ScaLAPACK pzgemm
 *   on the 2D-block distributed H(k), S(k), C(k).
 *   The result (nbands x nbands) is gathered to all processes in the pool.
 * - P_I_sub is computed by cal_PI_sub() which already does MPI_Allreduce.
 * - Subspace diagonalization (nbands x nbands) uses local LAPACK zhegvd.
 * - Psi rotation C_new = C_old . V uses ScaLAPACK pzgemm.
 *
 * @par Layout Conventions
 * - H_sub, S_sub, V: column-major nbands x nbands (LAPACK convention) after gathering
 * - desc_Eij local layout: nrow x ncol_bands (lld = nrow), NOT ncol_bands x ncol_bands
 * - C(k): 2D-block distributed (ScaLAPACK convention, ParaV->desc_wfc)
 * - P_I_sub[iat]: column-major nbands x nbands Hermitian
 */

#include "spin_constrain.h"
#include "source_base/timer.h"
#include "source_base/tool_title.h"
#include "source_base/global_variable.h"
#include "source_hsolver/diago_iter_assist.h"
#include "source_hsolver/kernels/hegvd_op.h"
#include "source_io/module_parameter/parameter.h"
#include "source_base/parallel_reduce.h"
#include "source_base/kernels/math_kernel_op.h"
#include "source_base/module_container/base/third_party/lapack.h"
#include "source_base/module_external/blas_connector.h"
#include "source_estate/elecstate_tools.h"
#include "source_estate/module_dm/cal_dm_psi.h"
#include "source_lcao/module_operator_lcao/dspin_lcao.h"

#ifdef __MPI
#include "source_base/module_external/scalapack_connector.h"
#include <mpi.h>
#endif

#include <complex>
#include <vector>
#include <cstring>
#include <cassert>
#include <algorithm>

namespace spinconstrain
{

template <>
void SpinConstrain<std::complex<double>>::free_lcao_subspace_cache()
{
    delete[] lcao_sub_h_save;
    delete[] lcao_sub_s_save;
    lcao_sub_h_save = nullptr;
    lcao_sub_s_save = nullptr;
    lcao_PI_sub_save_.clear();
    lcao_ekb_save_.clear();
    lcao_lambda_in_sub_.clear();
    lcao_subspace_initialized_ = false;
}

/**
 * @brief Gather a distributed nbands x nbands matrix to all processes in the pool.
 *
 * @details Each process holds its local portion of the ScaLAPACK-distributed matrix
 * (desc_Eij layout). The local storage has leading dimension = nrow (not ncol_bands).
 * Row indices map via local2global_row, column indices via local2global_col.
 *
 * @param local_data Local portion of the distributed matrix (lld = nrow)
 * @param ParaV Parallel orbitals distribution info
 * @param full_data Output: full nbands x nbands matrix on all processes (column-major)
 * @param nbands Number of bands
 */
static void gather_sub_matrix_to_all(
    const std::complex<double>* local_data,
    const Parallel_Orbitals* ParaV,
    std::complex<double>* full_data,
    int nbands)
{
#ifdef __MPI
    int nrow = ParaV->nrow;
    int ncol_bands = ParaV->ncol_bands;

    std::fill(full_data, full_data + nbands * nbands, std::complex<double>(0.0, 0.0));

    for (int j_local = 0; j_local < ncol_bands; j_local++)
    {
        int j_global = ParaV->local2global_col(j_local);
        if (j_global >= nbands) continue;
        for (int i_local = 0; i_local < nrow; i_local++)
        {
            int i_global = ParaV->local2global_row(i_local);
            if (i_global >= nbands) continue;
            full_data[i_global + j_global * nbands] = local_data[i_local + j_local * nrow];
        }
    }

    Parallel_Reduce::reduce_pool(full_data, nbands * nbands);
#else
    std::memcpy(full_data, local_data, sizeof(std::complex<double>) * nbands * nbands);
#endif
}

/**
 * @brief Scatter a full nbands x nbands matrix to local desc_Eij storage.
 *
 * @param full_data Full nbands x nbands matrix (column-major)
 * @param ParaV Parallel orbitals distribution info
 * @param local_data Output: local portion in desc_Eij layout (lld = nrow)
 * @param nbands Number of bands
 */
static void scatter_sub_matrix_to_local(
    const std::complex<double>* full_data,
    const Parallel_Orbitals* ParaV,
    std::complex<double>* local_data,
    int nbands)
{
#ifdef __MPI
    int nrow = ParaV->nrow;
    int ncol_bands = ParaV->ncol_bands;

    std::fill(local_data, local_data + nrow * ncol_bands, std::complex<double>(0.0, 0.0));

    for (int j_local = 0; j_local < ncol_bands; j_local++)
    {
        int j_global = ParaV->local2global_col(j_local);
        if (j_global >= nbands) continue;
        for (int i_local = 0; i_local < nrow; i_local++)
        {
            int i_global = ParaV->local2global_row(i_local);
            if (i_global >= nbands) continue;
            local_data[i_local + j_local * nrow] = full_data[i_global + j_global * nbands];
        }
    }
#else
    std::memcpy(local_data, full_data, sizeof(std::complex<double>) * nbands * nbands);
#endif
}

/**
 * @brief Compute H_sub = Ct.H.C and S_sub = Ct.S.C for LCAO subspace.
 */
template <>
void SpinConstrain<std::complex<double>>::calculate_lcao_sub_hs(
    void* hamilt,
    psi::Psi<std::complex<double>>& psi,
    const Parallel_Orbitals* ParaV,
    std::complex<double>* h_sub,
    std::complex<double>* s_sub,
    int ik, int nbands, int nlocal)
{
    ModuleBase::TITLE("SpinConstrain", "calculate_lcao_sub_hs");
    ModuleBase::timer::start("SpinConstrain", "calculate_lcao_sub_hs");

    hamilt::Hamilt<std::complex<double>>* hamilt_t
        = static_cast<hamilt::Hamilt<std::complex<double>>*>(hamilt);

    hamilt_t->updateHk(ik);

    hamilt::MatrixBlock<std::complex<double>> h_mat, s_mat;
    hamilt_t->matrix(h_mat, s_mat);

    psi.fix_k(ik);
    std::complex<double>* psi_ptr = psi.get_pointer();

    const std::complex<double> one = {1.0, 0.0};
    const std::complex<double> zero = {0.0, 0.0};

#ifdef __MPI
    int nloc_wfc = ParaV->nloc_wfc;
    std::vector<std::complex<double>> temp(nloc_wfc, zero);

    ScalapackConnector::gemm('N', 'N', nlocal, nbands, nlocal,
        one, h_mat.p, 1, 1, ParaV->desc,
        psi_ptr, 1, 1, ParaV->desc_wfc,
        zero, temp.data(), 1, 1, ParaV->desc_wfc);

    // desc_Eij has lld = nrow, so local buffer must be nrow * ncol_bands
    int lld_Eij = ParaV->nrow;
    int ncol_bands = ParaV->ncol_bands;
    int nloc_eij_local = lld_Eij * ncol_bands;

    std::vector<std::complex<double>> h_sub_local(nloc_eij_local, zero);

    ScalapackConnector::gemm('C', 'N', nbands, nbands, nlocal,
        one, psi_ptr, 1, 1, ParaV->desc_wfc,
        temp.data(), 1, 1, ParaV->desc_wfc,
        zero, h_sub_local.data(), 1, 1, ParaV->desc_Eij);

    gather_sub_matrix_to_all(h_sub_local.data(), ParaV, h_sub, nbands);

    // Repeat for S_sub
    std::vector<std::complex<double>> temp_s(nloc_wfc, zero);

    ScalapackConnector::gemm('N', 'N', nlocal, nbands, nlocal,
        one, s_mat.p, 1, 1, ParaV->desc,
        psi_ptr, 1, 1, ParaV->desc_wfc,
        zero, temp_s.data(), 1, 1, ParaV->desc_wfc);

    std::vector<std::complex<double>> s_sub_local(nloc_eij_local, zero);

    ScalapackConnector::gemm('C', 'N', nbands, nbands, nlocal,
        one, psi_ptr, 1, 1, ParaV->desc_wfc,
        temp_s.data(), 1, 1, ParaV->desc_wfc,
        zero, s_sub_local.data(), 1, 1, ParaV->desc_Eij);

    gather_sub_matrix_to_all(s_sub_local.data(), ParaV, s_sub, nbands);
#else
    // Serial case: use zgemm directly
    std::vector<std::complex<double>> temp(nlocal * nbands, zero);
    zgemm_("N", "N", &nlocal, &nbands, &nlocal,
        &one, h_mat.p, &nlocal, psi_ptr, &nlocal, &zero, temp.data(), &nlocal);

    zgemm_("C", "N", &nbands, &nbands, &nlocal,
        &one, psi_ptr, &nlocal, temp.data(), &nlocal, &zero, h_sub, &nbands);

    // S_sub
    std::vector<std::complex<double>> temp_s(nlocal * nbands, zero);
    zgemm_("N", "N", &nlocal, &nbands, &nlocal,
        &one, s_mat.p, &nlocal, psi_ptr, &nlocal, &zero, temp_s.data(), &nlocal);

    zgemm_("C", "N", &nbands, &nbands, &nlocal,
        &one, psi_ptr, &nlocal, temp_s.data(), &nlocal, &zero, s_sub, &nbands);
#endif

    ModuleBase::timer::end("SpinConstrain", "calculate_lcao_sub_hs");
}

template <>
void SpinConstrain<std::complex<double>>::calculate_delta_hcc_lcao(
    std::complex<double>* h_sub,
    const std::vector<std::vector<std::complex<double>>>& PI_sub,
    const ModuleBase::Vector3<double>* lambda,
    int nbands, int ik, bool full_update)
{
    ModuleBase::TITLE("SpinConstrain", "calculate_delta_hcc_lcao");
    ModuleBase::timer::start("SpinConstrain", "calculate_delta_hcc_lcao");

    const int nat = this->get_nat();
    std::vector<ModuleBase::Vector3<double>> actual_delta;
    const ModuleBase::Vector3<double>* effective_lambda = lambda;

    if (full_update)
    {
        actual_delta.resize(nat);
        for (int iat = 0; iat < nat; iat++)
        {
            actual_delta[iat] = lambda[iat] - this->lcao_lambda_in_sub_[iat];
        }
        effective_lambda = actual_delta.data();
    }

    if (this->npol_ == 2)
    {
        // nspin=4 (non-collinear)
        // Currently implemented: z-component (dominant term)
        // Full Pauli treatment TODO: requires 2nbands x 2nbands subspace structure
        for (int iat = 0; iat < nat; iat++)
        {
            if (PI_sub[iat].empty()) continue;

            const std::complex<double> coeff0(effective_lambda[iat][2], 0.0);
            const std::complex<double> coeff1(effective_lambda[iat][0], effective_lambda[iat][1]);
            const std::complex<double> coeff2(effective_lambda[iat][0], -effective_lambda[iat][1]);
            const std::complex<double> coeff3(-effective_lambda[iat][2], 0.0);

            const std::complex<double>* pi_ptr = PI_sub[iat].data();
            for (int ij = 0; ij < nbands * nbands; ij++)
            {
                h_sub[ij] += coeff0 * pi_ptr[ij];
            }
        }
    }
    else if (this->npol_ == 1)
    {
        // nspin=2 (collinear): H_sub += -spin_sign * lambda_z * P_I_sub
        // NOTE: The negative sign matches cal_coeff_lambda() in dspin_lcao.cpp:
        //   spin-up:   coeff = -delta_lambda_z
        //   spin-down: coeff = +delta_lambda_z
        int spin_sign = this->get_spin_sign(ik);

        for (int iat = 0; iat < nat; iat++)
        {
            if (PI_sub[iat].empty()) continue;

            const std::complex<double> coeff(-effective_lambda[iat][2] * spin_sign, 0.0);
            const std::complex<double>* pi_ptr = PI_sub[iat].data();
            for (int ij = 0; ij < nbands * nbands; ij++)
            {
                h_sub[ij] += coeff * pi_ptr[ij];
            }
        }
    }

    ModuleBase::timer::end("SpinConstrain", "calculate_delta_hcc_lcao");
}

template <>
void SpinConstrain<std::complex<double>>::rotate_psi_subspace_lcao(
    psi::Psi<std::complex<double>>& psi,
    const Parallel_Orbitals* ParaV,
    const std::vector<std::vector<std::complex<double>>>& vcc_all,
    int nbands, int nlocal, int nk)
{
    ModuleBase::TITLE("SpinConstrain", "rotate_psi_subspace_lcao");
    ModuleBase::timer::start("SpinConstrain", "rotate_psi_subspace_lcao");

    const std::complex<double> one = {1.0, 0.0};
    const std::complex<double> zero = {0.0, 0.0};

    for (int ik = 0; ik < nk; ik++)
    {
        psi.fix_k(ik);
        std::complex<double>* psi_ptr = psi.get_pointer();

#ifdef __MPI
        int nloc_wfc = ParaV->nloc_wfc;
        int lld_Eij = ParaV->nrow;
        int ncol_bands = ParaV->ncol_bands;
        int nloc_eij_local = lld_Eij * ncol_bands;

        // Scatter full V matrix (nbands x nbands) to local desc_Eij layout
        std::vector<std::complex<double>> v_local(nloc_eij_local, zero);
        scatter_sub_matrix_to_local(vcc_all[ik].data(), ParaV, v_local.data(), nbands);

        std::vector<std::complex<double>> temp(nloc_wfc, zero);

        // C_new = C_old . V (ScaLAPACK pzgemm)
        ScalapackConnector::gemm('N', 'N', nlocal, nbands, nbands,
            one, psi_ptr, 1, 1, ParaV->desc_wfc,
            v_local.data(), 1, 1, ParaV->desc_Eij,
            zero, temp.data(), 1, 1, ParaV->desc_wfc);

        std::memcpy(psi_ptr, temp.data(), sizeof(std::complex<double>) * nloc_wfc);
#else
        // Serial case
        std::vector<std::complex<double>> temp(nlocal * nbands);
        zgemm_("N", "N", &nlocal, &nbands, &nbands,
            &one, psi_ptr, &nlocal,
            vcc_all[ik].data(), &nbands,
            &zero, temp.data(), &nlocal);
        std::memcpy(psi_ptr, temp.data(), sizeof(std::complex<double>) * nlocal * nbands);
#endif
    }

    ModuleBase::timer::end("SpinConstrain", "rotate_psi_subspace_lcao");
}

template <>
void SpinConstrain<std::complex<double>>::cal_mi_lcao_subspace(
    const std::vector<std::vector<std::complex<double>>>& vcc_all,
    int nbands, int nk, int npol)
{
    ModuleBase::TITLE("SpinConstrain", "cal_mi_lcao_subspace");
    ModuleBase::timer::start("SpinConstrain", "cal_mi_lcao_subspace");

#ifdef __LCAO
    // ================================================================
    // DMR-based magnetic moment evaluation:
    //   1. Save original psi
    //   2. Rotate psi by subspace eigenvectors V(k)
    //   3. Build DMK and DMR from rotated psi
    //   4. Compute Mi via the standard cal_mi_lcao path (uses DMR)
    //   5. Restore original psi
    //
    // The trace formula (V†PV diagonal) was found to have systematic
    // bias compared to this DMR approach, so we revert to the more
    // accurate DMR-based evaluation.
    // ================================================================

    psi::Psi<std::complex<double>>* psi_t
        = static_cast<psi::Psi<std::complex<double>>*>(this->psi);

    const int nrow = this->ParaV->nrow;
    const int nloc_wfc = this->ParaV->nloc_wfc;

    // 1. Save original psi for all k-points
    std::vector<std::vector<std::complex<double>>> psi_save(nk);
    for (int ik = 0; ik < nk; ik++)
    {
        psi_t->fix_k(ik);
        std::complex<double>* ptr = psi_t->get_pointer();
        psi_save[ik].assign(ptr, ptr + nloc_wfc);
    }

    // 2. Rotate psi: C_new(k) = C_old(k) . V(k)
    this->rotate_psi_subspace_lcao(*psi_t, this->ParaV, vcc_all, nbands, nrow, nk);

    // 3. Build DMK and DMR from rotated psi
    elecstate::cal_dm_psi(this->ParaV, this->pelec->wg, *psi_t, *this->dm_);
    this->dm_->cal_DMR();

    // 4. Compute Mi using the standard LCAO path (reads DMR)
    this->cal_mi_lcao(0);

    // 5. Restore original psi
    for (int ik = 0; ik < nk; ik++)
    {
        psi_t->fix_k(ik);
        std::complex<double>* ptr = psi_t->get_pointer();
        std::memcpy(ptr, psi_save[ik].data(), sizeof(std::complex<double>) * nloc_wfc);
    }
#else
    this->zero_Mi();
    (void)vcc_all;
    (void)nbands;
    (void)nk;
    (void)npol;
#endif

    ModuleBase::timer::end("SpinConstrain", "cal_mi_lcao_subspace");
}

} // namespace spinconstrain