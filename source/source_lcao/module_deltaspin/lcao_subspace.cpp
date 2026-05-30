/**
 * @file lcao_subspace.cpp
 * @brief LCAO subspace diagonalization for DeltaSpin lambda loop optimization.
 *
 * @par Purpose
 * Implements the subspace diagonalization strategy for LCAO basis in DeltaSpin,
 * with mixed-precision support: hot GEMM operations can run in fp32 while
 * phase-sensitive operations (folding_HR) remain in fp64.
 *
 * @par Mixed Precision Strategy
 * - Input/Output: always complex<double> (fp64)
 * - Inner GEMM computations: complex<float> (fp32) when enabled
 * - Cast-down before GEMM, Cast-up after GEMM
 * - folding_HR always in fp64 (phase-sensitive)
 * - Controlled by subspace_exec_precision_ from GintPrecisionController
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
#include "source_lcao/module_hcontainer/hcontainer_funcs.h"

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
    lcao_PI_sub_diag_.clear();
    h_sub_local_buf_.clear();
    h_tmp_buf_.clear();
    s_tmp_buf_.clear();
    vcc_buf_.clear();
    s_copy_buf_.clear();
    eigenvalues_buf_.clear();
    lcao_ekb_save_.clear();
    lcao_lambda_in_sub_.clear();
    lcao_subspace_initialized_ = false;
}

void gather_sub_matrix_to_all(
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

void scatter_sub_matrix_to_local(
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

void extract_diagonal_from_local_block(
    const std::complex<double>* local_data,
    const Parallel_Orbitals* ParaV,
    double* diag,
    int nbands)
{
#ifdef __MPI
    int nrow = ParaV->nrow;
    int ncol_bands = ParaV->ncol_bands;
    std::fill(diag, diag + nbands, 0.0);
    for (int j_local = 0; j_local < ncol_bands; j_local++)
    {
        int j_global = ParaV->local2global_col(j_local);
        if (j_global >= nbands) continue;
        for (int i_local = 0; i_local < nrow; i_local++)
        {
            int i_global = ParaV->local2global_row(i_local);
            if (i_global == j_global && i_global < nbands)
            {
                diag[i_global] = local_data[i_local + j_local * nrow].real();
            }
        }
    }
    Parallel_Reduce::reduce_pool(diag, nbands);
#else
    for (int i = 0; i < nbands; i++)
    {
        diag[i] = local_data[i + i * nbands].real();
    }
#endif
}

static void cast_down_to_float(const std::complex<double>* src, std::complex<float>* dst, std::size_t n)
{
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096)
#endif
    for (std::size_t i = 0; i < n; ++i)
    {
        dst[i] = static_cast<std::complex<float>>(src[i]);
    }
}

static void cast_up_to_double(const std::complex<float>* src, std::complex<double>* dst, std::size_t n)
{
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096)
#endif
    for (std::size_t i = 0; i < n; ++i)
    {
        dst[i] = static_cast<std::complex<double>>(src[i]);
    }
}

static void gather_sub_matrix_to_all_float(
    const std::complex<float>* local_data,
    const Parallel_Orbitals* ParaV,
    std::complex<double>* full_data,
    int nbands)
{
    std::vector<std::complex<double>> local_db(ParaV->nrow * ParaV->ncol_bands);
    cast_up_to_double(local_data, local_db.data(), local_db.size());
    gather_sub_matrix_to_all(local_db.data(), ParaV, full_data, nbands);
}

static void scatter_sub_matrix_to_local_float(
    const std::complex<double>* full_data,
    const Parallel_Orbitals* ParaV,
    std::complex<float>* local_data,
    int nbands)
{
    int nrow = ParaV->nrow;
    int ncol_bands = ParaV->ncol_bands;
    int nloc_eij_local = nrow * ncol_bands;
#ifdef __MPI
    std::vector<std::complex<double>> local_db(nloc_eij_local, {0.0, 0.0});
    scatter_sub_matrix_to_local(full_data, ParaV, local_db.data(), nbands);
    cast_down_to_float(local_db.data(), local_data, nloc_eij_local);
#else
    std::vector<std::complex<double>> tmp(nbands * nbands);
    scatter_sub_matrix_to_local(full_data, ParaV, tmp.data(), nbands);
    cast_down_to_float(tmp.data(), local_data, nbands * nbands);
#endif
}

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

    const bool use_fp32 = (subspace_exec_precision_ == ModuleGint::GintPrecision::fp32);

    if (use_fp32)
    {
        ModuleBase::timer::start("SpinConstrain", "calc_sub_hs_fp32");
        const std::complex<float> one_f = {1.0f, 0.0f};
        const std::complex<float> zero_f = {0.0f, 0.0f};

#ifdef __MPI
        int nloc_wfc = ParaV->nloc_wfc;
        int lld_Eij = ParaV->nrow;
        int ncol_bands = ParaV->ncol_bands;
        int nloc_eij_local = lld_Eij * ncol_bands;

        std::vector<std::complex<float>> h_f(ParaV->nloc, zero_f);
        cast_down_to_float(h_mat.p, h_f.data(), ParaV->nloc);
        std::vector<std::complex<float>> s_f(ParaV->nloc, zero_f);
        cast_down_to_float(s_mat.p, s_f.data(), ParaV->nloc);
        std::vector<std::complex<float>> psi_f(nloc_wfc, zero_f);
        cast_down_to_float(psi_ptr, psi_f.data(), nloc_wfc);

        std::vector<std::complex<float>> temp_f(nloc_wfc, zero_f);
        ScalapackConnector::gemm('N', 'N', nlocal, nbands, nlocal,
            one_f, h_f.data(), 1, 1, ParaV->desc,
            psi_f.data(), 1, 1, ParaV->desc_wfc,
            zero_f, temp_f.data(), 1, 1, ParaV->desc_wfc);

        std::vector<std::complex<float>> h_sub_local_f(nloc_eij_local, zero_f);
        ScalapackConnector::gemm('C', 'N', nbands, nbands, nlocal,
            one_f, psi_f.data(), 1, 1, ParaV->desc_wfc,
            temp_f.data(), 1, 1, ParaV->desc_wfc,
            zero_f, h_sub_local_f.data(), 1, 1, ParaV->desc_Eij);

        gather_sub_matrix_to_all_float(h_sub_local_f.data(), ParaV, h_sub, nbands);

        std::vector<std::complex<float>> temp_s_f(nloc_wfc, zero_f);
        ScalapackConnector::gemm('N', 'N', nlocal, nbands, nlocal,
            one_f, s_f.data(), 1, 1, ParaV->desc,
            psi_f.data(), 1, 1, ParaV->desc_wfc,
            zero_f, temp_s_f.data(), 1, 1, ParaV->desc_wfc);

        std::vector<std::complex<float>> s_sub_local_f(nloc_eij_local, zero_f);
        ScalapackConnector::gemm('C', 'N', nbands, nbands, nlocal,
            one_f, psi_f.data(), 1, 1, ParaV->desc_wfc,
            temp_s_f.data(), 1, 1, ParaV->desc_wfc,
            zero_f, s_sub_local_f.data(), 1, 1, ParaV->desc_Eij);

        gather_sub_matrix_to_all_float(s_sub_local_f.data(), ParaV, s_sub, nbands);
#else
        const std::complex<double> one_d = {1.0, 0.0};
        const std::complex<double> zero_d = {0.0, 0.0};

        std::vector<std::complex<float>> h_f(nlocal * nlocal);
        cast_down_to_float(h_mat.p, h_f.data(), nlocal * nlocal);
        std::vector<std::complex<float>> s_f(nlocal * nlocal);
        cast_down_to_float(s_mat.p, s_f.data(), nlocal * nlocal);
        std::vector<std::complex<float>> psi_f(nlocal * nbands);
        cast_down_to_float(psi_ptr, psi_f.data(), nlocal * nbands);

        std::vector<std::complex<float>> temp_f(nlocal * nbands, zero_f);
        cgemm_("N", "N", &nlocal, &nbands, &nlocal,
            &one_f, h_f.data(), &nlocal, psi_f.data(), &nlocal, &zero_f, temp_f.data(), &nlocal);

        std::vector<std::complex<float>> h_sub_f(nbands * nbands, zero_f);
        cgemm_("C", "N", &nbands, &nbands, &nlocal,
            &one_f, psi_f.data(), &nlocal, temp_f.data(), &nlocal, &zero_f, h_sub_f.data(), &nbands);
        cast_up_to_double(h_sub_f.data(), h_sub, nbands * nbands);

        std::vector<std::complex<float>> temp_s_f(nlocal * nbands, zero_f);
        cgemm_("N", "N", &nlocal, &nbands, &nlocal,
            &one_f, s_f.data(), &nlocal, psi_f.data(), &nlocal, &zero_f, temp_s_f.data(), &nlocal);

        std::vector<std::complex<float>> s_sub_f(nbands * nbands, zero_f);
        cgemm_("C", "N", &nbands, &nbands, &nlocal,
            &one_f, psi_f.data(), &nlocal, temp_s_f.data(), &nlocal, &zero_f, s_sub_f.data(), &nbands);
        cast_up_to_double(s_sub_f.data(), s_sub, nbands * nbands);
#endif
        ModuleBase::timer::end("SpinConstrain", "calc_sub_hs_fp32");
    }
    else
    {
        const std::complex<double> one = {1.0, 0.0};
        const std::complex<double> zero = {0.0, 0.0};

#ifdef __MPI
        int nloc_wfc = ParaV->nloc_wfc;
        std::vector<std::complex<double>> temp(nloc_wfc, zero);

        ScalapackConnector::gemm('N', 'N', nlocal, nbands, nlocal,
            one, h_mat.p, 1, 1, ParaV->desc,
            psi_ptr, 1, 1, ParaV->desc_wfc,
            zero, temp.data(), 1, 1, ParaV->desc_wfc);

        int lld_Eij = ParaV->nrow;
        int ncol_bands = ParaV->ncol_bands;
        int nloc_eij_local = lld_Eij * ncol_bands;

        std::vector<std::complex<double>> h_sub_local(nloc_eij_local, zero);

        ScalapackConnector::gemm('C', 'N', nbands, nbands, nlocal,
            one, psi_ptr, 1, 1, ParaV->desc_wfc,
            temp.data(), 1, 1, ParaV->desc_wfc,
            zero, h_sub_local.data(), 1, 1, ParaV->desc_Eij);

        gather_sub_matrix_to_all(h_sub_local.data(), ParaV, h_sub, nbands);

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
        std::vector<std::complex<double>> temp(nlocal * nbands, zero);
        zgemm_("N", "N", &nlocal, &nbands, &nlocal,
            &one, h_mat.p, &nlocal, psi_ptr, &nlocal, &zero, temp.data(), &nlocal);

        zgemm_("C", "N", &nbands, &nbands, &nlocal,
            &one, psi_ptr, &nlocal, temp.data(), &nlocal, &zero, h_sub, &nbands);

        std::vector<std::complex<double>> temp_s(nlocal * nbands, zero);
        zgemm_("N", "N", &nlocal, &nbands, &nlocal,
            &one, s_mat.p, &nlocal, psi_ptr, &nlocal, &zero, temp_s.data(), &nlocal);

        zgemm_("C", "N", &nbands, &nbands, &nlocal,
            &one, psi_ptr, &nlocal, temp_s.data(), &nlocal, &zero, s_sub, &nbands);
#endif
    }

    ModuleBase::timer::end("SpinConstrain", "calculate_lcao_sub_hs");
}

template <>
void SpinConstrain<std::complex<double>>::calculate_PI_sub_from_hr(
    const hamilt::HContainer<double>* pre_hr_iat,
    psi::Psi<std::complex<double>>& psi,
    const Parallel_Orbitals* ParaV,
    const ModuleBase::Vector3<double>& kvec_d,
    std::complex<double>* PI_sub_local_out,
    int nbands, int nlocal)
{
    ModuleBase::TITLE("SpinConstrain", "calculate_PI_sub_from_hr");
    ModuleBase::timer::start("SpinConstrain", "calculate_PI_sub_from_hr");

    std::complex<double>* psi_ptr = psi.get_pointer();

    const bool use_fp32 = (subspace_exec_precision_ == ModuleGint::GintPrecision::fp32);

#ifdef __MPI
    int nloc = ParaV->nloc;
    std::vector<std::complex<double>> PI_k_local(nloc, {0.0, 0.0});
    hamilt::folding_HR(*pre_hr_iat, PI_k_local.data(), kvec_d, ParaV->nrow, 1);

    if (use_fp32)
    {
        ModuleBase::timer::start("SpinConstrain", "calc_PI_sub_fp32");
        const std::complex<float> one_f = {1.0f, 0.0f};
        const std::complex<float> zero_f = {0.0f, 0.0f};
        int nloc_wfc = ParaV->nloc_wfc;
        int lld_Eij = ParaV->nrow;
        int ncol_bands = ParaV->ncol_bands;
        int nloc_eij_local = lld_Eij * ncol_bands;

        std::vector<std::complex<float>> PI_k_f(nloc, zero_f);
        cast_down_to_float(PI_k_local.data(), PI_k_f.data(), nloc);
        std::vector<std::complex<float>> psi_f(nloc_wfc, zero_f);
        cast_down_to_float(psi_ptr, psi_f.data(), nloc_wfc);

        std::vector<std::complex<float>> temp_f(nloc_wfc, zero_f);
        ScalapackConnector::gemm('N', 'N', nlocal, nbands, nlocal,
            one_f, PI_k_f.data(), 1, 1, ParaV->desc,
            psi_f.data(), 1, 1, ParaV->desc_wfc,
            zero_f, temp_f.data(), 1, 1, ParaV->desc_wfc);

        std::vector<std::complex<float>> PI_sub_local_f(nloc_eij_local, zero_f);
        ScalapackConnector::gemm('C', 'N', nbands, nbands, nlocal,
            one_f, psi_f.data(), 1, 1, ParaV->desc_wfc,
            temp_f.data(), 1, 1, ParaV->desc_wfc,
            zero_f, PI_sub_local_f.data(), 1, 1, ParaV->desc_Eij);

        cast_up_to_double(PI_sub_local_f.data(), PI_sub_local_out, nloc_eij_local);
        ModuleBase::timer::end("SpinConstrain", "calc_PI_sub_fp32");
    }
    else
#endif
    {
        const std::complex<double> one = {1.0, 0.0};
        const std::complex<double> zero = {0.0, 0.0};

#ifdef __MPI
        int nloc_wfc = ParaV->nloc_wfc;
        std::vector<std::complex<double>> temp(nloc_wfc, zero);

        ScalapackConnector::gemm('N', 'N', nlocal, nbands, nlocal,
            one, PI_k_local.data(), 1, 1, ParaV->desc,
            psi_ptr, 1, 1, ParaV->desc_wfc,
            zero, temp.data(), 1, 1, ParaV->desc_wfc);

        int lld_Eij = ParaV->nrow;
        int ncol_bands = ParaV->ncol_bands;
        int nloc_eij_local = lld_Eij * ncol_bands;
        std::vector<std::complex<double>> PI_sub_local(nloc_eij_local, zero);

        ScalapackConnector::gemm('C', 'N', nbands, nbands, nlocal,
            one, psi_ptr, 1, 1, ParaV->desc_wfc,
            temp.data(), 1, 1, ParaV->desc_wfc,
            zero, PI_sub_local.data(), 1, 1, ParaV->desc_Eij);

        std::memcpy(PI_sub_local_out, PI_sub_local.data(), sizeof(std::complex<double>) * nloc_eij_local);
#else
        std::vector<std::complex<double>> PI_k(nlocal * nlocal, {0.0, 0.0});
        hamilt::folding_HR(*pre_hr_iat, PI_k.data(), kvec_d, nlocal, 1);

        std::vector<std::complex<double>> temp(nlocal * nbands, zero);
        zgemm_("N", "N", &nlocal, &nbands, &nlocal,
            &one, PI_k.data(), &nlocal, psi_ptr, &nlocal, &zero, temp.data(), &nlocal);

        zgemm_("C", "N", &nbands, &nbands, &nlocal,
            &one, psi_ptr, &nlocal, temp.data(), &nlocal, &zero, PI_sub_local_out, &nbands);
#endif
    }

    ModuleBase::timer::end("SpinConstrain", "calculate_PI_sub_from_hr");
}

template <>
void SpinConstrain<std::complex<double>>::calculate_delta_hcc_lcao(
    std::complex<double>* h_sub_local,
    const std::map<int, std::vector<std::complex<double>>>& PI_sub_local,
    const ModuleBase::Vector3<double>* lambda,
    int nbands, int ik, bool full_update,
    const Parallel_Orbitals* ParaV)
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

    const int nloc_eij = ParaV->nrow * ParaV->ncol_bands;

    if (this->npol_ == 2)
    {
        for (const auto& [iat, pi_local] : PI_sub_local)
        {
            const std::complex<double> coeff0(effective_lambda[iat][2], 0.0);
            const std::complex<double>* pi_ptr = pi_local.data();
            for (int ij = 0; ij < nloc_eij; ij++)
            {
                h_sub_local[ij] += coeff0 * pi_ptr[ij];
            }
        }
    }
    else if (this->npol_ == 1)
    {
        int spin_sign = this->get_spin_sign(ik);

        for (const auto& [iat, pi_local] : PI_sub_local)
        {
            const std::complex<double> coeff(effective_lambda[iat][2] * spin_sign, 0.0);
            const std::complex<double>* pi_ptr = pi_local.data();
            for (int ij = 0; ij < nloc_eij; ij++)
            {
                h_sub_local[ij] += coeff * pi_ptr[ij];
            }
        }
    }

    ModuleBase::timer::end("SpinConstrain", "calculate_delta_hcc_lcao");
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

    const int nn = nbands * nbands;

    if (this->npol_ == 2)
    {
        for (int iat = 0; iat < nat; iat++)
        {
            if (PI_sub[iat].empty()) continue;

            const std::complex<double> coeff0(effective_lambda[iat][2], 0.0);
            const std::complex<double>* pi_ptr = PI_sub[iat].data();
            for (int ij = 0; ij < nn; ij++)
            {
                h_sub[ij] += coeff0 * pi_ptr[ij];
            }
        }
    }
    else if (this->npol_ == 1)
    {
        int spin_sign = this->get_spin_sign(ik);

        for (int iat = 0; iat < nat; iat++)
        {
            if (PI_sub[iat].empty()) continue;

            const std::complex<double> coeff(effective_lambda[iat][2] * spin_sign, 0.0);
            const std::complex<double>* pi_ptr = PI_sub[iat].data();
            for (int ij = 0; ij < nn; ij++)
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

    const bool use_fp32 = (subspace_exec_precision_ == ModuleGint::GintPrecision::fp32);

    if (use_fp32)
    {
        ModuleBase::timer::start("SpinConstrain", "rotate_psi_fp32");
        const std::complex<float> one_f = {1.0f, 0.0f};
        const std::complex<float> zero_f = {0.0f, 0.0f};

        for (int ik = 0; ik < nk; ik++)
        {
            psi.fix_k(ik);
            std::complex<double>* psi_ptr = psi.get_pointer();

#ifdef __MPI
            int nloc_wfc = ParaV->nloc_wfc;
            int lld_Eij = ParaV->nrow;
            int ncol_bands = ParaV->ncol_bands;
            int nloc_eij_local = lld_Eij * ncol_bands;

            std::vector<std::complex<float>> v_local_f(nloc_eij_local, zero_f);
            scatter_sub_matrix_to_local_float(vcc_all[ik].data(), ParaV, v_local_f.data(), nbands);

            std::vector<std::complex<float>> psi_f(nloc_wfc, zero_f);
            cast_down_to_float(psi_ptr, psi_f.data(), nloc_wfc);

            std::vector<std::complex<float>> temp_f(nloc_wfc, zero_f);
            ScalapackConnector::gemm('N', 'N', nlocal, nbands, nbands,
                one_f, psi_f.data(), 1, 1, ParaV->desc_wfc,
                v_local_f.data(), 1, 1, ParaV->desc_Eij,
                zero_f, temp_f.data(), 1, 1, ParaV->desc_wfc);

            cast_up_to_double(temp_f.data(), psi_ptr, nloc_wfc);
#else
            std::vector<std::complex<float>> psi_f(nlocal * nbands);
            cast_down_to_float(psi_ptr, psi_f.data(), nlocal * nbands);
            std::vector<std::complex<float>> v_f(nbands * nbands);
            cast_down_to_float(vcc_all[ik].data(), v_f.data(), nbands * nbands);

            std::vector<std::complex<float>> temp_f(nlocal * nbands, zero_f);
            cgemm_("N", "N", &nlocal, &nbands, &nbands,
                &one_f, psi_f.data(), &nlocal,
                v_f.data(), &nbands,
                &zero_f, temp_f.data(), &nlocal);
            cast_up_to_double(temp_f.data(), psi_ptr, nlocal * nbands);
#endif
        }
        ModuleBase::timer::end("SpinConstrain", "rotate_psi_fp32");
    }
    else
    {
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

            std::vector<std::complex<double>> v_local(nloc_eij_local, zero);
            scatter_sub_matrix_to_local(vcc_all[ik].data(), ParaV, v_local.data(), nbands);

            std::vector<std::complex<double>> temp(nloc_wfc, zero);

            ScalapackConnector::gemm('N', 'N', nlocal, nbands, nbands,
                one, psi_ptr, 1, 1, ParaV->desc_wfc,
                v_local.data(), 1, 1, ParaV->desc_Eij,
                zero, temp.data(), 1, 1, ParaV->desc_wfc);

            std::memcpy(psi_ptr, temp.data(), sizeof(std::complex<double>) * nloc_wfc);
#else
            std::vector<std::complex<double>> temp(nlocal * nbands);
            zgemm_("N", "N", &nlocal, &nbands, &nbands,
                &one, psi_ptr, &nlocal,
                vcc_all[ik].data(), &nbands,
                &zero, temp.data(), &nlocal);
            std::memcpy(psi_ptr, temp.data(), sizeof(std::complex<double>) * nlocal * nbands);
#endif
        }
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
    psi::Psi<std::complex<double>>* psi_t
        = static_cast<psi::Psi<std::complex<double>>*>(this->psi);

    const int nrow = this->ParaV->nrow;
    const int nlocal = this->ParaV->get_global_row_size();
    const int nloc_wfc = this->ParaV->nloc_wfc;

    std::vector<std::vector<std::complex<double>>> psi_save(nk);
    for (int ik = 0; ik < nk; ik++)
    {
        psi_t->fix_k(ik);
        std::complex<double>* ptr = psi_t->get_pointer();
        psi_save[ik].assign(ptr, ptr + nloc_wfc);
    }

    this->rotate_psi_subspace_lcao(*psi_t, this->ParaV, vcc_all, nbands, nlocal, nk);

    if (this->subspace_exec_precision_ == ModuleGint::GintPrecision::fp32)
    {
        elecstate::cal_dm_psi_mixed(this->ParaV, this->pelec->wg, *psi_t, *this->dm_);
    }
    else
    {
        elecstate::cal_dm_psi(this->ParaV, this->pelec->wg, *psi_t, *this->dm_);
    }
    this->dm_->cal_DMR();

    this->cal_mi_lcao(0);

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
