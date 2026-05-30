#include "source_hsolver/kernels/bpcg_kernel_op.h"
#include "source_base/module_external/blas_connector.h"
#include "source_base/kernels/math_kernel_op.h"
#include "source_base/parallel_reduce.h"
#include <vector>
namespace hsolver
{

// ========== 优化后的 line_minimize_with_block_op ==========
template <typename T>
struct line_minimize_with_block_op<T, base_device::DEVICE_CPU>
{
    using Real = typename GetTypeReal<T>::type;
    void operator()(T* grad_out,
                    T* hgrad_out,
                    T* psi_out,
                    T* hpsi_out,
                    const int& n_basis,
                    const int& n_basis_max,
                    const int& n_band)
    {
        // 存储每个 band 的中间结果
        std::vector<Real> norms(n_band, 0.0);
        std::vector<Real> epsilo_0(n_band, 0.0);
        std::vector<Real> epsilo_1(n_band, 0.0);
        std::vector<Real> epsilo_2(n_band, 0.0);

        // ========== Phase 1: 并行计算 norm ==========
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 8)
#endif
        for (int band_idx = 0; band_idx < n_band; band_idx++)
        {
            auto A = reinterpret_cast<const Real*>(grad_out + band_idx * n_basis_max);
            norms[band_idx] = BlasConnector::dot(2 * n_basis, A, 1, A, 1);
        }

        // 归一化
        for (int i = 0; i < n_band; i++) {
            Parallel_Reduce::reduce_pool(norms[i]);
            norms[i] = 1.0 / sqrt(norms[i]);
        }

        // ========== Phase 2: 并行归一化并计算 epsilo ==========
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 8)
#endif
        for (int band_idx = 0; band_idx < n_band; band_idx++)
        {
            Real norm = norms[band_idx];

            for (int basis_idx = 0; basis_idx < n_basis; basis_idx++)
            {
                auto item = band_idx * n_basis_max + basis_idx;
                grad_out[item] *= norm;
                hgrad_out[item] *= norm;
                epsilo_0[band_idx] += std::real(hpsi_out[item] * std::conj(psi_out[item]));
                epsilo_1[band_idx] += std::real(grad_out[item] * std::conj(hpsi_out[item]));
                epsilo_2[band_idx] += std::real(grad_out[item] * std::conj(hgrad_out[item]));
            }
        }

        // 归一化 epsilo
        for (int i = 0; i < n_band; i++) {
            Parallel_Reduce::reduce_pool(epsilo_0[i]);
            Parallel_Reduce::reduce_pool(epsilo_1[i]);
            Parallel_Reduce::reduce_pool(epsilo_2[i]);
        }

        // ========== Phase 3: 并行更新 psi 和 hpsi ==========
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 8)
#endif
        for (int band_idx = 0; band_idx < n_band; band_idx++)
        {
            Real theta = 0.5 * std::abs(std::atan(2 * epsilo_1[band_idx] /
                                                  (epsilo_0[band_idx] - epsilo_2[band_idx])));
            Real cos_theta = std::cos(theta);
            Real sin_theta = std::sin(theta);

            for (int basis_idx = 0; basis_idx < n_basis; basis_idx++)
            {
                auto item = band_idx * n_basis_max + basis_idx;
                psi_out[item] = psi_out[item] * cos_theta + grad_out[item] * sin_theta;
                hpsi_out[item] = hpsi_out[item] * cos_theta + hgrad_out[item] * sin_theta;
            }
        }
    }
};

// ========== 优化后的 calc_grad_with_block_op ==========
template <typename T>
struct calc_grad_with_block_op<T, base_device::DEVICE_CPU>
{
    using Real = typename GetTypeReal<T>::type;
    void operator()(const Real* prec_in,
                    Real* err_out,
                    Real* beta_out,
                    T* psi_out,
                    T* hpsi_out,
                    T* grad_out,
                    T* grad_old_out,
                    const int& n_basis,
                    const int& n_basis_max,
                    const int& n_band)
    {
        // 存储每个 band 的中间结果
        std::vector<Real> norms(n_band, 0.0);
        std::vector<Real> epsilos(n_band, 0.0);
        std::vector<Real> errs(n_band, 0.0);
        std::vector<Real> betas(n_band, 0.0);

        // ========== Phase 1: 并行计算 norm ==========
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 8)
#endif
        for (int band_idx = 0; band_idx < n_band; band_idx++)
        {
            auto A = reinterpret_cast<const Real*>(psi_out + band_idx * n_basis_max);
            norms[band_idx] = BlasConnector::dot(2 * n_basis, A, 1, A, 1);
        }

        // 归一化
        for (int i = 0; i < n_band; i++) {
            Parallel_Reduce::reduce_pool(norms[i]);
            norms[i] = 1.0 / sqrt(norms[i]);
        }

        // ========== Phase 2: 并行归一化并计算 epsilo ==========
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 8)
#endif
        for (int band_idx = 0; band_idx < n_band; band_idx++)
        {
            Real norm = norms[band_idx];

            for (int basis_idx = 0; basis_idx < n_basis; basis_idx++)
            {
                auto item = band_idx * n_basis_max + basis_idx;
                psi_out[item] *= norm;
                hpsi_out[item] *= norm;
                epsilos[band_idx] += std::real(hpsi_out[item] * std::conj(psi_out[item]));
            }
        }

        // 归一化 epsilo
        for (int i = 0; i < n_band; i++) {
            Parallel_Reduce::reduce_pool(epsilos[i]);
        }

        // ========== Phase 3: 并行计算 err 和 beta ==========
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 8)
#endif
        for (int band_idx = 0; band_idx < n_band; band_idx++)
        {
            Real epsilo = epsilos[band_idx];
            T grad_1 = {0.0, 0.0};

            for (int basis_idx = 0; basis_idx < n_basis; basis_idx++)
            {
                auto item = band_idx * n_basis_max + basis_idx;
                grad_1 = hpsi_out[item] - epsilo * psi_out[item];
                Real grad_2 = std::norm(grad_1);
                errs[band_idx] += grad_2;
                betas[band_idx] += grad_2 / prec_in[basis_idx];
            }
        }

        // 归一化 err 和 beta
        for (int i = 0; i < n_band; i++) {
            Parallel_Reduce::reduce_pool(errs[i]);
            Parallel_Reduce::reduce_pool(betas[i]);
        }

        // ========== Phase 4: 并行计算最终梯度 ==========
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 8)
#endif
        for (int band_idx = 0; band_idx < n_band; band_idx++)
        {
            Real epsilo = epsilos[band_idx];
            Real beta = betas[band_idx];
            T grad_1 = {0.0, 0.0};

            for (int basis_idx = 0; basis_idx < n_basis; basis_idx++)
            {
                auto item = band_idx * n_basis_max + basis_idx;
                grad_1 = hpsi_out[item] - epsilo * psi_out[item];
                grad_out[item] = -grad_1 / prec_in[basis_idx]
                               + beta / beta_out[band_idx] * grad_old_out[item];
            }
            beta_out[band_idx] = beta;
            err_out[band_idx] = sqrt(errs[band_idx]);
        }
    }
};

template <typename T>
struct apply_eigenvalues_op<T, base_device::DEVICE_CPU>
{
    using Real = typename GetTypeReal<T>::type;
    void operator()(const int& nbase, const int& nbase_x, const int& notconv, T* result, const T* vectors, const Real* eigenvalues)
    {
        for (int m = 0; m < notconv; m++)
        {
            for (int idx = 0; idx < nbase; idx++)
            {
                result[m * nbase_x + idx] = eigenvalues[m] * vectors[m * nbase_x + idx];
            }
        }
    }
};

template <typename T>
struct precondition_op<T, base_device::DEVICE_CPU> {
    using Real = typename GetTypeReal<T>::type;
    void operator()(const int& dim,
                   T* psi_iter,
                   const int& nbase,
                   const int& notconv,
                   const Real* precondition,
                   const Real* eigenvalues)
    {
        std::vector<Real> pre(dim, 0.0);
        for (int m = 0; m < notconv; m++)
        {
            for (size_t i = 0; i < dim; i++)
            {
                Real x = std::abs(precondition[i] - eigenvalues[m]);
                pre[i] = 0.5 * (1.0 + x + sqrt(1 + (x - 1.0) * (x - 1.0)));
            }
            ModuleBase::vector_div_vector_op<T, base_device::DEVICE_CPU>()(
                                                             dim,
                                                             psi_iter + (nbase + m) * dim,
                                                             psi_iter + (nbase + m) * dim,
                                                             pre.data());
        }
    }
};

template <typename T>
struct normalize_op<T, base_device::DEVICE_CPU> {
    void operator()(const int& dim,
                   T* psi_iter,
                   const int& nbase,
                   const int& notconv,
                   typename GetTypeReal<T>::type* psi_norm)
    {
        using Real = typename GetTypeReal<T>::type;
        for (int m = 0; m < notconv; m++)
        {
            // Calculate norm using dot_real_op
            Real psi_m_norm = ModuleBase::dot_real_op<T, base_device::DEVICE_CPU>()(
                                                                dim,
                                                                psi_iter + (nbase + m) * dim,
                                                                psi_iter + (nbase + m) * dim,
                                                                true);
            assert(psi_m_norm > 0.0);
            psi_m_norm = sqrt(psi_m_norm);

            // Normalize using vector_div_constant_op
            ModuleBase::vector_div_constant_op<T, base_device::DEVICE_CPU>()(
                                                              dim,
                                                              psi_iter + (nbase + m) * dim,
                                                              psi_iter + (nbase + m) * dim,
                                                              psi_m_norm);
            if (psi_norm) {
                psi_norm[m] = psi_m_norm;
            }
        }
    }
};

template struct calc_grad_with_block_op<std::complex<float>, base_device::DEVICE_CPU>;
template struct line_minimize_with_block_op<std::complex<float>, base_device::DEVICE_CPU>;
template struct calc_grad_with_block_op<std::complex<double>, base_device::DEVICE_CPU>;
template struct line_minimize_with_block_op<std::complex<double>, base_device::DEVICE_CPU>;
template struct apply_eigenvalues_op<std::complex<float>, base_device::DEVICE_CPU>;
template struct apply_eigenvalues_op<std::complex<double>, base_device::DEVICE_CPU>;
template struct apply_eigenvalues_op<double, base_device::DEVICE_CPU>;
template struct precondition_op<std::complex<float>, base_device::DEVICE_CPU>;
template struct precondition_op<std::complex<double>, base_device::DEVICE_CPU>;
template struct precondition_op<double, base_device::DEVICE_CPU>;
template struct normalize_op<std::complex<float>, base_device::DEVICE_CPU>;
template struct normalize_op<std::complex<double>, base_device::DEVICE_CPU>;
template struct normalize_op<double, base_device::DEVICE_CPU>;
} // namespace hsolver
