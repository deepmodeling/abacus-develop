#include "source_hsolver/kernels/bpcg_kernel_op.h"
#include "source_base/module_external/blas_connector.h"
#include "source_base/kernels/math_kernel_op.h"
#include "source_base/parallel_reduce.h"
#include <vector>
namespace hsolver
{

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
        std::vector<Real> norms(n_band, 0.0);
        std::vector<Real> epsilo_0s(n_band, 0.0);
        std::vector<Real> epsilo_1s(n_band, 0.0);
        std::vector<Real> epsilo_2s(n_band, 0.0);

#ifdef _OPENMP
#pragma omp parallel if(n_band > 4)
#endif
        {
            // Step 1: compute norms for all bands
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
            for (int band_idx = 0; band_idx < n_band; band_idx++)
            {
                Real norm = 0.0;
                for (int basis_idx = 0; basis_idx < n_basis; basis_idx++)
                {
                    auto item = band_idx * n_basis_max + basis_idx;
                    norm += std::norm(grad_out[item]);
                }
                norms[band_idx] = norm;
            }

            // Step 2: reduce norms serially
#ifdef _OPENMP
#pragma omp single
#endif
            {
                for (int band_idx = 0; band_idx < n_band; band_idx++)
                {
                    Parallel_Reduce::reduce_pool(norms[band_idx]);
                    norms[band_idx] = 1.0 / sqrt(norms[band_idx]);
                }
            }

            // Step 3: normalize and compute epsilo for all bands
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
            for (int band_idx = 0; band_idx < n_band; band_idx++)
            {
                Real norm = norms[band_idx];
                Real epsilo_0 = 0.0, epsilo_1 = 0.0, epsilo_2 = 0.0;
                for (int basis_idx = 0; basis_idx < n_basis; basis_idx++)
                {
                    auto item = band_idx * n_basis_max + basis_idx;
                    grad_out[item] *= norm;
                    hgrad_out[item] *= norm;
                    epsilo_0 += std::real(hpsi_out[item] * std::conj(psi_out[item]));
                    epsilo_1 += std::real(grad_out[item] * std::conj(hpsi_out[item]));
                    epsilo_2 += std::real(grad_out[item] * std::conj(hgrad_out[item]));
                }
                epsilo_0s[band_idx] = epsilo_0;
                epsilo_1s[band_idx] = epsilo_1;
                epsilo_2s[band_idx] = epsilo_2;
            }

            // Step 4: reduce epsilos serially
#ifdef _OPENMP
#pragma omp single
#endif
            {
                for (int band_idx = 0; band_idx < n_band; band_idx++)
                {
                    Parallel_Reduce::reduce_pool(epsilo_0s[band_idx]);
                    Parallel_Reduce::reduce_pool(epsilo_1s[band_idx]);
                    Parallel_Reduce::reduce_pool(epsilo_2s[band_idx]);
                }
            }

            // Step 5: update psi and hpsi for all bands
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
            for (int band_idx = 0; band_idx < n_band; band_idx++)
            {
                Real theta = 0.5 * std::abs(std::atan(2 * epsilo_1s[band_idx] / (epsilo_0s[band_idx] - epsilo_2s[band_idx])));
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
    }
};

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
        std::vector<Real> norms(n_band, 0.0);
        std::vector<Real> epsilos(n_band, 0.0);
        std::vector<Real> errs(n_band, 0.0);
        std::vector<Real> betas(n_band, 0.0);

#ifdef _OPENMP
#pragma omp parallel if(n_band > 4)
#endif
        {
            // Step 1: compute norms for all bands
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
            for (int band_idx = 0; band_idx < n_band; band_idx++)
            {
                Real norm = 0.0;
                for (int basis_idx = 0; basis_idx < n_basis; basis_idx++)
                {
                    auto item = band_idx * n_basis_max + basis_idx;
                    norm += std::norm(psi_out[item]);
                }
                norms[band_idx] = norm;
            }

            // Step 2: reduce norms serially
#ifdef _OPENMP
#pragma omp single
#endif
            {
                for (int band_idx = 0; band_idx < n_band; band_idx++)
                {
                    Parallel_Reduce::reduce_pool(norms[band_idx]);
                    norms[band_idx] = 1.0 / sqrt(norms[band_idx]);
                }
            }

            // Step 3: normalize and compute epsilo for all bands
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
            for (int band_idx = 0; band_idx < n_band; band_idx++)
            {
                Real norm = norms[band_idx];
                Real epsilo = 0.0;
                for (int basis_idx = 0; basis_idx < n_basis; basis_idx++)
                {
                    auto item = band_idx * n_basis_max + basis_idx;
                    psi_out[item] *= norm;
                    hpsi_out[item] *= norm;
                    epsilo += std::real(hpsi_out[item] * std::conj(psi_out[item]));
                }
                epsilos[band_idx] = epsilo;
            }

            // Step 4: reduce epsilos serially
#ifdef _OPENMP
#pragma omp single
#endif
            {
                for (int band_idx = 0; band_idx < n_band; band_idx++)
                {
                    Parallel_Reduce::reduce_pool(epsilos[band_idx]);
                }
            }

            // Step 5: compute err and beta for all bands
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
            for (int band_idx = 0; band_idx < n_band; band_idx++)
            {
                Real epsilo = epsilos[band_idx];
                Real err = 0.0;
                Real beta = 0.0;
                for (int basis_idx = 0; basis_idx < n_basis; basis_idx++)
                {
                    auto item = band_idx * n_basis_max + basis_idx;
                    T grad_1 = hpsi_out[item] - epsilo * psi_out[item];
                    Real grad_2 = std::norm(grad_1);
                    err += grad_2;
                    beta += grad_2 / prec_in[basis_idx];
                }
                errs[band_idx] = err;
                betas[band_idx] = beta;
            }

            // Step 6: reduce errs and betas serially
#ifdef _OPENMP
#pragma omp single
#endif
            {
                for (int band_idx = 0; band_idx < n_band; band_idx++)
                {
                    Parallel_Reduce::reduce_pool(errs[band_idx]);
                    Parallel_Reduce::reduce_pool(betas[band_idx]);
                }
            }

            // Step 7: update grad and output err/beta for all bands
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
            for (int band_idx = 0; band_idx < n_band; band_idx++)
            {
                Real epsilo = epsilos[band_idx];
                Real beta = betas[band_idx];
                for (int basis_idx = 0; basis_idx < n_basis; basis_idx++)
                {
                    auto item = band_idx * n_basis_max + basis_idx;
                    T grad_1 = hpsi_out[item] - epsilo * psi_out[item];
                    grad_out[item] = -grad_1 / prec_in[basis_idx] + beta / beta_out[band_idx] * grad_old_out[item];
                }
                beta_out[band_idx] = beta;
                err_out[band_idx] = sqrt(errs[band_idx]);
            }
        }
    }
};

template <typename T>
struct apply_eigenvalues_op<T, base_device::DEVICE_CPU>
{
    using Real = typename GetTypeReal<T>::type;
    void operator()(const int& nbase, const int& nbase_x, const int& notconv, T* result, const T* vectors, const Real* eigenvalues)
    {
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static) if(notconv * nbase > 1024)
#endif
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
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(notconv > 4)
#endif
        for (int m = 0; m < notconv; m++)
        {
            std::vector<Real> pre(dim, 0.0);
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
        std::vector<Real> norms(notconv, 0.0);

#ifdef _OPENMP
#pragma omp parallel if(notconv > 4)
#endif
        {
            // Step 1: compute norms for all bands in parallel
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
            for (int m = 0; m < notconv; m++)
            {
                Real norm = 0.0;
                T* psi_m = psi_iter + (nbase + m) * dim;
                for (int i = 0; i < dim; i++)
                {
                    norm += std::norm(psi_m[i]);
                }
                norms[m] = norm;
            }

            // Step 2: reduce norms serially (MPI calls inside OpenMP must be serial)
#ifdef _OPENMP
#pragma omp single
#endif
            {
                for (int m = 0; m < notconv; m++)
                {
                    Parallel_Reduce::reduce_pool(norms[m]);
                    norms[m] = sqrt(norms[m]);
                }
            }

            // Step 3: normalize all bands in parallel
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
            for (int m = 0; m < notconv; m++)
            {
                Real psi_m_norm = norms[m];
                T* psi_m = psi_iter + (nbase + m) * dim;
                for (int i = 0; i < dim; i++)
                {
                    psi_m[i] /= psi_m_norm;
                }
                if (psi_norm) {
                    psi_norm[m] = psi_m_norm;
                }
            }
        }
    }
};

template <typename T>
struct refresh_hcc_scc_vcc_op<T, base_device::DEVICE_CPU>
{
    using Real = typename GetTypeReal<T>::type;
    void operator()(const int &n,
                  T *hcc,
                  T *scc,
                  T *vcc,
                  const int &ldh,
                  const Real *eigenvalue,
                  const T &one)
    {
#ifdef _OPENMP
#pragma omp parallel for collapse(1) schedule(static) if(n > 64)
#endif
        for (int i = 0; i < n; i++)
        {
            hcc[i * ldh + i] = eigenvalue[i];
            scc[i * ldh + i] = one;
            vcc[i * ldh + i] = one;
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
template struct refresh_hcc_scc_vcc_op<std::complex<float>, base_device::DEVICE_CPU>;
template struct refresh_hcc_scc_vcc_op<std::complex<double>, base_device::DEVICE_CPU>;
template struct refresh_hcc_scc_vcc_op<double, base_device::DEVICE_CPU>;
} // namespace hsolver