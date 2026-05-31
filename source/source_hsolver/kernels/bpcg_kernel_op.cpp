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
        // OpenMP parallelization over bands: each band accesses a disjoint memory region
        // [band_idx * n_basis_max, (band_idx+1) * n_basis_max), so no data races occur.
        // MPI collective calls (reduce_pool) use per-band scalar reductions executed serially
        // outside parallel regions to avoid thread-safety issues with MPI.
        // Static scheduling is used since each band has equal workload (n_basis).
        std::vector<Real> norms(n_band);
        std::vector<Real> epsilo_0_arr(n_band);
        std::vector<Real> epsilo_1_arr(n_band);
        std::vector<Real> epsilo_2_arr(n_band);

        // Phase 1: parallel computation of per-band norms via BLAS dot
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int band_idx = 0; band_idx < n_band; band_idx++)
        {
            auto A = reinterpret_cast<const Real*>(grad_out + band_idx * n_basis_max);
            norms[band_idx] = BlasConnector::dot(2 * n_basis, A, 1, A, 1);
        }

        // Phase 2: MPI reduction of norms across pool, then invert
        for (int band_idx = 0; band_idx < n_band; band_idx++)
        {
            Parallel_Reduce::reduce_pool(norms[band_idx]);
            norms[band_idx] = 1.0 / sqrt(norms[band_idx]);
        }

        // Phase 3: parallel normalization of grad/hgrad and accumulation of epsilons
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int band_idx = 0; band_idx < n_band; band_idx++)
        {
            Real eps_0 = 0.0, eps_1 = 0.0, eps_2 = 0.0;
            Real norm = norms[band_idx];
            for (int basis_idx = 0; basis_idx < n_basis; basis_idx++)
            {
                auto item = band_idx * n_basis_max + basis_idx;
                grad_out[item] *= norm;
                hgrad_out[item] *= norm;
                eps_0 += std::real(hpsi_out[item] * std::conj(psi_out[item]));
                eps_1 += std::real(grad_out[item] * std::conj(hpsi_out[item]));
                eps_2 += std::real(grad_out[item] * std::conj(hgrad_out[item]));
            }
            epsilo_0_arr[band_idx] = eps_0;
            epsilo_1_arr[band_idx] = eps_1;
            epsilo_2_arr[band_idx] = eps_2;
        }

        // Phase 4: MPI reduction of epsilons across pool
        for (int band_idx = 0; band_idx < n_band; band_idx++)
        {
            Parallel_Reduce::reduce_pool(epsilo_0_arr[band_idx]);
            Parallel_Reduce::reduce_pool(epsilo_1_arr[band_idx]);
            Parallel_Reduce::reduce_pool(epsilo_2_arr[band_idx]);
        }

        // Phase 5: parallel application of rotation to psi and hpsi
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int band_idx = 0; band_idx < n_band; band_idx++)
        {
            Real epsilo_0 = epsilo_0_arr[band_idx];
            Real epsilo_1 = epsilo_1_arr[band_idx];
            Real epsilo_2 = epsilo_2_arr[band_idx];
            Real theta = 0.5 * std::abs(std::atan(2 * epsilo_1 / (epsilo_0 - epsilo_2)));
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
        // OpenMP parallelization over bands: each band accesses a disjoint memory region
        // [band_idx * n_basis_max, (band_idx+1) * n_basis_max), so no data races occur.
        // MPI collective calls (reduce_pool) use per-band scalar reductions executed serially
        // outside parallel regions to avoid thread-safety issues with MPI.
        // Static scheduling is used since each band has equal workload (n_basis).
        std::vector<Real> norms(n_band);
        std::vector<Real> epsilo_arr(n_band);
        std::vector<Real> err_arr(n_band);
        std::vector<Real> beta_arr(n_band);

        // Phase 1: parallel computation of per-band norms via BLAS dot
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int band_idx = 0; band_idx < n_band; band_idx++)
        {
            auto A = reinterpret_cast<const Real*>(psi_out + band_idx * n_basis_max);
            norms[band_idx] = BlasConnector::dot(2 * n_basis, A, 1, A, 1);
        }

        // Phase 2: MPI reduction of norms across pool, then invert
        for (int band_idx = 0; band_idx < n_band; band_idx++)
        {
            Parallel_Reduce::reduce_pool(norms[band_idx]);
            norms[band_idx] = 1.0 / sqrt(norms[band_idx]);
        }

        // Phase 3: parallel normalization of psi/hpsi and accumulation of epsilo
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int band_idx = 0; band_idx < n_band; band_idx++)
        {
            Real eps = 0.0;
            Real norm = norms[band_idx];
            for (int basis_idx = 0; basis_idx < n_basis; basis_idx++)
            {
                auto item = band_idx * n_basis_max + basis_idx;
                psi_out[item] *= norm;
                hpsi_out[item] *= norm;
                eps += std::real(hpsi_out[item] * std::conj(psi_out[item]));
            }
            epsilo_arr[band_idx] = eps;
        }

        // Phase 4: MPI reduction of epsilons across pool
        for (int band_idx = 0; band_idx < n_band; band_idx++)
        {
            Parallel_Reduce::reduce_pool(epsilo_arr[band_idx]);
        }

        // Phase 5: parallel computation of per-band err and beta
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int band_idx = 0; band_idx < n_band; band_idx++)
        {
            Real err = 0.0;
            Real beta = 0.0;
            Real epsilo = epsilo_arr[band_idx];
            for (int basis_idx = 0; basis_idx < n_basis; basis_idx++)
            {
                auto item = band_idx * n_basis_max + basis_idx;
                T grad_1 = hpsi_out[item] - epsilo * psi_out[item];
                Real grad_2 = std::norm(grad_1);
                err += grad_2;
                beta += grad_2 / prec_in[basis_idx];
            }
            err_arr[band_idx] = err;
            beta_arr[band_idx] = beta;
        }

        // Phase 6: MPI reduction of err and beta across pool
        for (int band_idx = 0; band_idx < n_band; band_idx++)
        {
            Parallel_Reduce::reduce_pool(err_arr[band_idx]);
            Parallel_Reduce::reduce_pool(beta_arr[band_idx]);
        }

        // Phase 7: parallel update of gradient and write output arrays
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int band_idx = 0; band_idx < n_band; band_idx++)
        {
            Real epsilo = epsilo_arr[band_idx];
            Real beta = beta_arr[band_idx];
            for (int basis_idx = 0; basis_idx < n_basis; basis_idx++)
            {
                auto item = band_idx * n_basis_max + basis_idx;
                T grad_1 = hpsi_out[item] - epsilo * psi_out[item];
                grad_out[item] = -grad_1 / prec_in[basis_idx] + beta / beta_out[band_idx] * grad_old_out[item];
            }
            beta_out[band_idx] = beta;
            err_out[band_idx] = sqrt(err_arr[band_idx]);
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
#pragma omp parallel for collapse(1) schedule(static)
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