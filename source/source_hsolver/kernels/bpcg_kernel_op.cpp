#include "source_hsolver/kernels/bpcg_kernel_op.h"
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
                    T* sgrad_out,
                    T* psi_out,
                    T* hpsi_out,
                    T* spsi_out,
                    const int& n_basis,
                    const int& n_basis_max,
                    const int& n_band)
    {
        for (int band_idx = 0; band_idx < n_band; ++band_idx)
        {
            Real norm = 0.0;
            Real epsilo_0 = 0.0;
            Real epsilo_1 = 0.0;
            Real epsilo_2 = 0.0;
            for (int basis_idx = 0; basis_idx < n_basis; ++basis_idx)
            {
                const int item = band_idx * n_basis_max + basis_idx;
                norm += std::real(sgrad_out[item] * std::conj(grad_out[item]));
            }
            Parallel_Reduce::reduce_pool(norm);
            if (!(norm > 1.0e-20))
            {
                continue;
            }
            norm = 1.0 / std::sqrt(norm);
            for (int basis_idx = 0; basis_idx < n_basis; ++basis_idx)
            {
                const int item = band_idx * n_basis_max + basis_idx;
                grad_out[item] *= norm;
                hgrad_out[item] *= norm;
                sgrad_out[item] *= norm;
                epsilo_0 += std::real(hpsi_out[item] * std::conj(psi_out[item]));
                epsilo_1 += std::real(grad_out[item] * std::conj(hpsi_out[item]));
                epsilo_2 += std::real(grad_out[item] * std::conj(hgrad_out[item]));
            }
            Parallel_Reduce::reduce_pool(epsilo_0);
            Parallel_Reduce::reduce_pool(epsilo_1);
            Parallel_Reduce::reduce_pool(epsilo_2);
            Real theta = 0.5 * std::atan2(2 * epsilo_1, epsilo_0 - epsilo_2);
            // Choose the rotation associated with the lower Ritz value.
            const Real energy_delta = (epsilo_0 - epsilo_2) * std::cos(2.0 * theta)
                                    + 2.0 * epsilo_1 * std::sin(2.0 * theta);
            const Real energy_1 = 0.5 * (epsilo_0 + epsilo_2 + energy_delta);
            const Real energy_2 = 0.5 * (epsilo_0 + epsilo_2 - energy_delta);
            if (energy_1 > energy_2)
            {
                theta += 2.0 * std::atan(1.0);
            }
            const Real cos_theta = std::cos(theta);
            const Real sin_theta = std::sin(theta);
            for (int basis_idx = 0; basis_idx < n_basis; ++basis_idx)
            {
                const int item = band_idx * n_basis_max + basis_idx;
                psi_out[item] = psi_out[item] * cos_theta + grad_out[item] * sin_theta;
                hpsi_out[item] = hpsi_out[item] * cos_theta + hgrad_out[item] * sin_theta;
                spsi_out[item] = spsi_out[item] * cos_theta + sgrad_out[item] * sin_theta;
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
                    T* spsi_out,
                    T* grad_out,
                    T* grad_old_out,
                    const int& n_basis,
                    const int& n_basis_max,
                    const int& n_band)
    {
        for (int band_idx = 0; band_idx < n_band; ++band_idx)
        {
            Real norm = 0.0;
            for (int basis_idx = 0; basis_idx < n_basis; ++basis_idx)
            {
                const int item = band_idx * n_basis_max + basis_idx;
                norm += std::real(spsi_out[item] * std::conj(psi_out[item]));
            }
            Parallel_Reduce::reduce_pool(norm);
            norm = 1.0 / std::sqrt(norm);

            Real epsilo = 0.0;
            for (int basis_idx = 0; basis_idx < n_basis; ++basis_idx)
            {
                const int item = band_idx * n_basis_max + basis_idx;
                psi_out[item] *= norm;
                hpsi_out[item] *= norm;
                spsi_out[item] *= norm;
                epsilo += std::real(hpsi_out[item] * std::conj(psi_out[item]));
            }
            Parallel_Reduce::reduce_pool(epsilo);

            Real err = 0.0;
            Real beta = 0.0;
            for (int basis_idx = 0; basis_idx < n_basis; ++basis_idx)
            {
                const int item = band_idx * n_basis_max + basis_idx;
                const T residual = hpsi_out[item] - epsilo * spsi_out[item];
                const Real residual_norm = std::norm(residual);
                err += residual_norm;
                beta += residual_norm / prec_in[basis_idx];
            }
            Parallel_Reduce::reduce_pool(err);
            Parallel_Reduce::reduce_pool(beta);
            for (int basis_idx = 0; basis_idx < n_basis; ++basis_idx)
            {
                const int item = band_idx * n_basis_max + basis_idx;
                const T residual = hpsi_out[item] - epsilo * spsi_out[item];
                grad_out[item] = -residual / prec_in[basis_idx] + beta / beta_out[band_idx] * grad_old_out[item];
            }
            beta_out[band_idx] = beta;
            err_out[band_idx] = std::sqrt(err);
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
            ModuleBase::vector_div_vector_op<T, base_device::DEVICE_CPU, Real>()(
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
