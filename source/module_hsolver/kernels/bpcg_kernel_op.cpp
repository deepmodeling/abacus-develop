#include "module_hsolver/kernels/bpcg_kernel_op.h"
#include "module_base/blas_connector.h"
#include "module_base/kernels/math_kernel_op.h"
#include "module_base/parallel_reduce.h"
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
        for (int band_idx = 0; band_idx < n_band; band_idx++)
        {
            Real epsilo_0 = 0.0, epsilo_1 = 0.0, epsilo_2 = 0.0;
            Real theta = 0.0, cos_theta = 0.0, sin_theta = 0.0;
            auto A = reinterpret_cast<const Real*>(grad_out + band_idx * n_basis_max);
            Real norm = BlasConnector::dot(2 * n_basis, A, 1, A, 1);
            Parallel_Reduce::reduce_pool(norm);
            norm = 1.0 / sqrt(norm);
            for (int basis_idx = 0; basis_idx < n_basis; basis_idx++)
            {
                auto item = band_idx * n_basis_max + basis_idx;
                grad_out[item] *= norm;
                hgrad_out[item] *= norm;
                epsilo_0 += std::real(hpsi_out[item] * std::conj(psi_out[item]));
                epsilo_1 += std::real(grad_out[item] * std::conj(hpsi_out[item]));
                epsilo_2 += std::real(grad_out[item] * std::conj(hgrad_out[item]));
            }
            Parallel_Reduce::reduce_pool(epsilo_0);
            Parallel_Reduce::reduce_pool(epsilo_1);
            Parallel_Reduce::reduce_pool(epsilo_2);
            theta = 0.5 * std::abs(std::atan(2 * epsilo_1 / (epsilo_0 - epsilo_2)));
            cos_theta = std::cos(theta);
            sin_theta = std::sin(theta);
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
        for (int band_idx = 0; band_idx < n_band; band_idx++)
        {
            Real err = 0.0;
            Real beta = 0.0;
            Real epsilo = 0.0;
            Real grad_2 = {0.0};
            T grad_1 = {0.0, 0.0};
            auto A = reinterpret_cast<const Real*>(psi_out + band_idx * n_basis_max);
            Real norm = BlasConnector::dot(2 * n_basis, A, 1, A, 1);
            Parallel_Reduce::reduce_pool(norm);
            norm = 1.0 / sqrt(norm);
            for (int basis_idx = 0; basis_idx < n_basis; basis_idx++)
            {
                auto item = band_idx * n_basis_max + basis_idx;
                psi_out[item] *= norm;
                hpsi_out[item] *= norm;
                epsilo += std::real(hpsi_out[item] * std::conj(psi_out[item]));
            }
            Parallel_Reduce::reduce_pool(epsilo);
            for (int basis_idx = 0; basis_idx < n_basis; basis_idx++)
            {
                auto item = band_idx * n_basis_max + basis_idx;
                grad_1 = hpsi_out[item] - epsilo * psi_out[item];
                grad_2 = std::norm(grad_1);
                err += grad_2;
                beta += grad_2 / prec_in[basis_idx]; /// Mark here as we should div the prec?
            }
            Parallel_Reduce::reduce_pool(err);
            Parallel_Reduce::reduce_pool(beta);
            for (int basis_idx = 0; basis_idx < n_basis; basis_idx++)
            {
                auto item = band_idx * n_basis_max + basis_idx;
                grad_1 = hpsi_out[item] - epsilo * psi_out[item];
                grad_out[item] = -grad_1 / prec_in[basis_idx] + beta / beta_out[band_idx] * grad_old_out[item];
            }
            beta_out[band_idx] = beta;
            err_out[band_idx] = sqrt(err);
        }
    }
};

template struct calc_grad_with_block_op<std::complex<float>, base_device::DEVICE_CPU>;
template struct line_minimize_with_block_op<std::complex<float>, base_device::DEVICE_CPU>;
template struct calc_grad_with_block_op<std::complex<double>, base_device::DEVICE_CPU>;
template struct line_minimize_with_block_op<std::complex<double>, base_device::DEVICE_CPU>;
} // namespace hsolver