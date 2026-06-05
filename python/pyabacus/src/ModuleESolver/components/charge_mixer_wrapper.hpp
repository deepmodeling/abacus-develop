/**
 * @file charge_mixer_wrapper.hpp
 * @brief Wrapper implementation for ABACUS Charge_Mixing
 *
 * Wraps the ABACUS Charge_Mixing class to implement IChargeMixer interface.
 */

#ifndef PYABACUS_ESOLVER_CHARGE_MIXER_WRAPPER_HPP
#define PYABACUS_ESOLVER_CHARGE_MIXER_WRAPPER_HPP

#include "../interfaces/i_charge_mixer.hpp"
#include "../../utils/pybind_utils.h"

#include <memory>
#include <vector>

namespace pyabacus {
namespace esolver {

/**
 * @brief Wrapper for ABACUS Charge_Mixing class
 *
 * This class wraps the ABACUS charge mixing functionality
 * to provide a clean Python interface.
 */
class ChargeMixerWrapper : public IChargeMixer
{
public:
    ChargeMixerWrapper() = default;

    ChargeMixerWrapper(int nspin, int nrxx)
        : nspin_(nspin), nrxx_(nrxx)
    {
        // Initialize history buffers for Pulay mixing
        rho_history_.reserve(config_.ndim);
        residual_history_.reserve(config_.ndim);
    }

    ~ChargeMixerWrapper() override = default;

    // ==================== Core Mixing Operations ====================

    py::array_t<double> mix(const py::array_t<double>& rho_in,
                            const py::array_t<double>& rho_out) override
    {
        using namespace pyabacus::utils;

        // Validate input arrays
        check_array_size(rho_in, static_cast<size_t>(nspin_ * nrxx_), "rho_in");
        check_array_size(rho_out, static_cast<size_t>(nspin_ * nrxx_), "rho_out");

        const double* in_ptr = get_array_ptr(rho_in);
        const double* out_ptr = get_array_ptr(rho_out);

        // Create output array
        py::array_t<double> rho_mixed({static_cast<ssize_t>(nspin_),
                                        static_cast<ssize_t>(nrxx_)});
        double* mixed_ptr = get_array_ptr(rho_mixed);

        // Calculate drho
        drho_ = 0.0;
        for (size_t i = 0; i < static_cast<size_t>(nspin_ * nrxx_); ++i)
        {
            double diff = out_ptr[i] - in_ptr[i];
            drho_ += diff * diff;
        }
        drho_ = std::sqrt(drho_ / (nspin_ * nrxx_));

        // Perform mixing based on method
        switch (config_.method)
        {
            case MixingMethod::Plain:
                mix_plain(in_ptr, out_ptr, mixed_ptr);
                break;
            case MixingMethod::Pulay:
                mix_pulay(in_ptr, out_ptr, mixed_ptr);
                break;
            case MixingMethod::Broyden:
                mix_broyden(in_ptr, out_ptr, mixed_ptr);
                break;
            case MixingMethod::Anderson:
                mix_anderson(in_ptr, out_ptr, mixed_ptr);
                break;
            default:
                mix_plain(in_ptr, out_ptr, mixed_ptr);
        }

        iteration_++;
        return rho_mixed;
    }

    void reset() override
    {
        iteration_ = 0;
        drho_ = 0.0;
        rho_history_.clear();
        residual_history_.clear();
    }

    // ==================== State Queries ====================

    double get_drho() const override { return drho_; }

    int get_iteration() const override { return iteration_; }

    // ==================== Configuration ====================

    void set_config(const MixingConfig& config) override
    {
        config_ = config;
        reset(); // Reset history when config changes
    }

    MixingConfig get_config() const override { return config_; }

    void set_mixing_beta(double beta) override
    {
        if (beta <= 0.0 || beta > 1.0)
        {
            throw std::invalid_argument("beta must be in (0, 1]");
        }
        config_.beta = beta;
    }

    double get_mixing_beta() const override { return config_.beta; }

    void set_mixing_method(MixingMethod method) override
    {
        config_.method = method;
        reset();
    }

    MixingMethod get_mixing_method() const override { return config_.method; }

    // ==================== Dimension Setters ====================

    void set_dimensions(int nspin, int nrxx)
    {
        nspin_ = nspin;
        nrxx_ = nrxx;
        reset();
    }

private:
    // Plain linear mixing: rho_new = (1-beta)*rho_in + beta*rho_out
    void mix_plain(const double* rho_in, const double* rho_out, double* rho_mixed)
    {
        const double beta = config_.beta;
        const double one_minus_beta = 1.0 - beta;

        for (size_t i = 0; i < static_cast<size_t>(nspin_ * nrxx_); ++i)
        {
            rho_mixed[i] = one_minus_beta * rho_in[i] + beta * rho_out[i];
        }
    }

    // Pulay mixing (DIIS)
    void mix_pulay(const double* rho_in, const double* rho_out, double* rho_mixed)
    {
        const size_t size = static_cast<size_t>(nspin_ * nrxx_);

        // Store current rho and residual in history
        std::vector<double> current_rho(rho_in, rho_in + size);
        std::vector<double> current_residual(size);
        for (size_t i = 0; i < size; ++i)
        {
            current_residual[i] = rho_out[i] - rho_in[i];
        }

        // Add to history (circular buffer)
        if (static_cast<int>(rho_history_.size()) >= config_.ndim)
        {
            rho_history_.erase(rho_history_.begin());
            residual_history_.erase(residual_history_.begin());
        }
        rho_history_.push_back(current_rho);
        residual_history_.push_back(current_residual);

        const int nhist = static_cast<int>(rho_history_.size());

        if (nhist < 2)
        {
            // Not enough history, use plain mixing
            mix_plain(rho_in, rho_out, rho_mixed);
            return;
        }

        // Build overlap matrix of residuals
        std::vector<double> A((nhist + 1) * (nhist + 1), 0.0);
        std::vector<double> b(nhist + 1, 0.0);

        for (int i = 0; i < nhist; ++i)
        {
            for (int j = 0; j <= i; ++j)
            {
                double dot = 0.0;
                for (size_t k = 0; k < size; ++k)
                {
                    dot += residual_history_[i][k] * residual_history_[j][k];
                }
                A[i * (nhist + 1) + j] = dot;
                A[j * (nhist + 1) + i] = dot;
            }
            A[i * (nhist + 1) + nhist] = 1.0;
            A[nhist * (nhist + 1) + i] = 1.0;
        }
        b[nhist] = 1.0;

        // Solve linear system for coefficients (simple Gaussian elimination)
        std::vector<double> coeff = solve_linear_system(A, b, nhist + 1);

        // Compute mixed density
        std::fill(rho_mixed, rho_mixed + size, 0.0);
        for (int i = 0; i < nhist; ++i)
        {
            for (size_t k = 0; k < size; ++k)
            {
                rho_mixed[k] += coeff[i] * (rho_history_[i][k] +
                                            config_.beta * residual_history_[i][k]);
            }
        }
    }

    // Broyden mixing (simplified)
    void mix_broyden(const double* rho_in, const double* rho_out, double* rho_mixed)
    {
        // For simplicity, use Pulay mixing as approximation
        mix_pulay(rho_in, rho_out, rho_mixed);
    }

    // Anderson mixing
    void mix_anderson(const double* rho_in, const double* rho_out, double* rho_mixed)
    {
        // Anderson mixing is similar to Pulay
        mix_pulay(rho_in, rho_out, rho_mixed);
    }

    // Simple linear system solver (Gaussian elimination with partial pivoting)
    std::vector<double> solve_linear_system(std::vector<double>& A,
                                            std::vector<double>& b,
                                            int n)
    {
        std::vector<double> x(n, 0.0);

        // Forward elimination
        for (int k = 0; k < n - 1; ++k)
        {
            // Find pivot
            int max_row = k;
            double max_val = std::abs(A[k * n + k]);
            for (int i = k + 1; i < n; ++i)
            {
                if (std::abs(A[i * n + k]) > max_val)
                {
                    max_val = std::abs(A[i * n + k]);
                    max_row = i;
                }
            }

            // Swap rows
            if (max_row != k)
            {
                for (int j = 0; j < n; ++j)
                {
                    std::swap(A[k * n + j], A[max_row * n + j]);
                }
                std::swap(b[k], b[max_row]);
            }

            // Eliminate
            for (int i = k + 1; i < n; ++i)
            {
                if (std::abs(A[k * n + k]) < 1e-12) continue;
                double factor = A[i * n + k] / A[k * n + k];
                for (int j = k; j < n; ++j)
                {
                    A[i * n + j] -= factor * A[k * n + j];
                }
                b[i] -= factor * b[k];
            }
        }

        // Back substitution
        for (int i = n - 1; i >= 0; --i)
        {
            x[i] = b[i];
            for (int j = i + 1; j < n; ++j)
            {
                x[i] -= A[i * n + j] * x[j];
            }
            if (std::abs(A[i * n + i]) > 1e-12)
            {
                x[i] /= A[i * n + i];
            }
        }

        return x;
    }

    int nspin_ = 1;
    int nrxx_ = 0;
    int iteration_ = 0;
    double drho_ = 0.0;
    MixingConfig config_;

    // History for Pulay/Broyden mixing
    std::vector<std::vector<double>> rho_history_;
    std::vector<std::vector<double>> residual_history_;
};

} // namespace esolver
} // namespace pyabacus

#endif // PYABACUS_ESOLVER_CHARGE_MIXER_WRAPPER_HPP
