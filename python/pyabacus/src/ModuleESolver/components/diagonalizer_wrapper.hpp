/**
 * @file diagonalizer_wrapper.hpp
 * @brief Wrapper implementation for diagonalizers
 *
 * Wraps the ABACUS diagonalization solvers to implement IDiagonalizer interface.
 */

#ifndef PYABACUS_ESOLVER_DIAGONALIZER_WRAPPER_HPP
#define PYABACUS_ESOLVER_DIAGONALIZER_WRAPPER_HPP

#include "../interfaces/i_diagonalizer.hpp"
#include "../../utils/pybind_utils.h"
#include "../../hsolver/diago_adapter.hpp"

#include <memory>

namespace pyabacus {
namespace esolver {

/**
 * @brief Wrapper for diagonalization solvers
 *
 * This class wraps the various ABACUS diagonalization methods
 * to provide a unified Python interface.
 *
 * @tparam TK Type for k-space quantities
 */
template <typename TK>
class DiagonalizerWrapper : public IDiagonalizer<TK>
{
public:
    DiagonalizerWrapper() = default;

    DiagonalizerWrapper(int nbasis, int nbands)
        : nbasis_(nbasis), nbands_(nbands)
    {
    }

    ~DiagonalizerWrapper() override = default;

    // ==================== Direct Diagonalization ====================

    DiagResult<TK> diagonalize(int ik,
                               const py::array_t<TK>& Hk,
                               const py::array_t<TK>& Sk,
                               const py::array_t<TK>& psi_init) override
    {
        DiagResult<TK> result;

        // For direct diagonalization, we would use LAPACK/ScaLAPACK
        // This is a placeholder - actual implementation would call
        // the appropriate ABACUS solver

        // Create hpsi and spsi functions from matrices
        auto hpsi_func = [&Hk, this](const py::array_t<TK>& psi) -> py::array_t<TK> {
            return matrix_vector_multiply(Hk, psi);
        };

        auto spsi_func = [&Sk, this](const py::array_t<TK>& psi) -> py::array_t<TK> {
            return matrix_vector_multiply(Sk, psi);
        };

        // Use iterative method as fallback
        py::array_t<double> precond = compute_preconditioner(Hk);
        return diagonalize_iterative(ik, hpsi_func, spsi_func, psi_init, precond);
    }

    // ==================== Iterative Diagonalization ====================

    DiagResult<TK> diagonalize_iterative(
        int ik,
        std::function<py::array_t<TK>(const py::array_t<TK>&)> hpsi_func,
        std::function<py::array_t<TK>(const py::array_t<TK>&)> spsi_func,
        const py::array_t<TK>& psi_init,
        const py::array_t<double>& precond) override
    {
        DiagResult<TK> result;

        switch (config_.method)
        {
            case DiagMethod::Davidson:
                result = diagonalize_davidson(hpsi_func, psi_init, precond);
                break;
            case DiagMethod::DavSubspace:
                result = diagonalize_dav_subspace(hpsi_func, psi_init, precond);
                break;
            case DiagMethod::CG:
                result = diagonalize_cg(hpsi_func, psi_init, precond);
                break;
            default:
                result = diagonalize_davidson(hpsi_func, psi_init, precond);
        }

        return result;
    }

    // ==================== Configuration ====================

    void set_config(const DiagConfig& config) override
    {
        config_ = config;
    }

    DiagConfig get_config() const override { return config_; }

    void set_tolerance(double tol) override
    {
        config_.tolerance = tol;
    }

    void set_max_iterations(int max_iter) override
    {
        config_.max_iterations = max_iter;
    }

    // ==================== Dimension Queries ====================

    int get_nbasis() const override { return nbasis_; }

    int get_nbands() const override { return nbands_; }

    void set_nbands(int nbands) override { nbands_ = nbands; }

    void set_nbasis(int nbasis) { nbasis_ = nbasis; }

private:
    // Davidson diagonalization
    DiagResult<TK> diagonalize_davidson(
        std::function<py::array_t<TK>(const py::array_t<TK>&)> hpsi_func,
        const py::array_t<TK>& psi_init,
        const py::array_t<double>& precond)
    {
        DiagResult<TK> result;

        // Create Davidson adapter
        hsolver::PyDiagoDavidAdapter david(nbasis_, nbands_);

        // Set initial psi
        david.set_psi(psi_init);
        david.init_eigenvalue();

        // Convert preconditioner to vector
        std::vector<double> precond_vec(precond.data(), precond.data() + precond.size());

        // Create diag_ethr vector
        std::vector<double> diag_ethr(nbands_, config_.tolerance);

        // Create comm_info (single process for now)
        ::hsolver::diag_comm_info comm_info(0, 1);

        // Run diagonalization
        int niter = david.diag(
            hpsi_func,
            precond_vec,
            config_.dav_ndim,
            config_.tolerance,
            diag_ethr,
            config_.max_iterations,
            config_.use_paw,
            comm_info
        );

        // Get results
        result.psi = david.get_psi();
        result.eigenvalues = david.get_eigenvalue();
        result.iterations = niter;
        result.converged = (niter < config_.max_iterations);

        return result;
    }

    // Davidson-Subspace diagonalization
    DiagResult<TK> diagonalize_dav_subspace(
        std::function<py::array_t<TK>(const py::array_t<TK>&)> hpsi_func,
        const py::array_t<TK>& psi_init,
        const py::array_t<double>& precond)
    {
        DiagResult<TK> result;

        // Create DavSubspace adapter
        hsolver::PyDiagoDavSubspaceAdapter dav_sub(nbasis_, nbands_);

        // Set initial psi
        dav_sub.set_psi(psi_init);
        dav_sub.init_eigenvalue();

        // Convert preconditioner to vector
        std::vector<double> precond_vec(precond.data(), precond.data() + precond.size());

        // Create diag_ethr vector
        std::vector<double> diag_ethr(nbands_, config_.tolerance);

        // Create comm_info
        ::hsolver::diag_comm_info comm_info(0, 1);

        // Run diagonalization
        int niter = dav_sub.diag(
            hpsi_func,
            precond_vec,
            config_.dav_ndim,
            config_.tolerance,
            config_.max_iterations,
            false,  // need_subspace
            diag_ethr,
            true,   // scf_type
            comm_info,
            0,      // diag_subspace (LAPACK)
            1       // nb2d
        );

        // Get results
        result.psi = dav_sub.get_psi();
        result.eigenvalues = dav_sub.get_eigenvalue();
        result.iterations = niter;
        result.converged = (niter < config_.max_iterations);

        return result;
    }

    // CG diagonalization
    DiagResult<TK> diagonalize_cg(
        std::function<py::array_t<TK>(const py::array_t<TK>&)> hpsi_func,
        const py::array_t<TK>& psi_init,
        const py::array_t<double>& precond)
    {
        DiagResult<TK> result;

#ifdef __ENABLE_ATEN
        // Create CG adapter
        hsolver::PyDiagoCGAdapter cg(nbasis_, nbands_);

        // Set initial psi and preconditioner
        cg.set_psi(psi_init);
        cg.init_eig();
        cg.set_prec(precond);

        // Create diag_ethr vector
        std::vector<double> diag_ethr(nbands_, config_.tolerance);

        // Run diagonalization
        cg.diag(
            hpsi_func,
            config_.dav_ndim,
            config_.tolerance,
            diag_ethr,
            false,  // need_subspace
            true,   // scf_type
            config_.nproc_in_pool
        );

        // Get results
        result.psi = cg.get_psi();
        result.eigenvalues = cg.get_eig();
        result.converged = true;
#else
        // Fall back to Davidson if ATen not available
        result = diagonalize_davidson(hpsi_func, psi_init, precond);
#endif

        return result;
    }

    // Helper: Matrix-vector multiplication
    py::array_t<TK> matrix_vector_multiply(const py::array_t<TK>& matrix,
                                           const py::array_t<TK>& vec)
    {
        using namespace pyabacus::utils;

        auto mat_buf = matrix.request();
        auto vec_buf = vec.request();

        if (mat_buf.ndim != 2)
        {
            throw std::runtime_error("Matrix must be 2D");
        }

        const ssize_t nrow = mat_buf.shape[0];
        const ssize_t ncol = mat_buf.shape[1];
        const ssize_t nvec = (vec_buf.ndim == 1) ? 1 : vec_buf.shape[1];

        py::array_t<TK> result({nrow, nvec});
        auto res_buf = result.request();

        const TK* mat_ptr = static_cast<const TK*>(mat_buf.ptr);
        const TK* vec_ptr = static_cast<const TK*>(vec_buf.ptr);
        TK* res_ptr = static_cast<TK*>(res_buf.ptr);

        // Simple matrix-vector multiplication
        for (ssize_t i = 0; i < nrow; ++i)
        {
            for (ssize_t v = 0; v < nvec; ++v)
            {
                TK sum = TK(0);
                for (ssize_t j = 0; j < ncol; ++j)
                {
                    sum += mat_ptr[i * ncol + j] * vec_ptr[j * nvec + v];
                }
                res_ptr[i * nvec + v] = sum;
            }
        }

        return result;
    }

    // Helper: Compute diagonal preconditioner from Hamiltonian
    py::array_t<double> compute_preconditioner(const py::array_t<TK>& Hk)
    {
        auto buf = Hk.request();
        const ssize_t n = buf.shape[0];
        const TK* ptr = static_cast<const TK*>(buf.ptr);

        py::array_t<double> precond(n);
        double* prec_ptr = precond.mutable_data();

        for (ssize_t i = 0; i < n; ++i)
        {
            // Use diagonal elements as preconditioner
            TK diag = ptr[i * n + i];
            prec_ptr[i] = std::max(std::abs(diag), 1.0);
        }

        return precond;
    }

    int nbasis_ = 0;
    int nbands_ = 0;
    DiagConfig config_;
};

// Type aliases
using DiagonalizerWrapperGamma = DiagonalizerWrapper<double>;
using DiagonalizerWrapperMultiK = DiagonalizerWrapper<std::complex<double>>;

} // namespace esolver
} // namespace pyabacus

#endif // PYABACUS_ESOLVER_DIAGONALIZER_WRAPPER_HPP
