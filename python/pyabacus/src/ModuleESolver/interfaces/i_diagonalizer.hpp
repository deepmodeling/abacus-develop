/**
 * @file i_diagonalizer.hpp
 * @brief Abstract interface for eigenvalue solver (diagonalizer)
 *
 * Defines the interface for diagonalizing Hamiltonian matrices
 * to obtain eigenvalues and eigenvectors.
 */

#ifndef PYABACUS_ESOLVER_I_DIAGONALIZER_HPP
#define PYABACUS_ESOLVER_I_DIAGONALIZER_HPP

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/functional.h>

#include <complex>
#include <functional>
#include <string>
#include <vector>

namespace py = pybind11;

namespace pyabacus {
namespace esolver {

/**
 * @brief Diagonalization method types
 */
enum class DiagMethod
{
    Davidson,      ///< Davidson iterative method
    DavSubspace,   ///< Davidson with subspace rotation
    CG,            ///< Conjugate gradient
    LAPACK,        ///< Direct LAPACK diagonalization
    ScaLAPACK,     ///< Parallel ScaLAPACK
    ELPA           ///< ELPA eigensolver
};

/**
 * @brief Result of diagonalization
 */
template <typename TK>
struct DiagResult
{
    py::array_t<TK> psi;           ///< Eigenvectors
    py::array_t<double> eigenvalues; ///< Eigenvalues
    int iterations = 0;             ///< Number of iterations (for iterative methods)
    bool converged = false;         ///< Whether converged
    double residual = 0.0;          ///< Final residual
};

/**
 * @brief Configuration for diagonalization
 */
struct DiagConfig
{
    DiagMethod method = DiagMethod::Davidson;
    double tolerance = 1e-6;     ///< Convergence tolerance
    int max_iterations = 100;    ///< Maximum iterations
    int dav_ndim = 4;            ///< Davidson subspace dimension multiplier
    bool use_paw = false;        ///< Use PAW method
    int nproc_in_pool = 1;       ///< Number of processes in pool
};

/**
 * @brief Abstract interface for diagonalizer
 *
 * This interface defines the contract for eigenvalue solvers,
 * supporting both direct and iterative methods.
 *
 * @tparam TK Type for k-space quantities (double or complex<double>)
 */
template <typename TK>
class IDiagonalizer
{
public:
    virtual ~IDiagonalizer() = default;

    // ==================== Direct Diagonalization ====================

    /**
     * @brief Diagonalize H(k) and S(k) matrices directly
     * @param ik K-point index
     * @param Hk Hamiltonian matrix
     * @param Sk Overlap matrix
     * @param psi_init Initial guess for eigenvectors (optional)
     * @return Diagonalization result
     */
    virtual DiagResult<TK> diagonalize(int ik,
                                       const py::array_t<TK>& Hk,
                                       const py::array_t<TK>& Sk,
                                       const py::array_t<TK>& psi_init) = 0;

    // ==================== Iterative Diagonalization ====================

    /**
     * @brief Diagonalize using matrix-vector product functions
     *
     * For iterative methods that don't need explicit matrices.
     *
     * @param ik K-point index
     * @param hpsi_func Function computing H|psi>
     * @param spsi_func Function computing S|psi>
     * @param psi_init Initial guess for eigenvectors
     * @param precond Preconditioner (diagonal approximation to H)
     * @return Diagonalization result
     */
    virtual DiagResult<TK> diagonalize_iterative(
        int ik,
        std::function<py::array_t<TK>(const py::array_t<TK>&)> hpsi_func,
        std::function<py::array_t<TK>(const py::array_t<TK>&)> spsi_func,
        const py::array_t<TK>& psi_init,
        const py::array_t<double>& precond) = 0;

    // ==================== Configuration ====================

    /**
     * @brief Set diagonalization configuration
     * @param config Configuration
     */
    virtual void set_config(const DiagConfig& config) = 0;

    /**
     * @brief Get current configuration
     * @return Current configuration
     */
    virtual DiagConfig get_config() const = 0;

    /**
     * @brief Set convergence tolerance
     * @param tol Tolerance
     */
    virtual void set_tolerance(double tol) = 0;

    /**
     * @brief Set maximum iterations
     * @param max_iter Maximum iterations
     */
    virtual void set_max_iterations(int max_iter) = 0;

    // ==================== Dimension Queries ====================

    /**
     * @brief Get number of basis functions
     * @return Number of basis functions
     */
    virtual int get_nbasis() const = 0;

    /**
     * @brief Get number of bands to compute
     * @return Number of bands
     */
    virtual int get_nbands() const = 0;

    /**
     * @brief Set number of bands to compute
     * @param nbands Number of bands
     */
    virtual void set_nbands(int nbands) = 0;
};

// Type aliases for common use cases
using IDiagonalizerGamma = IDiagonalizer<double>;
using IDiagonalizerMultiK = IDiagonalizer<std::complex<double>>;

/**
 * @brief Convert DiagMethod enum to string
 */
inline std::string diag_method_to_string(DiagMethod method)
{
    switch (method)
    {
        case DiagMethod::Davidson: return "davidson";
        case DiagMethod::DavSubspace: return "dav_subspace";
        case DiagMethod::CG: return "cg";
        case DiagMethod::LAPACK: return "lapack";
        case DiagMethod::ScaLAPACK: return "scalapack";
        case DiagMethod::ELPA: return "elpa";
        default: return "unknown";
    }
}

/**
 * @brief Convert string to DiagMethod enum
 */
inline DiagMethod string_to_diag_method(const std::string& str)
{
    if (str == "davidson") return DiagMethod::Davidson;
    if (str == "dav_subspace") return DiagMethod::DavSubspace;
    if (str == "cg") return DiagMethod::CG;
    if (str == "lapack") return DiagMethod::LAPACK;
    if (str == "scalapack") return DiagMethod::ScaLAPACK;
    if (str == "elpa") return DiagMethod::ELPA;
    return DiagMethod::Davidson; // default
}

} // namespace esolver
} // namespace pyabacus

#endif // PYABACUS_ESOLVER_I_DIAGONALIZER_HPP
