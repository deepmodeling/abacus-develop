/**
 * @file i_hamiltonian_builder.hpp
 * @brief Abstract interface for Hamiltonian builder
 *
 * Defines the interface for building and accessing Hamiltonian matrices,
 * supporting both k-space and real-space representations.
 */

#ifndef PYABACUS_ESOLVER_I_HAMILTONIAN_BUILDER_HPP
#define PYABACUS_ESOLVER_I_HAMILTONIAN_BUILDER_HPP

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <complex>
#include <functional>

namespace py = pybind11;

namespace pyabacus {
namespace esolver {

/**
 * @brief Abstract interface for Hamiltonian builder
 *
 * This interface defines the contract for building Hamiltonian matrices
 * from charge density and providing access to H(k), S(k), H(R), S(R).
 *
 * @tparam TK Type for k-space quantities (double or complex<double>)
 * @tparam TR Type for real-space quantities (typically double)
 */
template <typename TK, typename TR = double>
class IHamiltonianBuilder
{
public:
    virtual ~IHamiltonianBuilder() = default;

    // ==================== Build/Update ====================

    /**
     * @brief Build Hamiltonian from charge density
     * @param rho Charge density array with shape (nspin, nrxx)
     */
    virtual void build_from_rho(const py::array_t<double>& rho) = 0;

    /**
     * @brief Update H(k) for a specific k-point
     * @param ik K-point index
     */
    virtual void update_Hk(int ik) = 0;

    /**
     * @brief Invalidate cached matrices (force rebuild)
     */
    virtual void invalidate() = 0;

    // ==================== K-space Matrix Access ====================

    /**
     * @brief Get H(k) matrix for specific k-point
     * @param ik K-point index
     * @return Hamiltonian matrix as numpy array
     */
    virtual py::array_t<TK> get_Hk(int ik) const = 0;

    /**
     * @brief Get S(k) overlap matrix for specific k-point
     * @param ik K-point index
     * @return Overlap matrix as numpy array
     */
    virtual py::array_t<TK> get_Sk(int ik) const = 0;

    // ==================== Real-space Matrix Access ====================

    /**
     * @brief Get H(R) in sparse format
     * @return Dictionary mapping (iat1, iat2, R) -> matrix
     */
    virtual py::dict get_HR() const = 0;

    /**
     * @brief Get S(R) in sparse format
     * @return Dictionary mapping (iat1, iat2, R) -> matrix
     */
    virtual py::dict get_SR() const = 0;

    // ==================== Matrix-Vector Products ====================

    /**
     * @brief Apply Hamiltonian to wave function: H|psi>
     * @param ik K-point index
     * @param psi_in Input wave function
     * @return H * psi_in
     */
    virtual py::array_t<TK> apply_H(int ik, const py::array_t<TK>& psi_in) const = 0;

    /**
     * @brief Apply overlap matrix to wave function: S|psi>
     * @param ik K-point index
     * @param psi_in Input wave function
     * @return S * psi_in
     */
    virtual py::array_t<TK> apply_S(int ik, const py::array_t<TK>& psi_in) const = 0;

    // ==================== Dimension Queries ====================

    /**
     * @brief Get number of basis functions
     * @return Number of basis functions
     */
    virtual int get_nbasis() const = 0;

    /**
     * @brief Get number of k-points
     * @return Number of k-points
     */
    virtual int get_nks() const = 0;

    /**
     * @brief Get local matrix dimensions (for 2D distribution)
     * @return Pair of (nrow, ncol)
     */
    virtual std::pair<int, int> get_local_dims() const = 0;

    /**
     * @brief Check if Hamiltonian data is valid
     * @return true if valid
     */
    virtual bool is_valid() const = 0;
};

// Type aliases for common use cases
using IHamiltonianBuilderGamma = IHamiltonianBuilder<double, double>;
using IHamiltonianBuilderMultiK = IHamiltonianBuilder<std::complex<double>, double>;

} // namespace esolver
} // namespace pyabacus

#endif // PYABACUS_ESOLVER_I_HAMILTONIAN_BUILDER_HPP
