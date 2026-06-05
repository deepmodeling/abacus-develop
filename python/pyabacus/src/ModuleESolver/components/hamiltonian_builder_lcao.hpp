/**
 * @file hamiltonian_builder_lcao.hpp
 * @brief LCAO Hamiltonian builder implementation
 *
 * Wraps the ABACUS LCAO Hamiltonian construction to implement IHamiltonianBuilder.
 */

#ifndef PYABACUS_ESOLVER_HAMILTONIAN_BUILDER_LCAO_HPP
#define PYABACUS_ESOLVER_HAMILTONIAN_BUILDER_LCAO_HPP

#include "../interfaces/i_hamiltonian_builder.hpp"
#include "../../utils/pybind_utils.h"

#include <memory>
#include <vector>
#include <map>
#include <tuple>

namespace pyabacus {
namespace esolver {

/**
 * @brief LCAO Hamiltonian builder implementation
 *
 * This class provides Hamiltonian construction for LCAO calculations,
 * wrapping the ABACUS HamiltLCAO functionality.
 *
 * @tparam TK Type for k-space quantities
 * @tparam TR Type for real-space quantities
 */
template <typename TK, typename TR = double>
class HamiltonianBuilderLCAO : public IHamiltonianBuilder<TK, TR>
{
public:
    HamiltonianBuilderLCAO() = default;

    HamiltonianBuilderLCAO(int nbasis, int nks, int nrow, int ncol)
        : nbasis_(nbasis), nks_(nks), nrow_(nrow), ncol_(ncol)
    {
        // Allocate storage for H(k) and S(k) matrices
        Hk_matrices_.resize(nks);
        Sk_matrices_.resize(nks);
        for (int ik = 0; ik < nks; ++ik)
        {
            Hk_matrices_[ik].resize(nrow * ncol, TK(0));
            Sk_matrices_[ik].resize(nrow * ncol, TK(0));
        }
    }

    ~HamiltonianBuilderLCAO() override = default;

    // ==================== Build/Update ====================

    void build_from_rho(const py::array_t<double>& rho) override
    {
        // Store rho for potential reconstruction
        // In full implementation, this would:
        // 1. Compute Hartree potential from rho
        // 2. Compute XC potential from rho
        // 3. Update H(R) matrices
        // 4. Invalidate H(k) cache

        using namespace pyabacus::utils;
        auto buf = rho.request();

        if (buf.ndim != 2)
        {
            throw std::runtime_error("rho must be 2D array with shape (nspin, nrxx)");
        }

        nspin_ = static_cast<int>(buf.shape[0]);
        nrxx_ = static_cast<int>(buf.shape[1]);

        // Store rho
        rho_data_.assign(static_cast<const double*>(buf.ptr),
                         static_cast<const double*>(buf.ptr) + nspin_ * nrxx_);

        // Mark H(k) as needing update
        hk_valid_.assign(nks_, false);
        valid_ = true;
    }

    void update_Hk(int ik) override
    {
        if (ik < 0 || ik >= nks_)
        {
            throw std::out_of_range("k-point index out of range");
        }

        // In full implementation, this would:
        // 1. Fourier transform H(R) to H(k) for k-point ik
        // 2. Store result in Hk_matrices_[ik]

        hk_valid_[ik] = true;
    }

    void invalidate() override
    {
        hk_valid_.assign(nks_, false);
        valid_ = false;
    }

    // ==================== K-space Matrix Access ====================

    py::array_t<TK> get_Hk(int ik) const override
    {
        validate_kpoint(ik);

        return utils::numpy_from_ptr_copy_2d(
            Hk_matrices_[ik].data(),
            static_cast<size_t>(nrow_),
            static_cast<size_t>(ncol_)
        );
    }

    py::array_t<TK> get_Sk(int ik) const override
    {
        validate_kpoint(ik);

        return utils::numpy_from_ptr_copy_2d(
            Sk_matrices_[ik].data(),
            static_cast<size_t>(nrow_),
            static_cast<size_t>(ncol_)
        );
    }

    // ==================== Real-space Matrix Access ====================

    py::dict get_HR() const override
    {
        py::dict result;
        // In full implementation, would return H(R) data
        // Format: {(iat1, iat2, (R0, R1, R2)): matrix}
        return result;
    }

    py::dict get_SR() const override
    {
        py::dict result;
        // In full implementation, would return S(R) data
        return result;
    }

    // ==================== Matrix-Vector Products ====================

    py::array_t<TK> apply_H(int ik, const py::array_t<TK>& psi_in) const override
    {
        validate_kpoint(ik);

        // H * psi
        return matrix_multiply(Hk_matrices_[ik], psi_in);
    }

    py::array_t<TK> apply_S(int ik, const py::array_t<TK>& psi_in) const override
    {
        validate_kpoint(ik);

        // S * psi
        return matrix_multiply(Sk_matrices_[ik], psi_in);
    }

    // ==================== Dimension Queries ====================

    int get_nbasis() const override { return nbasis_; }

    int get_nks() const override { return nks_; }

    std::pair<int, int> get_local_dims() const override
    {
        return {nrow_, ncol_};
    }

    bool is_valid() const override { return valid_; }

    // ==================== Data Setters (for testing/compatibility) ====================

    void set_Hk_data(int ik, const TK* data, int nrow, int ncol)
    {
        if (ik < 0 || ik >= nks_)
        {
            throw std::out_of_range("k-point index out of range");
        }

        nrow_ = nrow;
        ncol_ = ncol;

        if (static_cast<int>(Hk_matrices_[ik].size()) != nrow * ncol)
        {
            Hk_matrices_[ik].resize(nrow * ncol);
        }

        std::copy(data, data + nrow * ncol, Hk_matrices_[ik].begin());
        hk_valid_[ik] = true;
        valid_ = true;
    }

    void set_Sk_data(int ik, const TK* data, int nrow, int ncol)
    {
        if (ik < 0 || ik >= nks_)
        {
            throw std::out_of_range("k-point index out of range");
        }

        if (static_cast<int>(Sk_matrices_[ik].size()) != nrow * ncol)
        {
            Sk_matrices_[ik].resize(nrow * ncol);
        }

        std::copy(data, data + nrow * ncol, Sk_matrices_[ik].begin());
    }

    void set_dimensions(int nbasis, int nks, int nrow, int ncol)
    {
        nbasis_ = nbasis;
        nks_ = nks;
        nrow_ = nrow;
        ncol_ = ncol;

        Hk_matrices_.resize(nks);
        Sk_matrices_.resize(nks);
        hk_valid_.resize(nks, false);

        for (int ik = 0; ik < nks; ++ik)
        {
            Hk_matrices_[ik].resize(nrow * ncol, TK(0));
            Sk_matrices_[ik].resize(nrow * ncol, TK(0));
        }
    }

private:
    void validate_kpoint(int ik) const
    {
        if (!valid_)
        {
            throw std::runtime_error("Hamiltonian not built. Call build_from_rho first.");
        }
        if (ik < 0 || ik >= nks_)
        {
            throw std::out_of_range("k-point index out of range");
        }
    }

    py::array_t<TK> matrix_multiply(const std::vector<TK>& matrix,
                                    const py::array_t<TK>& vec) const
    {
        auto vec_buf = vec.request();
        const TK* vec_ptr = static_cast<const TK*>(vec_buf.ptr);

        const ssize_t nvec = (vec_buf.ndim == 1) ? 1 :
                             (vec_buf.ndim == 2) ? vec_buf.shape[1] : 1;
        const ssize_t vec_rows = (vec_buf.ndim == 1) ? vec_buf.shape[0] : vec_buf.shape[0];

        if (vec_rows != ncol_)
        {
            throw std::runtime_error("Vector dimension mismatch");
        }

        py::array_t<TK> result({static_cast<ssize_t>(nrow_), nvec});
        auto res_buf = result.request();
        TK* res_ptr = static_cast<TK*>(res_buf.ptr);

        // Matrix-vector multiplication
        for (int i = 0; i < nrow_; ++i)
        {
            for (ssize_t v = 0; v < nvec; ++v)
            {
                TK sum = TK(0);
                for (int j = 0; j < ncol_; ++j)
                {
                    sum += matrix[i * ncol_ + j] * vec_ptr[j * nvec + v];
                }
                res_ptr[i * nvec + v] = sum;
            }
        }

        return result;
    }

    int nbasis_ = 0;
    int nks_ = 0;
    int nrow_ = 0;
    int ncol_ = 0;
    int nspin_ = 1;
    int nrxx_ = 0;
    bool valid_ = false;

    std::vector<std::vector<TK>> Hk_matrices_;
    std::vector<std::vector<TK>> Sk_matrices_;
    std::vector<bool> hk_valid_;
    std::vector<double> rho_data_;
};

// Type aliases
using HamiltonianBuilderLCAOGamma = HamiltonianBuilderLCAO<double, double>;
using HamiltonianBuilderLCAOMultiK = HamiltonianBuilderLCAO<std::complex<double>, double>;

} // namespace esolver
} // namespace pyabacus

#endif // PYABACUS_ESOLVER_HAMILTONIAN_BUILDER_LCAO_HPP
