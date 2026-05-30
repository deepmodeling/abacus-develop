/**
 * @file diago_traits.hpp
 * @brief Type traits and storage policies for diagonalizer adapters
 *
 * This header defines:
 * - Storage policies for different memory management strategies
 * - Solver traits for different diagonalization algorithms
 */

#ifndef PYABACUS_HSOLVER_DIAGO_TRAITS_HPP
#define PYABACUS_HSOLVER_DIAGO_TRAITS_HPP

#include <complex>
#include <memory>
#include <vector>
#include <functional>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>

#include "source_hsolver/diago_david.h"
#include "source_hsolver/diago_dav_subspace.h"
#include "source_hsolver/diago_cg.h"
#include "source_base/module_device/memory_op.h"

#ifdef __ENABLE_ATEN
#include <ATen/core/tensor.h>
#include <ATen/core/tensor_map.h>
#include <ATen/core/tensor_types.h>
#endif

namespace py = pybind11;

namespace pyabacus {
namespace hsolver {

// ============================================================================
// Storage Policies
// ============================================================================

/**
 * @brief Storage policy using raw pointers
 *
 * Used by DiagoDavid and DiagoDavSubspace which work with raw arrays.
 */
template <typename T>
class RawPointerStorage
{
public:
    using value_type = T;

    RawPointerStorage() = default;

    void allocate(int nbasis, int nband)
    {
        nbasis_ = nbasis;
        nband_ = nband;
        psi_ = new T[nbasis * nband];
        eigenvalue_ = new double[nband];
    }

    ~RawPointerStorage()
    {
        cleanup();
    }

    // Move semantics
    RawPointerStorage(RawPointerStorage&& other) noexcept
        : psi_(other.psi_)
        , eigenvalue_(other.eigenvalue_)
        , nbasis_(other.nbasis_)
        , nband_(other.nband_)
    {
        other.psi_ = nullptr;
        other.eigenvalue_ = nullptr;
    }

    RawPointerStorage& operator=(RawPointerStorage&& other) noexcept
    {
        if (this != &other)
        {
            cleanup();
            psi_ = other.psi_;
            eigenvalue_ = other.eigenvalue_;
            nbasis_ = other.nbasis_;
            nband_ = other.nband_;
            other.psi_ = nullptr;
            other.eigenvalue_ = nullptr;
        }
        return *this;
    }

    // Disable copy
    RawPointerStorage(const RawPointerStorage&) = delete;
    RawPointerStorage& operator=(const RawPointerStorage&) = delete;

    void set_psi(py::array_t<T> psi_in)
    {
        if (static_cast<size_t>(psi_in.size()) != static_cast<size_t>(nbasis_ * nband_))
        {
            throw std::runtime_error("psi_in size mismatch");
        }
        for (size_t i = 0; i < static_cast<size_t>(nbasis_ * nband_); ++i)
        {
            psi_[i] = psi_in.at(i);
        }
    }

    py::array_t<T> get_psi() const
    {
        py::array_t<T> psi_out(nband_ * nbasis_);
        py::buffer_info buf = psi_out.request();
        T* ptr = static_cast<T*>(buf.ptr);
        for (size_t i = 0; i < static_cast<size_t>(nband_ * nbasis_); ++i)
        {
            ptr[i] = psi_[i];
        }
        return psi_out;
    }

    void init_eigenvalue()
    {
        for (int i = 0; i < nband_; ++i)
        {
            eigenvalue_[i] = 0.0;
        }
    }

    py::array_t<double> get_eigenvalue() const
    {
        py::array_t<double> eig_out(nband_);
        py::buffer_info buf = eig_out.request();
        double* ptr = static_cast<double*>(buf.ptr);
        for (int i = 0; i < nband_; ++i)
        {
            ptr[i] = eigenvalue_[i];
        }
        return eig_out;
    }

    // Accessors for solver
    T* psi_ptr() { return psi_; }
    double* eigenvalue_ptr() { return eigenvalue_; }
    int nbasis() const { return nbasis_; }
    int nband() const { return nband_; }

private:
    void cleanup()
    {
        if (psi_ != nullptr)
        {
            delete[] psi_;
            psi_ = nullptr;
        }
        if (eigenvalue_ != nullptr)
        {
            delete[] eigenvalue_;
            eigenvalue_ = nullptr;
        }
    }

    T* psi_ = nullptr;
    double* eigenvalue_ = nullptr;
    int nbasis_ = 0;
    int nband_ = 0;
};

#ifdef __ENABLE_ATEN
/**
 * @brief Storage policy using ATen tensors
 *
 * Used by DiagoCG which works with ATen tensor interface.
 */
template <typename T>
class TensorStorage
{
public:
    using value_type = T;

    TensorStorage() = default;

    void allocate(int dim, int num_eigs)
    {
        dim_ = dim;
        num_eigs_ = num_eigs;
        // Tensors are allocated lazily
    }

    ~TensorStorage()
    {
        cleanup();
    }

    // Move semantics
    TensorStorage(TensorStorage&& other) noexcept
        : psi_(other.psi_)
        , eig_(other.eig_)
        , prec_(other.prec_)
        , dim_(other.dim_)
        , num_eigs_(other.num_eigs_)
    {
        other.psi_ = nullptr;
        other.eig_ = nullptr;
        other.prec_ = nullptr;
    }

    TensorStorage& operator=(TensorStorage&& other) noexcept
    {
        if (this != &other)
        {
            cleanup();
            psi_ = other.psi_;
            eig_ = other.eig_;
            prec_ = other.prec_;
            dim_ = other.dim_;
            num_eigs_ = other.num_eigs_;
            other.psi_ = nullptr;
            other.eig_ = nullptr;
            other.prec_ = nullptr;
        }
        return *this;
    }

    // Disable copy
    TensorStorage(const TensorStorage&) = delete;
    TensorStorage& operator=(const TensorStorage&) = delete;

    void set_psi(py::array_t<T> psi_in)
    {
        py::buffer_info buf = psi_in.request();
        T* ptr = static_cast<T*>(buf.ptr);

        psi_ = new ct::TensorMap(
            ptr,
            ct::DataType::DT_COMPLEX_DOUBLE,
            ct::DeviceType::CpuDevice,
            ct::TensorShape({num_eigs_, dim_})
        );
    }

    py::array_t<T> get_psi() const
    {
        if (psi_ == nullptr)
        {
            throw std::runtime_error("psi is not initialized");
        }
        py::array_t<T> psi_out({num_eigs_, dim_});
        py::buffer_info buf = psi_out.request();
        T* ptr = static_cast<T*>(buf.ptr);
        T* psi_ptr = psi_->data<T>();
        std::copy(psi_ptr, psi_ptr + psi_->NumElements(), ptr);
        return psi_out;
    }

    void init_eigenvalue()
    {
        eig_ = new ct::Tensor(ct::DataType::DT_DOUBLE, {num_eigs_});
        eig_->zero();
    }

    py::array_t<double> get_eigenvalue() const
    {
        if (eig_ == nullptr)
        {
            throw std::runtime_error("eigenvalue is not initialized");
        }
        py::array_t<double> eig_out(eig_->NumElements());
        py::buffer_info buf = eig_out.request();
        double* ptr = static_cast<double*>(buf.ptr);
        double* eig_ptr = eig_->data<double>();
        std::copy(eig_ptr, eig_ptr + eig_->NumElements(), ptr);
        return eig_out;
    }

    void set_preconditioner(py::array_t<double> prec_in)
    {
        py::buffer_info buf = prec_in.request();
        double* ptr = static_cast<double*>(buf.ptr);

        prec_ = new ct::TensorMap(
            ptr,
            ct::DataType::DT_DOUBLE,
            ct::DeviceType::CpuDevice,
            ct::TensorShape({dim_})
        );
    }

    // Accessors for solver
    ct::Tensor* psi_tensor() { return psi_; }
    ct::Tensor* eig_tensor() { return eig_; }
    ct::Tensor* prec_tensor() { return prec_; }
    int dim() const { return dim_; }
    int num_eigs() const { return num_eigs_; }

private:
    void cleanup()
    {
        if (psi_ != nullptr)
        {
            delete psi_;
            psi_ = nullptr;
        }
        if (eig_ != nullptr)
        {
            delete eig_;
            eig_ = nullptr;
        }
        if (prec_ != nullptr)
        {
            delete prec_;
            prec_ = nullptr;
        }
    }

    ct::Tensor* psi_ = nullptr;
    ct::Tensor* eig_ = nullptr;
    ct::Tensor* prec_ = nullptr;
    int dim_ = 0;
    int num_eigs_ = 0;
};
#endif // __ENABLE_ATEN

// ============================================================================
// Solver Traits
// ============================================================================

/**
 * @brief Traits for DiagoDavid solver
 */
struct DiagoDavidTraits
{
    using T = std::complex<double>;
    using SolverType = ::hsolver::DiagoDavid<T, base_device::DEVICE_CPU>;
    using StoragePolicy = RawPointerStorage<T>;

    static constexpr const char* name = "diago_david";
    static constexpr bool uses_f_style = true;  // Column-major arrays
    static constexpr bool has_preconditioner = true;

    // Memory synchronization operation
    using syncmem_op = base_device::memory::synchronize_memory_op<
        T, base_device::DEVICE_CPU, base_device::DEVICE_CPU>;
};

/**
 * @brief Traits for DiagoDavSubspace solver
 */
struct DiagoDavSubspaceTraits
{
    using T = std::complex<double>;
    using SolverType = ::hsolver::Diago_DavSubspace<T, base_device::DEVICE_CPU>;
    using StoragePolicy = RawPointerStorage<T>;

    static constexpr const char* name = "diago_dav_subspace";
    static constexpr bool uses_f_style = true;  // Column-major arrays
    static constexpr bool has_preconditioner = true;

    // Memory synchronization operation
    using syncmem_op = base_device::memory::synchronize_memory_op<
        T, base_device::DEVICE_CPU, base_device::DEVICE_CPU>;
};

#ifdef __ENABLE_ATEN
/**
 * @brief Traits for DiagoCG solver
 */
struct DiagoCGTraits
{
    using T = std::complex<double>;
    using SolverType = ::hsolver::DiagoCG<T, base_device::DEVICE_CPU>;
    using StoragePolicy = TensorStorage<T>;

    static constexpr const char* name = "diago_cg";
    static constexpr bool uses_f_style = false;  // Row-major arrays
    static constexpr bool has_preconditioner = true;

    // Memory synchronization operation for tensor interface
    using syncmem_op = base_device::memory::synchronize_memory_op<
        T, base_device::DEVICE_CPU, base_device::DEVICE_CPU>;
};
#endif // __ENABLE_ATEN

// ============================================================================
// Helper Functions for Creating HPsi/SPsi Lambdas
// ============================================================================

/**
 * @brief Create hpsi_func lambda for raw pointer interface (F-style)
 *
 * Wraps a Python callable to work with ABACUS raw pointer interface.
 * Handles array layout conversion between Python (row-major) and
 * ABACUS (column-major for Davidson methods).
 */
template <typename T>
auto make_hpsi_func_fstyle(
    std::function<py::array_t<T>(py::array_t<T>)> mm_op)
{
    return [mm_op](T* psi_in, T* hpsi_out, const int ld_psi, const int nvec) {
        // Create F-style numpy array (column-major)
        py::array_t<T, py::array::f_style> psi({ld_psi, nvec});
        py::buffer_info buf = psi.request();
        T* ptr = static_cast<T*>(buf.ptr);
        std::copy(psi_in, psi_in + nvec * ld_psi, ptr);

        // Call Python function
        py::array_t<T, py::array::f_style> hpsi = mm_op(psi);

        // Copy result back
        py::buffer_info hpsi_buf = hpsi.request();
        T* hpsi_ptr = static_cast<T*>(hpsi_buf.ptr);
        std::copy(hpsi_ptr, hpsi_ptr + nvec * ld_psi, hpsi_out);
    };
}

/**
 * @brief Create spsi_func lambda for raw pointer interface (identity)
 *
 * For non-orthogonal basis, S*psi = psi (identity operation).
 */
template <typename Traits>
auto make_spsi_func_identity()
{
    using T = typename Traits::T;
    using syncmem_op = typename Traits::syncmem_op;

    return [](const T* psi_in, T* spsi_out, const int nrow, const int nbands) {
        syncmem_op()(spsi_out, psi_in, static_cast<size_t>(nbands * nrow));
    };
}

#ifdef __ENABLE_ATEN
/**
 * @brief Create hpsi_func lambda for tensor interface
 */
template <typename T>
auto make_hpsi_func_tensor(
    std::function<py::array_t<T>(py::array_t<T>)> mm_op)
{
    return [mm_op](const ct::Tensor& psi_in, ct::Tensor& hpsi_out) {
        const auto ndim = psi_in.shape().ndim();
        REQUIRES_OK(ndim <= 2, "dims of psi_in should be less than or equal to 2");
        const int nvec = ndim == 1 ? 1 : psi_in.shape().dim_size(0);
        const int ld_psi = ndim == 1 ? psi_in.NumElements() : psi_in.shape().dim_size(1);

        // Create numpy array (row-major for CG)
        py::array_t<T> psi({ld_psi, nvec});
        py::buffer_info buf = psi.request();
        T* ptr = static_cast<T*>(buf.ptr);
        std::copy(psi_in.data<T>(), psi_in.data<T>() + nvec * ld_psi, ptr);

        // Call Python function
        py::array_t<T> hpsi = mm_op(psi);

        // Copy result back
        py::buffer_info hpsi_buf = hpsi.request();
        T* hpsi_ptr = static_cast<T*>(hpsi_buf.ptr);
        std::copy(hpsi_ptr, hpsi_ptr + nvec * ld_psi, hpsi_out.data<T>());
    };
}

/**
 * @brief Create spsi_func lambda for tensor interface (identity)
 */
template <typename Traits>
auto make_spsi_func_tensor_identity()
{
    using T = typename Traits::T;
    using syncmem_op = typename Traits::syncmem_op;

    return [](const ct::Tensor& psi_in, ct::Tensor& spsi_out) {
        const auto ndim = psi_in.shape().ndim();
        REQUIRES_OK(ndim <= 2, "dims of psi_in should be less than or equal to 2");
        const int nrow = ndim == 1 ? psi_in.NumElements() : psi_in.shape().dim_size(1);
        const int nbands = ndim == 1 ? 1 : psi_in.shape().dim_size(0);
        syncmem_op()(
            spsi_out.data<T>(),
            psi_in.data<T>(),
            static_cast<size_t>(nrow * nbands)
        );
    };
}
#endif // __ENABLE_ATEN

} // namespace hsolver
} // namespace pyabacus

#endif // PYABACUS_HSOLVER_DIAGO_TRAITS_HPP
