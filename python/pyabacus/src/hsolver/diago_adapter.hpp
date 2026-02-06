/**
 * @file diago_adapter.hpp
 * @brief Template adapter for diagonalizer wrappers
 *
 * This header provides a unified template adapter that wraps different
 * diagonalization solvers (David, DavSubspace, CG) with a common interface.
 */

#ifndef PYABACUS_HSOLVER_DIAGO_ADAPTER_HPP
#define PYABACUS_HSOLVER_DIAGO_ADAPTER_HPP

#include <complex>
#include <memory>
#include <functional>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>

#include "diago_traits.hpp"
#include "../utils/pybind_utils.h"

namespace py = pybind11;

namespace pyabacus {
namespace hsolver {

// ============================================================================
// PyDiagoDavid Adapter
// ============================================================================

/**
 * @brief Adapter for DiagoDavid solver
 */
class PyDiagoDavidAdapter
{
public:
    using Traits = DiagoDavidTraits;
    using T = typename Traits::T;
    using SolverType = typename Traits::SolverType;

    PyDiagoDavidAdapter(int nbasis, int nband)
        : nbasis_(nbasis), nband_(nband)
    {
        storage_.allocate(nbasis, nband);
    }

    PyDiagoDavidAdapter(const PyDiagoDavidAdapter&) = delete;
    PyDiagoDavidAdapter& operator=(const PyDiagoDavidAdapter&) = delete;

    PyDiagoDavidAdapter(PyDiagoDavidAdapter&& other) noexcept
        : storage_(std::move(other.storage_))
        , nbasis_(other.nbasis_)
        , nband_(other.nband_)
    {
    }

    void set_psi(py::array_t<T> psi_in)
    {
        storage_.set_psi(psi_in);
    }

    py::array_t<T> get_psi() const
    {
        return storage_.get_psi();
    }

    void init_eigenvalue()
    {
        storage_.init_eigenvalue();
    }

    py::array_t<double> get_eigenvalue() const
    {
        return storage_.get_eigenvalue();
    }

    int diag(
        std::function<py::array_t<T>(py::array_t<T>)> mm_op,
        std::vector<double>& precond_vec,
        int dav_ndim,
        double tol,
        std::vector<double>& diag_ethr,
        int max_iter,
        bool use_paw,
        ::hsolver::diag_comm_info comm_info)
    {
        auto hpsi_func = make_hpsi_func_fstyle<T>(mm_op);
        auto spsi_func = make_spsi_func_identity<Traits>();

        solver_ = std::make_unique<SolverType>(
            precond_vec.data(),
            nband_,
            nbasis_,
            dav_ndim,
            use_paw,
            comm_info
        );

        return solver_->diag(
            hpsi_func,
            spsi_func,
            nbasis_,
            storage_.psi_ptr(),
            storage_.eigenvalue_ptr(),
            diag_ethr,
            max_iter
        );
    }

private:
    RawPointerStorage<T> storage_;
    std::unique_ptr<SolverType> solver_;
    int nbasis_;
    int nband_;
};

// ============================================================================
// PyDiagoDavSubspace Adapter
// ============================================================================

/**
 * @brief Adapter for DiagoDavSubspace solver
 */
class PyDiagoDavSubspaceAdapter
{
public:
    using Traits = DiagoDavSubspaceTraits;
    using T = typename Traits::T;
    using SolverType = typename Traits::SolverType;

    PyDiagoDavSubspaceAdapter(int nbasis, int nband)
        : nbasis_(nbasis), nband_(nband)
    {
        storage_.allocate(nbasis, nband);
    }

    PyDiagoDavSubspaceAdapter(const PyDiagoDavSubspaceAdapter&) = delete;
    PyDiagoDavSubspaceAdapter& operator=(const PyDiagoDavSubspaceAdapter&) = delete;

    PyDiagoDavSubspaceAdapter(PyDiagoDavSubspaceAdapter&& other) noexcept
        : storage_(std::move(other.storage_))
        , nbasis_(other.nbasis_)
        , nband_(other.nband_)
    {
    }

    void set_psi(py::array_t<T> psi_in)
    {
        storage_.set_psi(psi_in);
    }

    py::array_t<T> get_psi() const
    {
        return storage_.get_psi();
    }

    void init_eigenvalue()
    {
        storage_.init_eigenvalue();
    }

    py::array_t<double> get_eigenvalue() const
    {
        return storage_.get_eigenvalue();
    }

    int diag(
        std::function<py::array_t<T>(py::array_t<T>)> mm_op,
        std::vector<double>& precond_vec,
        int dav_ndim,
        double tol,
        int max_iter,
        bool need_subspace,
        std::vector<double>& diag_ethr,
        bool scf_type,
        ::hsolver::diag_comm_info comm_info,
        int diag_subspace,
        int nb2d)
    {
        auto hpsi_func = make_hpsi_func_fstyle<T>(mm_op);
        auto spsi_func = make_spsi_func_identity<Traits>();

        solver_ = std::make_unique<SolverType>(
            precond_vec,
            nband_,
            nbasis_,
            dav_ndim,
            tol,
            max_iter,
            comm_info,
            diag_subspace,
            nb2d
        );

        return solver_->diag(
            hpsi_func,
            spsi_func,
            storage_.psi_ptr(),
            nbasis_,
            storage_.eigenvalue_ptr(),
            diag_ethr,
            scf_type
        );
    }

private:
    RawPointerStorage<T> storage_;
    std::unique_ptr<SolverType> solver_;
    int nbasis_;
    int nband_;
};

#ifdef __ENABLE_ATEN
// ============================================================================
// PyDiagoCG Adapter
// ============================================================================

/**
 * @brief Adapter for DiagoCG solver
 */
class PyDiagoCGAdapter
{
public:
    using Traits = DiagoCGTraits;
    using T = typename Traits::T;
    using SolverType = typename Traits::SolverType;

    PyDiagoCGAdapter(int dim, int num_eigs)
        : dim_(dim), num_eigs_(num_eigs)
    {
        storage_.allocate(dim, num_eigs);
    }

    PyDiagoCGAdapter(const PyDiagoCGAdapter&) = delete;
    PyDiagoCGAdapter& operator=(const PyDiagoCGAdapter&) = delete;

    PyDiagoCGAdapter(PyDiagoCGAdapter&& other) noexcept
        : storage_(std::move(other.storage_))
        , dim_(other.dim_)
        , num_eigs_(other.num_eigs_)
    {
    }

    void set_psi(py::array_t<T> psi_in)
    {
        storage_.set_psi(psi_in);
    }

    py::array_t<T> get_psi() const
    {
        return storage_.get_psi();
    }

    void init_eig()
    {
        storage_.init_eigenvalue();
    }

    py::array_t<double> get_eig() const
    {
        return storage_.get_eigenvalue();
    }

    void set_prec(py::array_t<double> prec_in)
    {
        storage_.set_preconditioner(prec_in);
    }

    void diag(
        std::function<py::array_t<T>(py::array_t<T>)> mm_op,
        int diag_ndim,
        double tol,
        const std::vector<double>& diag_ethr,
        bool need_subspace,
        bool scf_type,
        int nproc_in_pool = 1)
    {
        const std::string basis_type = "pw";
        const std::string calculation = scf_type ? "scf" : "nscf";

        auto hpsi_func = make_hpsi_func_tensor<T>(mm_op);
        auto spsi_func = make_spsi_func_tensor_identity<Traits>();
        auto subspace_func = [](const ct::Tensor& psi_in, ct::Tensor& psi_out, const bool S_orth) {
            // Do nothing - placeholder
        };

        solver_ = std::make_unique<SolverType>(
            basis_type,
            calculation,
            need_subspace,
            subspace_func,
            tol,
            diag_ndim,
            nproc_in_pool
        );

        solver_->diag(
            hpsi_func,
            spsi_func,
            *storage_.psi_tensor(),
            *storage_.eig_tensor(),
            diag_ethr,
            *storage_.prec_tensor()
        );
    }

private:
    TensorStorage<T> storage_;
    std::unique_ptr<SolverType> solver_;
    int dim_;
    int num_eigs_;
};
#endif // __ENABLE_ATEN

// ============================================================================
// Backward Compatibility Aliases
// ============================================================================

namespace py_hsolver_compat {

using PyDiagoDavid = PyDiagoDavidAdapter;
using PyDiagoDavSubspace = PyDiagoDavSubspaceAdapter;

#ifdef __ENABLE_ATEN
using PyDiagoCG = PyDiagoCGAdapter;
#endif

} // namespace py_hsolver_compat

} // namespace hsolver
} // namespace pyabacus

#endif // PYABACUS_HSOLVER_DIAGO_ADAPTER_HPP
