/**
 * @file py_hsolver.cpp
 * @brief Python bindings for HSolver diagonalization methods
 *
 * This file provides pybind11 bindings for the diagonalization solvers
 * using the unified adapter template approach.
 */

#include <complex>
#include <functional>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/numpy.h>

#include "source_hsolver/diago_dav_subspace.h"
#include "source_base/kernels/math_kernel_op.h"
#include "source_base/module_device/types.h"

#include "diago_adapter.hpp"

namespace py = pybind11;
using namespace pybind11::literals;

using namespace pyabacus::hsolver;

void bind_hsolver(py::module& m)
{
    // Bind diag_comm_info struct
    py::class_<hsolver::diag_comm_info>(m, "diag_comm_info")
        .def(py::init<const int, const int>(), "rank"_a, "nproc"_a)
        .def_readonly("rank", &hsolver::diag_comm_info::rank)
        .def_readonly("nproc", &hsolver::diag_comm_info::nproc);

    // Bind PyDiagoDavSubspace using adapter
    py::class_<PyDiagoDavSubspaceAdapter>(m, "diago_dav_subspace")
        .def(py::init<int, int>(), R"pbdoc(
            Constructor of diago_dav_subspace, a class for diagonalizing
            a linear operator using the Davidson-Subspace Method.

            This class serves as a backend computation class. The interface
            for invoking this class is a function defined in _hsolver.py,
            which uses this class to perform the calculations.

            Parameters
            ----------
            nbasis : int
                The number of basis functions.
            nband : int
                The number of bands to be calculated.
        )pbdoc", "nbasis"_a, "nband"_a)
        .def("diag", &PyDiagoDavSubspaceAdapter::diag, R"pbdoc(
            Diagonalize the linear operator using the Davidson-Subspace Method.

            Parameters
            ----------
            mm_op : Callable[[NDArray[np.complex128]], NDArray[np.complex128]],
                The operator to be diagonalized, which is a function that takes a matrix as input
                and returns a matrix mv_op(X) = H * X as output.
            precond_vec : np.ndarray
                The preconditioner vector.
            dav_ndim : int
                The number of vectors, which is a multiple of the number of
                eigenvectors to be calculated.
            tol : double
                The tolerance for the convergence.
            max_iter : int
                The maximum number of iterations.
            need_subspace : bool
                Whether to use the subspace function.
            diag_ethr : List[float] | None, optional
                The list of thresholds of bands, by default None.
            scf_type : bool
                Whether to use the SCF type, which is used to determine the
                convergence criterion.
                If true, it indicates a self-consistent field (SCF) calculation,
                where the initial precision of eigenvalue calculation can be coarse.
                If false, it indicates a non-self-consistent field (non-SCF) calculation,
                where high precision in eigenvalue calculation is required from the start.
            comm_info : diag_comm_info
                The communicator information.
            diago_subspace : int
                The method to solve the generalized eigenvalue problem.
                0: LAPACK, 1: Gen-ELPA, 2: ScaLAPACK
            nb2d : int
                The block size in 2d block cyclic distribution if use elpa or scalapack.
        )pbdoc",
        "mm_op"_a,
        "precond_vec"_a,
        "dav_ndim"_a,
        "tol"_a,
        "max_iter"_a,
        "need_subspace"_a,
        "diag_ethr"_a,
        "scf_type"_a,
        "comm_info"_a,
        "diago_subspace"_a,
        "nb2d"_a)
        .def("set_psi", &PyDiagoDavSubspaceAdapter::set_psi, R"pbdoc(
            Set the initial guess of the eigenvectors, i.e. the wave functions.
        )pbdoc", "psi_in"_a)
        .def("get_psi", &PyDiagoDavSubspaceAdapter::get_psi, R"pbdoc(
            Get the eigenvectors.
        )pbdoc")
        .def("init_eigenvalue", &PyDiagoDavSubspaceAdapter::init_eigenvalue, R"pbdoc(
            Initialize the eigenvalues as zero.
        )pbdoc")
        .def("get_eigenvalue", &PyDiagoDavSubspaceAdapter::get_eigenvalue, R"pbdoc(
            Get the eigenvalues.
        )pbdoc");

    // Bind PyDiagoDavid using adapter
    py::class_<PyDiagoDavidAdapter>(m, "diago_david")
        .def(py::init<int, int>(), R"pbdoc(
            Constructor of diago_david, a class for diagonalizing
            a linear operator using the Davidson Method.

            This class serves as a backend computation class. The interface
            for invoking this class is a function defined in _hsolver.py,
            which uses this class to perform the calculations.

            Parameters
            ----------
            nbasis : int
                The number of basis functions.
            nband : int
                The number of bands to be calculated.
        )pbdoc", "nbasis"_a, "nband"_a)
        .def("diag", &PyDiagoDavidAdapter::diag, R"pbdoc(
            Diagonalize the linear operator using the Davidson Method.

            Parameters
            ----------
            mm_op : Callable[[NDArray[np.complex128]], NDArray[np.complex128]],
                The operator to be diagonalized, which is a function that takes a matrix as input
                and returns a matrix mv_op(X) = H * X as output.
            precond_vec : np.ndarray
                The preconditioner vector.
            dav_ndim : int
                The number of vectors, which is a multiple of the number of
                eigenvectors to be calculated.
            tol : double
                The tolerance for the convergence.
            diag_ethr: np.ndarray
                The tolerance vector.
            max_iter : int
                The maximum number of iterations.
        )pbdoc",
        "mm_op"_a,
        "precond_vec"_a,
        "dav_ndim"_a,
        "tol"_a,
        "diag_ethr"_a,
        "max_iter"_a,
        "comm_info"_a)
        .def("set_psi", &PyDiagoDavidAdapter::set_psi, R"pbdoc(
            Set the initial guess of the eigenvectors, i.e. the wave functions.
        )pbdoc", "psi_in"_a)
        .def("get_psi", &PyDiagoDavidAdapter::get_psi, R"pbdoc(
            Get the eigenvectors.
        )pbdoc")
        .def("init_eigenvalue", &PyDiagoDavidAdapter::init_eigenvalue, R"pbdoc(
            Initialize the eigenvalues as zero.
        )pbdoc")
        .def("get_eigenvalue", &PyDiagoDavidAdapter::get_eigenvalue, R"pbdoc(
            Get the eigenvalues.
        )pbdoc");

#ifdef __ENABLE_ATEN
    // Bind PyDiagoCG using adapter (only when ATen is available)
    py::class_<PyDiagoCGAdapter>(m, "diago_cg")
        .def(py::init<int, int>(), R"pbdoc(
            Constructor of diago_cg, a class for diagonalizing
            a linear operator using the Conjugate Gradient Method.

            This class serves as a backend computation class. The interface
            for invoking this class is a function defined in _hsolver.py,
            which uses this class to perform the calculations.
        )pbdoc")
        .def("diag",
             &PyDiagoCGAdapter::diag,
             R"pbdoc(
            Diagonalize the linear operator using the Conjugate Gradient Method.

            Parameters
            ----------
            mm_op : Callable[[NDArray[np.complex128]], NDArray[np.complex128]],
                The operator to be diagonalized, which is a function that takes a matrix as input
                and returns a matrix mv_op(X) = H * X as output.
            max_iter : int
                The maximum number of iterations.
            tol : double
                The tolerance for the convergence.
            need_subspace : bool
                Whether to use the subspace function.
            scf_type : bool
                Whether to use the SCF type, which is used to determine the
                convergence criterion.
        )pbdoc",
        "mm_op"_a,
        "max_iter"_a,
        "tol"_a,
        "diag_ethr"_a,
        "need_subspace"_a,
        "scf_type"_a,
        "nproc_in_pool"_a)
        .def("init_eig", &PyDiagoCGAdapter::init_eig, R"pbdoc(
            Initialize the eigenvalues.
        )pbdoc")
        .def("get_eig", &PyDiagoCGAdapter::get_eig, R"pbdoc(
            Get the eigenvalues.
        )pbdoc")
        .def("set_psi", &PyDiagoCGAdapter::set_psi, R"pbdoc(
            Set the eigenvectors.
        )pbdoc", "psi_in"_a)
        .def("get_psi", &PyDiagoCGAdapter::get_psi, R"pbdoc(
            Get the eigenvectors.
        )pbdoc")
        .def("set_prec", &PyDiagoCGAdapter::set_prec, R"pbdoc(
            Set the preconditioner.
        )pbdoc", "prec_in"_a);
#else
    // Provide stub binding when ATen is not available
    // This allows the module to load but will raise an error if used
    m.def("diago_cg_available", []() { return false; },
          "Check if diago_cg is available (requires ATen)");
#endif
}

PYBIND11_MODULE(_hsolver_pack, m)
{
    m.doc() = "Submodule for pyabacus: hsolver";

    bind_hsolver(m);
}
