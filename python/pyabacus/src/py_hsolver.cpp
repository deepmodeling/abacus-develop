#include <complex>
#include <functional>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/numpy.h>

#include "module_hsolver/diago_dav_subspace.h"
#include "module_hsolver/kernels/math_kernel_op.h"
#include "module_base/module_device/types.h"

#include "./py_diago_dav_subspace.hpp"
#include "./py_diago_david.hpp"

namespace py = pybind11;
using namespace pybind11::literals;

void bind_hsolver(py::module& m)
{
    py::module hsolver = m.def_submodule("hsolver");

    py::class_<hsolver::diag_comm_info>(hsolver, "diag_comm_info")
        .def(py::init<const int, const int>(), "rank"_a, "nproc"_a)
        .def_readonly("rank", &hsolver::diag_comm_info::rank)
        .def_readonly("nproc", &hsolver::diag_comm_info::nproc);

    py::class_<py_hsolver::PyDiagoDavSubspace>(hsolver, "diago_dav_subspace")
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
        .def("diag", &py_hsolver::PyDiagoDavSubspace::diag, R"pbdoc(
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
            diag_ethr : list[float]
                A list of float values indicating the thresholds of each band for the diagonalization,
                meaning that the corresponding eigenvalue is to be calculated.
            scf_type : bool
                Whether to use the SCF type, which is used to determine the
                convergence criterion.
                If true, it indicates a self-consistent field (SCF) calculation,
                where the initial precision of eigenvalue calculation can be coarse.
                If false, it indicates a non-self-consistent field (non-SCF) calculation,
                where high precision in eigenvalue calculation is required from the start.
        )pbdoc", 
        "mm_op"_a, 
        "precond_vec"_a, 
        "dav_ndim"_a, 
        "tol"_a, 
        "max_iter"_a, 
        "need_subspace"_a, 
        "diag_ethr"_a, 
        "scf_type"_a, 
        "comm_info"_a)
        .def("set_psi", &py_hsolver::PyDiagoDavSubspace::set_psi, R"pbdoc(
            Set the initial guess of the eigenvectors, i.e. the wave functions.
        )pbdoc", "psi_in"_a)
        .def("get_psi", &py_hsolver::PyDiagoDavSubspace::get_psi, R"pbdoc(
            Get the eigenvectors.
        )pbdoc")
        .def("init_eigenvalue", &py_hsolver::PyDiagoDavSubspace::init_eigenvalue, R"pbdoc(
            Initialize the eigenvalues as zero.
        )pbdoc")
        .def("get_eigenvalue", &py_hsolver::PyDiagoDavSubspace::get_eigenvalue, R"pbdoc(
            Get the eigenvalues.        
        )pbdoc");

    py::class_<py_hsolver::PyDiagoDavid>(hsolver, "diago_david")
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
        .def("diag", &py_hsolver::PyDiagoDavid::diag, R"pbdoc(
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
            max_iter : int
                The maximum number of iterations.
            use_paw : bool
                Whether to use the projector augmented wave method.
        )pbdoc", 
        "mm_op"_a, 
        "precond_vec"_a, 
        "dav_ndim"_a, 
        "tol"_a, 
        "max_iter"_a, 
        "use_paw"_a, 
        "comm_info"_a)
        .def("set_psi", &py_hsolver::PyDiagoDavid::set_psi, R"pbdoc(
            Set the initial guess of the eigenvectors, i.e. the wave functions.
        )pbdoc", "psi_in"_a)
        .def("get_psi", &py_hsolver::PyDiagoDavid::get_psi, R"pbdoc(
            Get the eigenvectors.
        )pbdoc")
        .def("init_eigenvalue", &py_hsolver::PyDiagoDavid::init_eigenvalue, R"pbdoc(
            Initialize the eigenvalues as zero.
        )pbdoc")
        .def("get_eigenvalue", &py_hsolver::PyDiagoDavid::get_eigenvalue, R"pbdoc(
            Get the eigenvalues.        
        )pbdoc");
}
