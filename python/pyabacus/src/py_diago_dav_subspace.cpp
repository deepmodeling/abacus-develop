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

namespace py = pybind11;
using namespace pybind11::literals;

void bind_diago_dav_subspace(py::module& m)
{
    py::module hsolver = m.def_submodule("hsolver");

    py::class_<hsolver::diag_comm_info>(hsolver, "diag_comm_info")
        .def(py::init<const int, const int>(), "rank"_a, "nproc"_a)
        .def_readonly("rank", &hsolver::diag_comm_info::rank)
        .def_readonly("nproc", &hsolver::diag_comm_info::nproc);

    py::class_<py_hsolver::PyDiagoDavSubspace>(hsolver, "diago_dav_subspace")
        .def(py::init<int, int>())
        .def("diag", &py_hsolver::PyDiagoDavSubspace::diag)
        .def("set_psi", &py_hsolver::PyDiagoDavSubspace::set_psi)
        .def("get_psi", &py_hsolver::PyDiagoDavSubspace::get_psi)
        .def("init_eigenvalue", &py_hsolver::PyDiagoDavSubspace::init_eigenvalue)
        .def("get_eigenvalue", &py_hsolver::PyDiagoDavSubspace::get_eigenvalue);
}
