#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

void bind_numerical_radial(py::module& m);
void bind_math_base(py::module& m);

PYBIND11_MODULE(_core, m)
{
    // bind_numerical_radial(m);
    bind_math_base(m);
}