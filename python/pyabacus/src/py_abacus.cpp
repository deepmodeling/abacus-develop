#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

void bind_base_math(py::module& m);
void bind_m_nao(py::module& m);

PYBIND11_MODULE(_core, m)
{
    m.doc() = "Python extension for ABACUS built with pybind11 and scikit-build.";
    bind_base_math(m);
    bind_m_nao(m);
}