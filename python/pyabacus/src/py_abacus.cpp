#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

void bind_base_math(py::module& m);
void bind_m_nao(py::module& m);

PYBIND11_MODULE(_core, m)
{
    bind_base_math(m);
    bind_m_nao(m);
}