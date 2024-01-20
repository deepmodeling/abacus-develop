#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "module_base/math_sphbes.h"

namespace py = pybind11;
using namespace pybind11::literals;
template <typename... Args>
using overload_cast_ = pybind11::detail::overload_cast_impl<Args...>;

void bind_math_base(py::module& m)
{
    py::module module_base = m.def_submodule("ModuleBase");

    py::class_<ModuleBase::Sphbes>(module_base, "Sphbes")
        .def(py::init<>())
        .def_static("sphbesj", overload_cast_<const int, const double>()(&ModuleBase::Sphbes::sphbesj), "l"_a, "x"_a)
        .def_static("dsphbesj", overload_cast_<const int, const double>()(&ModuleBase::Sphbes::dsphbesj), "l"_a, "x"_a)
        .def_static("sphbesj",
                    [](const int n, py::array_t<double> r, const double q, const int l, py::array_t<double> jl) {
                        py::buffer_info r_info = r.request();
                        if (r_info.ndim != 1)
                        {
                            throw std::runtime_error("r array must be 1-dimensional");
                        }
                        py::buffer_info jl_info = jl.request();
                        if (jl_info.ndim != 1)
                        {
                            throw std::runtime_error("jl array must be 1-dimensional");
                        }
                        ModuleBase::Sphbes::sphbesj(n,
                                                    static_cast<const double* const>(r_info.ptr),
                                                    q,
                                                    l,
                                                    static_cast<double* const>(jl_info.ptr));
                    })
        .def_static("dsphbesj",
                    [](const int n, py::array_t<double> r, const double q, const int l, py::array_t<double> djl) {
                        py::buffer_info r_info = r.request();
                        if (r_info.ndim != 1)
                        {
                            throw std::runtime_error("r array must be 1-dimensional");
                        }
                        py::buffer_info djl_info = djl.request();
                        if (djl_info.ndim != 1)
                        {
                            throw std::runtime_error("djl array must be 1-dimensional");
                        }
                        ModuleBase::Sphbes::dsphbesj(n,
                                                     static_cast<const double* const>(r_info.ptr),
                                                     q,
                                                     l,
                                                     static_cast<double* const>(djl_info.ptr));
                    })
        .def_static("sphbes_zeros", [](const int l, const int n, py::array_t<double> zeros) {
            py::buffer_info zeros_info = zeros.request();
            if (zeros_info.ndim != 1)
            {
                throw std::runtime_error("zeros array must be 1-dimensional");
            }
            ModuleBase::Sphbes::sphbes_zeros(l, n, static_cast<double* const>(zeros_info.ptr));
        });
}