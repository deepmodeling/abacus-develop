#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "module_base/math_sphbes.h"
#include "module_base/math_integral.h"

namespace py = pybind11;
using namespace pybind11::literals;
template <typename... Args>
using overload_cast_ = pybind11::detail::overload_cast_impl<Args...>;

void bind_math_base(py::module& m)
{
    py::module module_base = m.def_submodule("ModuleBase");

    // python binding for class Sphbes
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

    // python binding for class Integral
    py::class_<ModuleBase::Integral>(module_base, "Integral")
        .def(py::init<>())
        .def_static("Simpson_Integral", [](const int mesh, py::array_t<double> func, py::array_t<double> rab, double asum) {
            py::buffer_info func_info = func.request();
            if (func_info.ndim != 1)
            {
                throw std::runtime_error("func array must be 1-dimensional");
            }
            py::buffer_info rab_info = rab.request();
            if (rab.ndim() != 1)
            {
                throw std::runtime_error("rab array must be 1-dimensional");
            }

            double isum = asum;
            ModuleBase::Integral::Simpson_Integral(mesh,
                                                    static_cast<const double* const>(func_info.ptr),
                                                    static_cast<const double* const>(rab_info.ptr),
                                                    isum);
            return isum;
        })
        .def_static("Simpson_Integral", [](const int mesh, py::array_t<double> func, const double dr, double asum){
            py::buffer_info func_info = func.request();
            if (func_info.ndim != 1)
            {
                throw std::runtime_error("func array must be 1-dimensional");
            }
                
                double isum = asum;
                ModuleBase::Integral::Simpson_Integral(mesh,
                                                        static_cast<const double* const>(func_info.ptr),
                                                        dr,
                                                        isum);
                return isum;
        })
        .def_static("Simpson_Integral_0toall", [](const int mesh, py::array_t<double> func, py::array_t<double> rab, py::array_t<double> asum){
            py::buffer_info func_info = func.request();
            if (func_info.ndim != 1)
            {
                throw std::runtime_error("func array must be 1-dimensional");
            }
            py::buffer_info rab_info = rab.request();
            if (rab.ndim() != 1)
            {
                throw std::runtime_error("rab array must be 1-dimensional");
            }
            py::buffer_info asum_info = asum.request();
            if (asum.ndim() != 1)
            {
                throw std::runtime_error("asum array must be 1-dimensional");
            }
            ModuleBase::Integral::Simpson_Integral_0toall(mesh,
                                                            static_cast<const double* const>(func_info.ptr),
                                                            static_cast<const double* const>(rab_info.ptr),
                                                            static_cast<double* const>(asum_info.ptr));
        })
        .def_static("Simpson_Integral_alltoinf", [](const int mesh, py::array_t<double> func, py::array_t<double> rab, py::array_t<double> asum){
            py::buffer_info func_info = func.request();
            if (func_info.ndim != 1)
            {
                throw std::runtime_error("func array must be 1-dimensional");
            }
            py::buffer_info rab_info = rab.request();
            if (rab.ndim() != 1)
            {
                throw std::runtime_error("rab array must be 1-dimensional");
            }
            py::buffer_info asum_info = asum.request();
            if (asum.ndim() != 1)
            {
                throw std::runtime_error("asum array must be 1-dimensional");
            }
            ModuleBase::Integral::Simpson_Integral_alltoinf(mesh,
                                                            static_cast<const double* const>(func_info.ptr),
                                                            static_cast<const double* const>(rab_info.ptr),
                                                            static_cast<double* const>(asum_info.ptr));
        })
        .def_static("Simpson_Integral_alltoinf", [](const int mesh, py::array_t<double> func, py::array_t<double> rab, py::array_t<double> asum){
            py::buffer_info func_info = func.request();
            if (func_info.ndim != 1)
            {
                throw std::runtime_error("func array must be 1-dimensional");
            }
            py::buffer_info rab_info = rab.request();
            if (rab.ndim() != 1)
            {
                throw std::runtime_error("rab array must be 1-dimensional");
            }
            py::buffer_info asum_info = asum.request();
            if (asum.ndim() != 1)
            {
                throw std::runtime_error("asum array must be 1-dimensional");
            }
            ModuleBase::Integral::Simpson_Integral_alltoinf(mesh,
                                                            static_cast<const double* const>(func_info.ptr),
                                                            static_cast<const double* const>(rab_info.ptr),
                                                            static_cast<double* const>(asum_info.ptr));
        })
        .def_static("simpson", [](const int n, py::array_t<double> f, const double dx){
            py::buffer_info f_info = f.request();
            if (f_info.ndim != 1)
            {
                throw std::runtime_error("f array must be 1-dimensional");
            }
            return ModuleBase::Integral::simpson(n,
                                                    static_cast<const double* const>(f_info.ptr),
                                                    dx);
        })
        .def_static("simpson", [](const int n, py::array_t<double> f, py::array_t<double> h){
            py::buffer_info f_info = f.request();
            if (f_info.ndim != 1)
            {
                throw std::runtime_error("f array must be 1-dimensional");
            }
            py::buffer_info h_info = h.request();
            if (h.ndim() != 1)
            {
                throw std::runtime_error("h array must be 1-dimensional");
            }
            return ModuleBase::Integral::simpson(n,
                                                    static_cast<const double* const>(f_info.ptr),
                                                    static_cast<const double* const>(h_info.ptr));
        })
        .def_static("Gauss_Legendre_grid_and_weight", [](const int n, py::array_t<double> x, py::array_t<double> w){
            py::buffer_info x_info = x.request();
            if (x_info.ndim != 1)
            {
                throw std::runtime_error("x array must be 1-dimensional");
            }
            py::buffer_info w_info = w.request();
            if (w.ndim() != 1)
            {
                throw std::runtime_error("w array must be 1-dimensional");
            }
            ModuleBase::Integral::Gauss_Legendre_grid_and_weight(n,
                                                                    static_cast<double*>(x_info.ptr),
                                                                    static_cast<double*>(w_info.ptr));
        })
        .def_static("Gauss_Legendre_grid_and_weight", [](const double xmin, const double xmax, const int n, py::array_t<double> x, py::array_t<double> w){
            py::buffer_info x_info = x.request();
            if (x_info.ndim != 1)
            {
                throw std::runtime_error("x array must be 1-dimensional");
            }
            py::buffer_info w_info = w.request();
            if (w.ndim() != 1)
            {
                throw std::runtime_error("w array must be 1-dimensional");
            }
            ModuleBase::Integral::Gauss_Legendre_grid_and_weight(xmin,
                                                                    xmax,
                                                                    n,
                                                                    static_cast<double*>(x_info.ptr),
                                                                    static_cast<double*>(w_info.ptr));
        });
}