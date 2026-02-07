#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "source_base/math_sphbes.h"
#include "source_base/math_integral.h"
#include "source_base/spherical_bessel_transformer.h"

#include "../utils/pybind_utils.h"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace pyabacus::utils;

template <typename... Args>
using overload_cast_ = pybind11::detail::overload_cast_impl<Args...>;

void bind_base_math(py::module& m)
{
    // python binding for class Sphbes
    py::class_<ModuleBase::Sphbes>(m, "Sphbes")
        .def(py::init<>())
        .def_static("sphbesj", overload_cast_<const int, const double>()(&ModuleBase::Sphbes::sphbesj), "l"_a, "x"_a)
        .def_static("dsphbesj", overload_cast_<const int, const double>()(&ModuleBase::Sphbes::dsphbesj), "l"_a, "x"_a)
        .def_static("sphbesj",
                    [](const int n, py::array_t<double> r, const double q, const int l, py::array_t<double> jl) {
                        check_1d_array(r, "r");
                        check_1d_array(jl, "jl");
                        ModuleBase::Sphbes::sphbesj(n,
                                                    get_array_ptr(r),
                                                    q,
                                                    l,
                                                    get_array_ptr(jl));
                    })
        .def_static("dsphbesj",
                    [](const int n, py::array_t<double> r, const double q, const int l, py::array_t<double> djl) {
                        check_1d_array(r, "r");
                        check_1d_array(djl, "djl");
                        ModuleBase::Sphbes::dsphbesj(n,
                                                     get_array_ptr(r),
                                                     q,
                                                     l,
                                                     get_array_ptr(djl));
                    })
        .def_static("sphbes_zeros", [](const int l, const int n, py::array_t<double> zeros) {
            check_1d_array(zeros, "zeros");
            ModuleBase::Sphbes::sphbes_zeros(l, n, get_array_ptr(zeros));
        });

    // python binding for class Integral
    py::class_<ModuleBase::Integral>(m, "Integral")
        .def(py::init<>())
        .def_static("Simpson_Integral", [](const int mesh, py::array_t<double> func, py::array_t<double> rab, double asum) {
            check_1d_array(func, "func");
            check_1d_array(rab, "rab");

            double isum = asum;
            ModuleBase::Integral::Simpson_Integral(mesh,
                                                    get_array_ptr(func),
                                                    get_array_ptr(rab),
                                                    isum);
            return isum;
        })
        .def_static("Simpson_Integral", [](const int mesh, py::array_t<double> func, const double dr, double asum){
            check_1d_array(func, "func");

            double isum = asum;
            ModuleBase::Integral::Simpson_Integral(mesh,
                                                    get_array_ptr(func),
                                                    dr,
                                                    isum);
            return isum;
        })
        .def_static("Simpson_Integral_0toall", [](const int mesh, py::array_t<double> func, py::array_t<double> rab, py::array_t<double> asum){
            check_1d_array(func, "func");
            check_1d_array(rab, "rab");
            check_1d_array(asum, "asum");
            ModuleBase::Integral::Simpson_Integral_0toall(mesh,
                                                            get_array_ptr(func),
                                                            get_array_ptr(rab),
                                                            get_array_ptr(asum));
        })
        .def_static("Simpson_Integral_alltoinf", [](const int mesh, py::array_t<double> func, py::array_t<double> rab, py::array_t<double> asum){
            check_1d_array(func, "func");
            check_1d_array(rab, "rab");
            check_1d_array(asum, "asum");
            ModuleBase::Integral::Simpson_Integral_alltoinf(mesh,
                                                            get_array_ptr(func),
                                                            get_array_ptr(rab),
                                                            get_array_ptr(asum));
        })
        .def_static("simpson", [](const int n, py::array_t<double> f, const double dx){
            check_1d_array(f, "f");
            return ModuleBase::Integral::simpson(n,
                                                    get_array_ptr(f),
                                                    dx);
        })
        .def_static("simpson", [](const int n, py::array_t<double> f, py::array_t<double> h){
            check_1d_array(f, "f");
            check_1d_array(h, "h");
            return ModuleBase::Integral::simpson(n,
                                                    get_array_ptr(f),
                                                    get_array_ptr(h));
        })
        .def_static("Gauss_Legendre_grid_and_weight", [](const int n, py::array_t<double> x, py::array_t<double> w){
            check_1d_array(x, "x");
            check_1d_array(w, "w");
            ModuleBase::Integral::Gauss_Legendre_grid_and_weight(n,
                                                                    get_array_ptr(x),
                                                                    get_array_ptr(w));
        })
        .def_static("Gauss_Legendre_grid_and_weight", [](const double xmin, const double xmax, const int n, py::array_t<double> x, py::array_t<double> w){
            check_1d_array(x, "x");
            check_1d_array(w, "w");
            ModuleBase::Integral::Gauss_Legendre_grid_and_weight(xmin,
                                                                    xmax,
                                                                    n,
                                                                    get_array_ptr(x),
                                                                    get_array_ptr(w));
        });
    py::class_<ModuleBase::SphericalBesselTransformer>(m, "SphericalBesselTransformer")
        .def(py::init<>());
}

PYBIND11_MODULE(_base_pack, m)
{
    m.doc() = "Submodule for pyabacus: ModuleBase";

    bind_base_math(m);
}
