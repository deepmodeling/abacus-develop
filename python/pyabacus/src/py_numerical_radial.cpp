#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "source_basis/module_nao/numerical_radial.h"

#include "utils/pybind_utils.h"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace pyabacus::utils;

template <typename... Args>
using overload_cast_ = pybind11::detail::overload_cast_impl<Args...>;

void bind_numerical_radial(py::module& m)
{
    // Create the submodule for NumericalRadial
    py::module m_numerical_radial = m.def_submodule("NumericalRadial");

    py::class_<NumericalRadial>(m_numerical_radial, "NumericalRadial")
        .def(py::init<>())
        .def(
            "build",
            [](NumericalRadial& self,
               const int l,
               const bool for_r_space,
               const int ngrid,
               py::array_t<double> grid,
               py::array_t<double> value,
               const int p = 0,
               const int izeta = 0,
               const std::string symbol = "",
               const int itype,
               const bool init_sbt = true) {
                check_1d_array(grid, "grid");
                check_1d_array(value, "value");
                check_same_size(grid, value, "grid", "value");

                self.build(l,
                           for_r_space,
                           static_cast<int>(grid.size()),
                           get_array_ptr(grid),
                           get_array_ptr(value),
                           p,
                           izeta,
                           symbol,
                           itype,
                           init_sbt);
            },
            "l"_a,
            "for_r_space"_a,
            "ngrid"_a,
            "grid"_a,
            "value"_a,
            "p"_a = 0,
            "izeta"_a = 0,
            "symbol"_a = "",
            "itype"_a = 0,
            "init_sbt"_a = true)
        // leave set_transformer for future since no wrapper for Transformer yet
        .def(
            "set_grid",
            [](NumericalRadial& self,
               const bool for_r_space,
               const int ngrid,
               py::array_t<double> grid,
               const char mode) {
                check_1d_array(grid, "grid");

                self.set_grid(for_r_space, ngrid, get_array_ptr(grid), mode);
            },
            "for_r_space"_a,
            "ngrid"_a,
            "grid"_a,
            "mode"_a = 'i')
        .def("set_uniform_grid",
             &NumericalRadial::set_uniform_grid,
             "for_r_space"_a,
             "ngrid"_a,
             "cutoff"_a,
             "mode"_a = 'i',
             "enable_fft"_a = false)
        .def(
            "set_value",
            [](NumericalRadial& self, const bool for_r_space, py::array_t<double> value, const int p) {
                check_1d_array(value, "value");

                self.set_value(for_r_space, get_array_ptr(value), p);
            },
            "for_r_space"_a,
            "value"_a,
            "p"_a)
        .def("wipe", &NumericalRadial::wipe, "r_space"_a = true, "k_space"_a = true)
        .def(
            "radtab",
            [](NumericalRadial& self,
               const char op,
               NumericalRadial& ket,
               const int l,
               py::array_t<double> table,
               const int nr_tab,
               const double rmax_tab,
               const bool deriv) {
                check_1d_array(table, "table");

                self.radtab(op, ket, l, get_array_ptr(table), nr_tab, rmax_tab, deriv);
            },
            "op"_a,
            "ket"_a,
            "l"_a,
            "table"_a,
            "nr_tab"_a,
            "rmax_tab"_a,
            "deriv"_a = false)
        .def("normalize", &NumericalRadial::normalize, "for_r_space"_a = true)
        // Getters
        .def_property_readonly("symbol", &NumericalRadial::symbol)
        .def_property_readonly("itype", &NumericalRadial::itype)
        .def_property_readonly("izeta", &NumericalRadial::izeta)
        .def_property_readonly("l", &NumericalRadial::l)
        .def_property_readonly("nr", &NumericalRadial::nr)
        .def_property_readonly("nk", &NumericalRadial::nk)
        .def_property_readonly("rcut", &NumericalRadial::rcut)
        .def_property_readonly("kcut", &NumericalRadial::kcut)
        .def_property_readonly("rgrid",
                               [](NumericalRadial& self) {
                                   return numpy_from_ptr_copy(self.rgrid(), static_cast<size_t>(self.nr()));
                               })
        .def_property_readonly("kgrid",
                               [](NumericalRadial& self) {
                                   return numpy_from_ptr_copy(self.kgrid(), static_cast<size_t>(self.nk()));
                               })
        .def_property_readonly("rvalue",
                               [](NumericalRadial& self) {
                                   return numpy_from_ptr_copy(self.rvalue(), static_cast<size_t>(self.nr()));
                               })
        .def_property_readonly("kvalue",
                               [](NumericalRadial& self) {
                                   return numpy_from_ptr_copy(self.kvalue(), static_cast<size_t>(self.nk()));
                               })
        .def_property_readonly("pr", &NumericalRadial::pr)
        .def_property_readonly("pk", &NumericalRadial::pk)
        .def_property_readonly("is_fft_compliant", overload_cast_<>()(&NumericalRadial::is_fft_compliant, py::const_));
}
