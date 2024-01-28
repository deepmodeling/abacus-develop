#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "module_basis/module_nao/numerical_radial.h"

namespace py = pybind11;
using namespace pybind11::literals;
template <typename... Args>
using overload_cast_ = pybind11::detail::overload_cast_impl<Args...>;

PYBIND11_MODULE(_core, m)
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
                py::buffer_info grid_info = grid.request();
                py::buffer_info value_info = value.request();
                if (grid_info.ndim != 1 || value_info.ndim != 1)
                {
                    throw std::runtime_error("Input arrays must be 1-dimensional");
                }
                if (grid_info.shape[0] != value_info.shape[0])
                {
                    throw std::runtime_error("Grid and value arrays must have the same size");
                }

                self.build(l,
                           for_r_space,
                           grid_info.shape[0],
                           static_cast<double*>(grid_info.ptr),
                           static_cast<double*>(value_info.ptr),
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
                py::buffer_info grid_info = grid.request();
                if (grid_info.ndim != 1)
                {
                    throw std::runtime_error("Input array must be 1-dimensional");
                }

                self.set_grid(for_r_space, ngrid, static_cast<double*>(grid_info.ptr), mode);
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
             "enable+fft"_a = false)
        .def(
            "set_value",
            [](NumericalRadial& self, const bool for_r_space, py::array_t<double> value, const int p) {
                py::buffer_info value_info = value.request();
                if (value_info.ndim != 1)
                {
                    throw std::runtime_error("Input array must be 1-dimensional");
                }

                self.set_value(for_r_space, static_cast<double*>(value_info.ptr), p);
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
                py::buffer_info table_info = table.request();
                if (table_info.ndim != 1)
                {
                    throw std::runtime_error("Table array must be 1-dimensional");
                }

                self.radtab(op, ket, l, static_cast<double*>(table_info.ptr), nr_tab, rmax_tab, deriv);
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
                                   const double* rgrid = self.rgrid();
                                   return py::array_t<double>({self.nr()}, rgrid);
                               })
        .def_property_readonly("kgrid",
                               [](NumericalRadial& self) {
                                   const double* kgrid = self.kgrid();
                                   return py::array_t<double>({self.nk()}, kgrid);
                               })
        .def_property_readonly("rvalue",
                               [](NumericalRadial& self) {
                                   const double* rvalue = self.rvalue();
                                   return py::array_t<double>({self.nr()}, rvalue);
                               })
        .def_property_readonly("kvalue",
                               [](NumericalRadial& self) {
                                   const double* kvalue = self.kvalue();
                                   return py::array_t<double>({self.nk()}, kvalue);
                               })
        .def_property_readonly("pr", &NumericalRadial::pr)
        .def_property_readonly("pk", &NumericalRadial::pk)
        .def_property_readonly("is_fft_compliant", overload_cast_<>()(&NumericalRadial::is_fft_compliant, py::const_))
        // leave transformer for future
        .def_property_readonly("rgrid", overload_cast_<int>()(&NumericalRadial::rgrid, py::const_))
        .def_property_readonly("kgrid", overload_cast_<int>()(&NumericalRadial::kgrid, py::const_))
        .def_property_readonly("rvalue", overload_cast_<int>()(&NumericalRadial::rvalue, py::const_))
        .def_property_readonly("kvalue", overload_cast_<int>()(&NumericalRadial::kvalue, py::const_));
}
