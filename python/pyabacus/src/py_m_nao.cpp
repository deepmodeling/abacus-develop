#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "module_basis/module_nao/radial_collection.h"
#include "module_basis/module_nao/two_center_integrator.h"
#include "module_base/vector3.h"

namespace py = pybind11;
using namespace pybind11::literals;
template <typename... Args>
using overload_cast_ = pybind11::detail::overload_cast_impl<Args...>;

void bind_m_nao(py::module& m)
{
    // Create the submodule for Module NAO
    py::module m_nao = m.def_submodule("ModuleNAO");

    // Bind the RadialCollection class
    py::class_<RadialCollection>(m_nao, "RadialCollection")
        .def(py::init<>())
        .def("build", [](RadialCollection& self, int nfile, const py::list &file_list, char ftype){
             std::vector<std::string> files;
             files.reserve(nfile);
            for (auto file : file_list)
            {
                files.push_back(file.cast<std::string>());
            }
            self.build(nfile, files.data(), ftype);
        }, "nfile"_a, "file_list"_a, "ftype"_a = '\0')
        .def("set_transformer", &RadialCollection::set_transformer, "sbt"_a, "update"_a = 0)
        .def("set_uniform_grid", &RadialCollection::set_uniform_grid, "for_r_space"_a, "ngrid"_a, "cutoff"_a, "mode"_a = 'i', "enable_fft"_a = false)
        .def("set_grid", [](RadialCollection& self, const bool for_r_space, const int ngrid, py::array_t<double> grid, const char mode = 'i'){
            py::buffer_info grid_info = grid.request();
            if (grid_info.ndim != 1)
            {
                throw std::runtime_error("Input array must be 1-dimensional");
            }
            self.set_grid(for_r_space, ngrid, static_cast<double*>(grid_info.ptr), mode);
        }, "for_r_space"_a, "ngrid"_a, "grid"_a, "mode"_a = 'i')
        // Getters
        .def("symbol", &RadialCollection::symbol, "itype"_a)
        .def_property_readonly("ntype", &RadialCollection::ntype)
        .def("lmax", overload_cast_<const int>()(&RadialCollection::lmax, py::const_), "itype"_a)
        .def_property_readonly("lmax", overload_cast_<>()(&RadialCollection::lmax, py::const_))
        .def("rcut_max", overload_cast_<const int>()(&RadialCollection::rcut_max, py::const_), "itype"_a)
        .def_property_readonly("rcut_max", overload_cast_<>()(&RadialCollection::rcut_max, py::const_))
        .def("nzeta", &RadialCollection::nzeta, "itype"_a, "l"_a)
        .def("nzeta_max", overload_cast_<const int>()(&RadialCollection::nzeta_max, py::const_), "itype"_a)
        .def_property_readonly("nzeta_max", overload_cast_<>()(&RadialCollection::nzeta_max, py::const_))
        .def("nchi", overload_cast_<const int>()(&RadialCollection::nchi, py::const_), "itype"_a)
        .def_property_readonly("nchi", overload_cast_<>()(&RadialCollection::nchi, py::const_));
    //Bind the TwoCenterIntegrator class
    py::class_<TwoCenterIntegrator>(m_nao, "TwoCenterIntegrator")
        .def(py::init<>())
        .def("tabulate", &TwoCenterIntegrator::tabulate, "bra"_a, "ket"_a, "op"_a, "nr"_a, "cutoff"_a)
        .def("calculate", [](TwoCenterIntegrator& self, const int itype1, const int l1, const int izeta1, const int m1, const int itype2, const int l2, const int izeta2,  const int m2, py::array_t<double> pvR, bool cal_grad = false)
        {
            py::buffer_info pvR_info = pvR.request();
            if (pvR_info.size != 3)
            {
                throw std::runtime_error("Radial part must have 3 elements");
            }
            double* cvR = static_cast<double*>(pvR_info.ptr);
            ModuleBase::Vector3<double> vR(cvR[0], cvR[1], cvR[2]);
            double out[1] = {0.0};
            double* grad_out = nullptr;
            if (cal_grad)
            {
                grad_out = new double[3];
            }
            self.calculate(itype1, l1, izeta1, m1, itype2, l2, izeta2, m2, vR, out, grad_out);
            py::array_t<double> out_array({1}, out);
            if (cal_grad)
            {
                py::array_t<double> grad_out_array({3}, grad_out);
                return py::make_tuple(out_array, grad_out_array);
            }
            else
            {
                py::array_t<double> grad_out_array({0});
                return py::make_tuple(out_array, grad_out_array);
            }
            
        }, "itype1"_a, "l1"_a, "izeta1"_a, "m1"_a, "itype2"_a, "l2"_a, "izeta2"_a, "m2"_a, "pvR"_a, "cal_grad"_a = false)
        .def("snap", [](TwoCenterIntegrator& self, const int itype1, const int l1, const int izeta1, const int m1, const int itype2, py::array_t<double> pvR, const bool deriv){
            py::buffer_info pvR_info = pvR.request();
            if (pvR_info.size != 3)
            {
                throw std::runtime_error("Radial part must have 3 elements");
            }
            double* cvR = static_cast<double*>(pvR_info.ptr);
            ModuleBase::Vector3<double> vR(cvR[0], cvR[1], cvR[2]);
            // TODO: check deriv & out memory allocation
            std::vector<std::vector<double>> out;
            self.snap(itype1, l1, izeta1, m1, itype2, vR, deriv, out);
            return out;
        }, "itype1"_a, "l1"_a, "izeta1"_a, "m1"_a, "itype2"_a, "pvR"_a, "deriv"_a = false);
        
}